# Cell-GPS searcher_findee_score_gpu.py
"""GPU (PyTorch) backend for the searcher-findee nearest-neighbor distance kernel.

The whole Cell-GPS cophenetic / StructureMap pipeline is built on a single heavy
computation: the *directed inter-group average nearest-neighbor distance* matrix
(rows = "Searcher" clusters, columns = "Findee" clusters), in which every element is
the mean over the source-cluster cells of the distance to the nearest target-cluster
cell. On CPU this is done with ``sklearn.neighbors.NearestNeighbors`` (one 1-NN model
per target cluster). For large cell counts that nearest-neighbor step dominates the
runtime; everything downstream (``linkage`` / ``cophenet`` on the small k x k cluster
matrix) is negligible and stays on CPU.

This module moves only that nearest-neighbor step to the GPU via batched
``torch.cdist``. It is intentionally numerically equivalent to the CPU path:

* ``dtype="float64"`` (default) reproduces the scikit-learn result to ~1e-12 and is
  what the equivalence tests assert against. This is the contract: GPU == CPU.
* ``dtype="float32"`` is faster on the GPU at the cost of ~1e-5 relative deviation; use
  it only when that tolerance is acceptable for the downstream clustering.

The exact (non-matmul) Euclidean formula is used on purpose
(``compute_mode="donot_use_mm_for_euclid_dist"``); the matmul expansion
``sqrt(x^2 + y^2 - 2xy)`` loses precision for distant points and would break the
"identical to CPU" guarantee.
"""

from __future__ import annotations

from typing import Dict, Iterable, Optional, Union

import numpy as np
import pandas as pd
import torch

__all__ = [
    "compute_searcher_findee_distance_matrix_from_df_gpu",
    "nearest_cluster_distance_columns_gpu",
]

_DTYPES = {
    "float64": torch.float64,
    "double": torch.float64,
    "float32": torch.float32,
    "float": torch.float32,
}


def _resolve_dtype(dtype: Union[str, "torch.dtype"]) -> "torch.dtype":
    if isinstance(dtype, torch.dtype):
        return dtype
    try:
        return _DTYPES[str(dtype).lower()]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported dtype {dtype!r}; use 'float64' (exact) or 'float32' (fast)."
        ) from exc


def _nearest_distance_to_set(
    query: "torch.Tensor",
    ref: "torch.Tensor",
    max_elements: int,
) -> "torch.Tensor":
    """Distance from each ``query`` row to the nearest ``ref`` row (Euclidean).

    The (chunk x len(ref)) distance block is capped at ``max_elements`` entries so the
    computation stays within a memory budget for large inputs. Uses the exact Euclidean
    formula so results match scikit-learn's ``NearestNeighbors`` to floating-point
    precision.
    """
    n_query = query.shape[0]
    n_ref = ref.shape[0]
    out = torch.empty(n_query, device=query.device, dtype=query.dtype)
    chunk = max(1, int(max_elements // max(1, n_ref)))
    for start in range(0, n_query, chunk):
        end = min(n_query, start + chunk)
        block = torch.cdist(
            query[start:end],
            ref,
            p=2.0,
            compute_mode="donot_use_mm_for_euclid_dist",
        )
        out[start:end] = block.min(dim=1).values
    return out


def nearest_cluster_distance_columns_gpu(
    coords: np.ndarray,
    masks: Dict[object, np.ndarray],
    unique_clusters: Iterable,
    device: str = "cuda",
    dtype: Union[str, "torch.dtype"] = "float64",
    max_memory_gb: float = 4.0,
) -> np.ndarray:
    """Per-cell nearest-distance to each cluster as an ``(n_cells, n_clusters)`` array.

    Column ``j`` holds, for every cell, the distance to the nearest cell belonging to
    cluster ``list(unique_clusters)[j]``; a cluster with no cells yields an all-NaN
    column (mirroring the CPU behaviour). This is the shared GPU primitive used by every
    searcher-findee / cophenetic entry point.
    """
    torch_dtype = _resolve_dtype(dtype)
    coords_t = torch.as_tensor(
        np.ascontiguousarray(coords), device=device, dtype=torch_dtype
    )
    n_cells = coords_t.shape[0]
    unique_clusters = list(unique_clusters)
    n_clusters = len(unique_clusters)

    bytes_per = torch.empty(0, dtype=torch_dtype).element_size()
    max_elements = max(1, int((max_memory_gb * (1024 ** 3)) / bytes_per))

    result = np.full((n_cells, n_clusters), np.nan, dtype=np.float64)
    for j, c in enumerate(unique_clusters):
        mask_c = np.asarray(masks[c], dtype=bool)
        if not mask_c.any():
            continue  # cluster with no cells -> column stays NaN
        ref = coords_t[torch.as_tensor(mask_c, device=device)]
        mind = _nearest_distance_to_set(coords_t, ref, max_elements)
        result[:, j] = mind.detach().to("cpu", dtype=torch.float64).numpy()
    return result


def compute_searcher_findee_distance_matrix_from_df_gpu(
    df: pd.DataFrame,
    x_col: str = "x",
    y_col: str = "y",
    z_col: Optional[str] = None,
    celltype_col: str = "celltype",
    device: str = "cuda",
    dtype: Union[str, "torch.dtype"] = "float64",
    max_memory_gb: float = 4.0,
) -> pd.DataFrame:
    """GPU version of ``compute_searcher_findee_distance_matrix_from_df``.

    Returns the directed inter-cluster average nearest-neighbor distance matrix
    (n_clusters x n_clusters; rows = Searcher, columns = Findee). Numerically equivalent
    to the CPU implementation: with ``dtype="float64"`` the two agree to ~1e-12.

    Parameters mirror the CPU function, plus GPU controls:
        device: torch device, e.g. "cuda", "cuda:0", or "cpu" (used for testing).
        dtype: "float64" (exact, default) or "float32" (faster, ~1e-5 deviation).
        max_memory_gb: memory budget for each distance block; larger = fewer batches.
    """
    required_cols = {x_col, y_col, celltype_col}
    if z_col is not None:
        required_cols.add(z_col)
    if not required_cols.issubset(df.columns):
        raise ValueError(f"DataFrame must contain the following columns: {required_cols}")

    # Same category handling as the CPU kernel so row/column labels match exactly.
    clusters = df[celltype_col].astype("category")
    clusters = clusters.cat.remove_unused_categories()
    unique_clusters = clusters.cat.categories

    coord_cols = [x_col, y_col] + ([z_col] if z_col is not None else [])
    coords = df[coord_cols].values

    masks = {c: (clusters == c).values for c in unique_clusters}
    nearest = nearest_cluster_distance_columns_gpu(
        coords,
        masks,
        unique_clusters,
        device=device,
        dtype=dtype,
        max_memory_gb=max_memory_gb,
    )

    df_nearest_cluster_dist = pd.DataFrame(
        nearest, index=df.index, columns=unique_clusters
    )
    distance_matrix = df_nearest_cluster_dist.groupby(clusters, observed=False).mean()
    distance_matrix = distance_matrix.dropna(axis=1, how="all")
    return distance_matrix
