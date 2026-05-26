from __future__ import annotations

import importlib
import os
from pathlib import Path
from typing import Optional

import pandas as pd


def _load_scanpy_module():
    try:
        return importlib.import_module("scanpy")
    except ImportError as exc:
        raise ImportError(
            "scanpy is required for Xenium normalization. Install Cell-GPS[xenium] to use this option."
        ) from exc


def _load_pyxenium_io_module():
    try:
        return importlib.import_module("pyXenium.io")
    except ImportError as exc:
        raise ImportError(
            "pyXenium>=0.4.3 is required for Xenium loading. Install Cell-GPS[xenium] "
            "or install pyXenium in the active environment."
        ) from exc


def _normalize_if_requested(adata, normalize: bool):
    if not normalize:
        return adata
    sc = _load_scanpy_module()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    return adata


def _resolve_bundle_path(
    folder: Path,
    explicit_path: Optional[str | os.PathLike[str]],
    pattern: str,
) -> str | os.PathLike[str]:
    if explicit_path is not None:
        return explicit_path

    matches = sorted(folder.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"Could not find a file matching {pattern!r} under {folder}")
    if len(matches) > 1:
        names = ", ".join(str(path.name) for path in matches[:5])
        raise FileExistsError(f"Expected one file matching {pattern!r} under {folder}, found: {names}")
    return matches[0]


def _resolve_cluster_column(adata, preferred: str | None = None) -> str | None:
    candidates: list[str] = []
    if preferred:
        candidates.append(preferred)
    candidates.extend(["Cluster", "cluster", "Clusters"])
    candidates.extend(str(column) for column in adata.obs.columns if str(column).startswith("cluster__"))

    seen: set[str] = set()
    for column in candidates:
        if column in seen:
            continue
        seen.add(column)
        if column in adata.obs.columns:
            return column
    return None


def _ensure_spatial_coordinates(adata) -> None:
    if "spatial" in adata.obsm:
        return
    for x_col, y_col in (
        ("x_centroid", "y_centroid"),
        ("cell_x_centroid", "cell_y_centroid"),
        ("cell_centroid_x", "cell_centroid_y"),
        ("centroid_x", "centroid_y"),
        ("x", "y"),
    ):
        if x_col in adata.obs.columns and y_col in adata.obs.columns:
            adata.obsm["spatial"] = adata.obs[[x_col, y_col]].to_numpy(dtype=float)
            return
    raise ValueError("Xenium cell coordinates were not found in pyXenium output.")


def _coerce_cellgps_xenium_adata(
    adata,
    *,
    normalize: bool,
    cluster_col: str | None = None,
):
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)

    if "cell_id" not in adata.obs.columns:
        adata.obs["cell_id"] = adata.obs_names.astype(str)
    else:
        adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)

    resolved_cluster_col = _resolve_cluster_column(adata, preferred=cluster_col)
    if resolved_cluster_col is None:
        adata.obs["Cluster"] = pd.Series(pd.NA, index=adata.obs.index, dtype="string")
    else:
        adata.obs[resolved_cluster_col] = adata.obs[resolved_cluster_col].astype("string")
        adata.obs["Cluster"] = adata.obs[resolved_cluster_col].astype("string")

    _ensure_spatial_coordinates(adata)
    if adata.raw is None:
        adata.raw = adata.copy()
    return _normalize_if_requested(adata, normalize=normalize)


def load_xenium_table_bundle(
    folder: str | os.PathLike[str],
    *,
    cells_path: Optional[str | os.PathLike[str]] = None,
    cell_groups_path: Optional[str | os.PathLike[str]] = None,
    feature_matrix_path: Optional[str | os.PathLike[str]] = None,
    normalize: bool = False,
    cluster_col: str = "Clusters",
    cell_id_col: str = "Barcode",
    x_col: str = "x_centroid",
    y_col: str = "y_centroid",
):
    """
    Load a Xenium table bundle through ``pyXenium.io.read_xenium``.

    The returned object keeps the requested cluster labels in
    ``adata.obs[cluster_col]`` and mirrors them into ``adata.obs["Cluster"]`` for
    backward compatibility with the existing Cell-GPS API.
    """
    folder_path = Path(folder)
    px_io = _load_pyxenium_io_module()

    adata = px_io.read_xenium(
        str(folder_path),
        as_="anndata",
        prefer="h5",
        include_transcripts=False,
        include_boundaries=False,
        include_images=False,
        clusters_relpath=None,
        cluster_column_name=cluster_col,
        cells_path=_resolve_bundle_path(folder_path, cells_path, "cells.parquet"),
        cell_groups_path=_resolve_bundle_path(folder_path, cell_groups_path, "*_cell_groups.csv"),
        feature_matrix_path=_resolve_bundle_path(
            folder_path,
            feature_matrix_path,
            "cell_feature_matrix.h5",
        ),
        cell_id_col=cell_id_col,
        cluster_col=cluster_col,
        x_col=x_col,
        y_col=y_col,
    )
    return _coerce_cellgps_xenium_adata(
        adata,
        normalize=normalize,
        cluster_col=cluster_col,
    )


def load_xenium_data(folder: str, normalize: bool = True):
    """Load and preprocess a Xenium run through ``pyXenium.io.read_xenium``."""
    px_io = _load_pyxenium_io_module()
    adata = px_io.read_xenium(
        folder,
        as_="anndata",
        include_transcripts=False,
        include_boundaries=False,
        include_images=False,
        cluster_column_name="cluster",
    )
    return _coerce_cellgps_xenium_adata(adata, normalize=normalize, cluster_col="cluster")


_coerce_sfplot_xenium_adata = _coerce_cellgps_xenium_adata
