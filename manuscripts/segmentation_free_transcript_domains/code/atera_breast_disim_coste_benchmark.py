#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import psutil
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.sparse import csc_matrix, csr_matrix, issparse
from scipy.sparse.linalg import svds
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors


DEFAULT_INPUT_DIR = Path("/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast")
DEFAULT_OUTPUT_DIR = Path("/data/taobo.hu/atera_breast_disim_coste_benchmark")


@dataclass
class TimedResult:
    name: str
    seconds: float
    peak_rss_gb: float
    payload: Any


class PeakMemoryTracker:
    def __init__(self, interval_seconds: float = 0.2) -> None:
        self.interval_seconds = interval_seconds
        self.process = psutil.Process()
        self.peak_rss = 0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def __enter__(self) -> "PeakMemoryTracker":
        self._thread = threading.Thread(target=self._sample, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2)
        self._sample_once()

    def _sample(self) -> None:
        while not self._stop.is_set():
            self._sample_once()
            self._stop.wait(self.interval_seconds)

    def _sample_once(self) -> None:
        rss = 0
        try:
            rss += self.process.memory_info().rss
            for child in self.process.children(recursive=True):
                try:
                    rss += child.memory_info().rss
                except psutil.Error:
                    continue
        except psutil.Error:
            return
        self.peak_rss = max(self.peak_rss, int(rss))


def timed(name: str, func, *args, **kwargs) -> TimedResult:
    t0 = time.perf_counter()
    with PeakMemoryTracker() as tracker:
        payload = func(*args, **kwargs)
    return TimedResult(
        name=name,
        seconds=time.perf_counter() - t0,
        peak_rss_gb=tracker.peak_rss / (1024**3),
        payload=payload,
    )


def _decode(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def _git_revision(repo_root: Path | None) -> str | None:
    if repo_root is None:
        return None
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(repo_root),
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
    except Exception:
        return None
    return result.stdout.strip() or None


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def normalize_matrix(df: pd.DataFrame) -> pd.DataFrame:
    values = df.to_numpy(dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return pd.DataFrame(np.zeros_like(values), index=df.index, columns=df.columns)
    vmin = float(np.nanmin(values[finite]))
    vmax = float(np.nanmax(values[finite]))
    if vmax <= vmin:
        out = np.zeros_like(values, dtype=float)
    else:
        out = (values - vmin) / (vmax - vmin)
    out[~finite] = np.nan
    if out.shape[0] == out.shape[1]:
        np.fill_diagonal(out, 0.0)
    return pd.DataFrame(out, index=df.index, columns=df.columns)


def fill_numeric_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    m = matrix.astype(float).copy()
    if not m.isna().any().any():
        return m
    col_means = m.mean(axis=0)
    global_mean = float(np.nanmean(m.to_numpy(dtype=float)))
    if math.isnan(global_mean):
        global_mean = 0.0
    return m.fillna(col_means.fillna(global_mean))


def matrix_stats(matrix: pd.DataFrame) -> dict[str, Any]:
    values = matrix.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    out: dict[str, Any] = {
        "n_rows": int(matrix.shape[0]),
        "n_cols": int(matrix.shape[1]),
        "finite_fraction": float(finite.size / values.size) if values.size else 0.0,
    }
    if finite.size:
        out.update(
            {
                "min": float(np.min(finite)),
                "median": float(np.median(finite)),
                "mean": float(np.mean(finite)),
                "max": float(np.max(finite)),
            }
        )
    if matrix.shape[0] == matrix.shape[1]:
        arr = np.nan_to_num(values, nan=np.nanmean(finite) if finite.size else 0.0)
        out["asymmetry_rmse"] = float(np.sqrt(np.mean((arr - arr.T) ** 2)))
        out["diagonal_mean"] = float(np.mean(np.diag(arr)))
    return out


def load_xenium_metadata(input_dir: Path) -> dict[str, Any]:
    payload: dict[str, Any] = {"input_dir": str(input_dir)}
    experiment = input_dir / "experiment.xenium"
    if experiment.exists():
        data = json.loads(experiment.read_text(encoding="utf-8"))
        for key in (
            "run_name",
            "region_name",
            "panel_name",
            "panel_num_targets_predesigned",
            "panel_num_targets_custom",
            "num_cells",
            "transcripts_per_cell",
            "chemistry_version",
            "analysis_sw_version",
        ):
            payload[key] = data.get(key)
    metrics = input_dir / "metrics_summary.csv"
    if metrics.exists():
        frame = pd.read_csv(metrics)
        if not frame.empty:
            row = frame.iloc[0]
            for key in (
                "num_cells_detected",
                "median_genes_per_cell",
                "median_transcripts_per_cell",
                "total_high_quality_decoded_transcripts",
            ):
                if key in row.index and pd.notna(row[key]):
                    value = row[key]
                    payload[key] = int(value) if float(value).is_integer() else float(value)
    return payload


def load_cell_table(
    cells_path: Path,
    groups_path: Path,
    *,
    max_cells: int | None,
    random_state: int,
) -> pd.DataFrame:
    cells = pd.read_csv(
        cells_path,
        usecols=["cell_id", "x_centroid", "y_centroid", "transcript_counts", "total_counts"],
    )
    groups = pd.read_csv(groups_path, usecols=["cell_id", "group", "color"])
    merged = cells.merge(groups, on="cell_id", how="inner", validate="one_to_one")
    merged = merged.rename(columns={"x_centroid": "x", "y_centroid": "y", "group": "celltype"})
    merged = merged.dropna(subset=["x", "y", "celltype"]).reset_index(drop=True)
    if max_cells is not None and max_cells > 0 and len(merged) > max_cells:
        merged = merged.sample(n=max_cells, random_state=random_state).sort_index().reset_index(drop=True)
    return merged


def compute_celltype_searcher_findee_matrix(cells: pd.DataFrame) -> pd.DataFrame:
    labels = cells["celltype"].astype("category")
    categories = labels.cat.categories.to_list()
    coords = cells[["x", "y"]].to_numpy(dtype=float)
    nearest = np.empty((len(cells), len(categories)), dtype=np.float32)
    for idx, celltype in enumerate(categories):
        mask = (labels == celltype).to_numpy()
        model = NearestNeighbors(n_neighbors=1, algorithm="auto")
        model.fit(coords[mask])
        dist, _ = model.kneighbors(coords)
        nearest[:, idx] = dist[:, 0].astype(np.float32)
    nearest_frame = pd.DataFrame(nearest, columns=categories)
    distance = nearest_frame.groupby(labels, observed=True).mean()
    return distance.loc[categories, categories]


def read_feature_summaries(matrix_h5: Path) -> pd.DataFrame:
    import h5py

    with h5py.File(matrix_h5, "r") as handle:
        feature_names = _decode(handle["matrix/features/name"][:])
        feature_ids = _decode(handle["matrix/features/id"][:])
        feature_types = _decode(handle["matrix/features/feature_type"][:])
        indices = np.asarray(handle["matrix/indices"][:], dtype=np.int64)
        data = np.asarray(handle["matrix/data"][:], dtype=np.float64)
        n_features = int(handle["matrix/shape"][0])
    counts = np.bincount(indices, weights=data, minlength=n_features)
    nonzero_cells = np.bincount(indices, minlength=n_features)
    return pd.DataFrame(
        {
            "feature_index": np.arange(n_features, dtype=int),
            "feature_id": feature_ids,
            "gene": feature_names,
            "feature_type": feature_types,
            "total_counts": counts,
            "nonzero_entries": nonzero_cells,
        }
    )


def load_h5_selected_cell_gene_matrix(
    matrix_h5: Path,
    selected_feature_indices: np.ndarray,
) -> tuple[csr_matrix, list[str], list[str]]:
    import h5py

    with h5py.File(matrix_h5, "r") as handle:
        shape = tuple(int(x) for x in handle["matrix/shape"][:])
        data = np.asarray(handle["matrix/data"][:], dtype=np.float32)
        indices = np.asarray(handle["matrix/indices"][:], dtype=np.int64)
        indptr = np.asarray(handle["matrix/indptr"][:], dtype=np.int64)
        barcodes = _decode(handle["matrix/barcodes"][:])
        names = _decode(handle["matrix/features/name"][:])
    feature_by_cell = csc_matrix((data, indices, indptr), shape=shape)
    sub = feature_by_cell[selected_feature_indices, :].T.tocsr()
    selected_names = [names[int(i)] for i in selected_feature_indices]
    return sub, selected_names, barcodes


def build_spatial_bins(cells_path: Path, barcodes: list[str], grid_bins: int) -> tuple[np.ndarray, np.ndarray, int]:
    cells = pd.read_csv(cells_path, usecols=["cell_id", "x_centroid", "y_centroid"])
    cells = cells.set_index("cell_id")
    aligned = cells.reindex(barcodes)
    if aligned[["x_centroid", "y_centroid"]].isna().any().any():
        missing = int(aligned["x_centroid"].isna().sum())
        raise ValueError(f"{missing} matrix barcodes are missing from cells.csv.gz")
    coords = aligned[["x_centroid", "y_centroid"]].to_numpy(dtype=np.float32)
    x = coords[:, 0]
    y = coords[:, 1]
    x_edges = np.linspace(float(x.min()), float(x.max()), grid_bins + 1)
    y_edges = np.linspace(float(y.min()), float(y.max()), grid_bins + 1)
    bx = np.clip(np.searchsorted(x_edges, x, side="right") - 1, 0, grid_bins - 1)
    by = np.clip(np.searchsorted(y_edges, y, side="right") - 1, 0, grid_bins - 1)
    raw_bin = bx * grid_bins + by
    unique_bins, inverse = np.unique(raw_bin, return_inverse=True)
    n_bins = len(unique_bins)
    sums = np.zeros((n_bins, 2), dtype=np.float64)
    counts = np.bincount(inverse, minlength=n_bins).astype(np.float64)
    np.add.at(sums[:, 0], inverse, coords[:, 0])
    np.add.at(sums[:, 1], inverse, coords[:, 1])
    bin_coords = (sums / counts[:, None]).astype(np.float32)
    return inverse.astype(np.int64), bin_coords, n_bins


def aggregate_gene_expression_to_bins(
    cell_gene: csr_matrix,
    cell_to_bin: np.ndarray,
    n_bins: int,
) -> csr_matrix:
    n_cells = int(cell_gene.shape[0])
    bin_by_cell = csr_matrix(
        (np.ones(n_cells, dtype=np.float32), (cell_to_bin, np.arange(n_cells))),
        shape=(n_bins, n_cells),
    )
    return (bin_by_cell @ cell_gene).tocsr()


def compute_gene_bin_directed_matrix(
    matrix_h5: Path,
    cells_path: Path,
    *,
    max_genes: int,
    grid_bins: int,
    target_batch_size: int,
    use_gpu: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    features = read_feature_summaries(matrix_h5)
    gene_features = features.loc[features["feature_type"] == "Gene Expression"].copy()
    gene_features = gene_features.sort_values(["total_counts", "nonzero_entries"], ascending=False)
    panel_gene_count = int(len(gene_features))
    if max_genes > 0:
        gene_features = gene_features.head(max_genes)
    selected_indices = gene_features["feature_index"].to_numpy(dtype=np.int64)
    cell_gene, selected_names, barcodes = load_h5_selected_cell_gene_matrix(matrix_h5, selected_indices)
    cell_to_bin, bin_coords, n_bins = build_spatial_bins(cells_path, barcodes, grid_bins)
    bin_gene = aggregate_gene_expression_to_bins(cell_gene, cell_to_bin, n_bins)
    gene_sums = np.asarray(bin_gene.sum(axis=0)).ravel().astype(np.float32)
    keep = gene_sums > 0
    if not np.all(keep):
        bin_gene = bin_gene[:, keep]
        gene_sums = gene_sums[keep]
        gene_features = gene_features.loc[keep].copy()
        selected_names = [name for name, ok in zip(selected_names, keep) if bool(ok)]
    values = _compute_weighted_centroid_directed_matrix(
        bin_gene,
        bin_coords,
        gene_sums,
        target_batch_size=target_batch_size,
        use_gpu=use_gpu,
    )
    matrix = pd.DataFrame(values, index=selected_names, columns=selected_names)
    metadata = {
        "panel_gene_count": panel_gene_count,
        "benchmark_gene_count": int(len(selected_names)),
        "grid_bins_per_axis": int(grid_bins),
        "occupied_bins": int(n_bins),
        "target_batch_size": int(target_batch_size),
        "gene_selection": "top_total_counts",
        "feature_summary": {
            "selected_total_counts_min": float(gene_features["total_counts"].min()),
            "selected_total_counts_median": float(gene_features["total_counts"].median()),
            "selected_total_counts_max": float(gene_features["total_counts"].max()),
        },
    }
    selected = gene_features[["feature_index", "feature_id", "gene", "total_counts", "nonzero_entries"]].copy()
    return matrix, selected, metadata


def _compute_weighted_centroid_directed_matrix(
    bin_gene: csr_matrix,
    bin_coords: np.ndarray,
    gene_sums: np.ndarray,
    *,
    target_batch_size: int,
    use_gpu: bool,
) -> np.ndarray:
    weights = bin_gene.astype(np.float32)
    if issparse(weights):
        dense_weights = weights.toarray().astype(np.float32, copy=False)
    else:
        dense_weights = np.asarray(weights, dtype=np.float32)
    centers = (dense_weights.T @ bin_coords.astype(np.float32)) / gene_sums[:, None]
    n_genes = dense_weights.shape[1]

    if use_gpu:
        try:
            import torch

            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            if device.type == "cuda":
                _log(f"building directed matrix on CUDA for {n_genes} genes")
                w = torch.as_tensor(dense_weights, dtype=torch.float32, device=device)
                coords = torch.as_tensor(bin_coords, dtype=torch.float32, device=device)
                center_t = torch.as_tensor(centers, dtype=torch.float32, device=device)
                denom = torch.as_tensor(gene_sums, dtype=torch.float32, device=device).view(-1, 1)
                out = torch.empty((n_genes, n_genes), dtype=torch.float32, device=device)
                wt = w.transpose(0, 1).contiguous()
                for start in range(0, n_genes, target_batch_size):
                    stop = min(start + target_batch_size, n_genes)
                    dist = torch.cdist(coords, center_t[start:stop])
                    out[:, start:stop] = (wt @ dist) / denom
                    if start == 0 or stop == n_genes or (start // target_batch_size) % 20 == 0:
                        _log(f"directed matrix columns {start:,}-{stop:,} / {n_genes:,}")
                return out.cpu().numpy()
        except Exception as exc:
            print(f"[warn] GPU path failed; falling back to NumPy: {exc}", flush=True)

    _log(f"building directed matrix on CPU for {n_genes} genes")
    out = np.empty((n_genes, n_genes), dtype=np.float32)
    wt = dense_weights.T
    denom = gene_sums[:, None]
    for start in range(0, n_genes, target_batch_size):
        stop = min(start + target_batch_size, n_genes)
        diff = bin_coords[:, None, :] - centers[None, start:stop, :]
        dist = np.sqrt(np.sum(diff * diff, axis=2, dtype=np.float32))
        out[:, start:stop] = (wt @ dist) / denom
        if start == 0 or stop == n_genes or (start // target_batch_size) % 20 == 0:
            _log(f"directed matrix columns {start:,}-{stop:,} / {n_genes:,}")
    return out


def compute_coste_cophenetic(matrix: pd.DataFrame, *, method: str, metric: str) -> dict[str, Any]:
    m = fill_numeric_matrix(matrix)
    if m.shape[0] < 2 or m.shape[1] < 2:
        return {
            "row_cophenetic": pd.DataFrame(np.zeros((m.shape[0], m.shape[0])), index=m.index, columns=m.index),
            "col_cophenetic": pd.DataFrame(np.zeros((m.shape[1], m.shape[1])), index=m.columns, columns=m.columns),
            "row_cophenetic_corr": None,
            "col_cophenetic_corr": None,
        }
    row_dist = pdist(m.to_numpy(dtype=float), metric=metric)
    col_dist = pdist(m.T.to_numpy(dtype=float), metric=metric)
    row_linkage = linkage(row_dist, method=method)
    col_linkage = linkage(col_dist, method=method)
    row_corr, row_coph = cophenet(row_linkage, row_dist)
    col_corr, col_coph = cophenet(col_linkage, col_dist)
    row_coph_df = normalize_matrix(pd.DataFrame(squareform(row_coph), index=m.index, columns=m.index))
    col_coph_df = normalize_matrix(pd.DataFrame(squareform(col_coph), index=m.columns, columns=m.columns))
    return {
        "row_cophenetic": row_coph_df,
        "col_cophenetic": col_coph_df,
        "row_cophenetic_corr": float(row_corr),
        "col_cophenetic_corr": float(col_corr),
    }


def distance_to_affinity(
    matrix: pd.DataFrame,
    *,
    transform: str,
    zero_diagonal: bool,
) -> pd.DataFrame:
    m = fill_numeric_matrix(matrix)
    values = m.to_numpy(dtype=np.float64)
    finite = values[np.isfinite(values)]
    positive = finite[finite > 0]
    scale = float(np.median(positive)) if positive.size else 1.0
    if transform == "exp":
        affinity = np.exp(-values / max(scale, 1e-12))
    elif transform == "inverse":
        affinity = 1.0 / (1.0 + values / max(scale, 1e-12))
    elif transform == "linear":
        norm = normalize_matrix(m).to_numpy(dtype=np.float64)
        affinity = 1.0 - norm
    else:
        raise ValueError(f"Unsupported affinity transform: {transform}")
    affinity[~np.isfinite(affinity)] = 0.0
    if zero_diagonal and affinity.shape[0] == affinity.shape[1]:
        np.fill_diagonal(affinity, 0.0)
    return pd.DataFrame(affinity, index=m.index, columns=m.columns)


def _row_normalize(x: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(x, axis=1)
    norms[norms == 0] = 1.0
    return x / norms[:, None]


def run_disim(
    distance_matrix: pd.DataFrame,
    *,
    n_clusters: int,
    n_components: int,
    affinity_transform: str,
    zero_diagonal: bool,
    random_state: int,
    return_affinity: bool = True,
) -> dict[str, Any]:
    affinity = distance_to_affinity(
        distance_matrix,
        transform=affinity_transform,
        zero_diagonal=zero_diagonal,
    )
    a = affinity.to_numpy(dtype=np.float64)
    out_degree = a.sum(axis=1)
    in_degree = a.sum(axis=0)
    tau = float(a.sum() / max(a.shape[0], 1))
    regularized = a / np.sqrt((out_degree + tau)[:, None] * (in_degree + tau)[None, :])
    max_components = min(regularized.shape) - 1
    k = int(max(1, min(n_components, max_components)))
    if min(regularized.shape) <= 256 or k >= min(regularized.shape) - 1:
        u, singular_values, vt = np.linalg.svd(regularized, full_matrices=False)
        u = u[:, :k]
        singular_values = singular_values[:k]
        v = vt[:k, :].T
    else:
        u, singular_values, vt = svds(regularized, k=k, which="LM")
        order = np.argsort(singular_values)[::-1]
        singular_values = singular_values[order]
        u = u[:, order]
        v = vt[order, :].T
    row_raw_norms = np.linalg.norm(u, axis=1)
    col_raw_norms = np.linalg.norm(v, axis=1)
    row_embedding = _row_normalize(u)
    col_embedding = _row_normalize(v)
    row_clusters = min(n_clusters, row_embedding.shape[0])
    col_clusters = min(n_clusters, col_embedding.shape[0])
    row_labels = KMeans(n_clusters=row_clusters, n_init=50, random_state=random_state).fit_predict(row_embedding)
    col_labels = KMeans(n_clusters=col_clusters, n_init=50, random_state=random_state).fit_predict(col_embedding)
    result: dict[str, Any] = {
        "row_embedding": pd.DataFrame(row_embedding, index=distance_matrix.index),
        "col_embedding": pd.DataFrame(col_embedding, index=distance_matrix.columns),
        "row_raw_norm": pd.Series(row_raw_norms, index=distance_matrix.index, name="row_raw_norm"),
        "col_raw_norm": pd.Series(col_raw_norms, index=distance_matrix.columns, name="col_raw_norm"),
        "assignments": pd.DataFrame(
            {
                "node": list(distance_matrix.index) + list(distance_matrix.columns),
                "role": ["row"] * len(distance_matrix.index) + ["column"] * len(distance_matrix.columns),
                "cluster": np.concatenate([row_labels, col_labels]).astype(int),
            }
        ),
        "singular_values": singular_values.astype(float).tolist(),
        "tau": tau,
        "n_components": k,
        "row_cluster_count": int(row_clusters),
        "col_cluster_count": int(col_clusters),
    }
    if return_affinity:
        result["affinity"] = affinity
    if distance_matrix.shape[0] == distance_matrix.shape[1] and list(distance_matrix.index) == list(distance_matrix.columns):
        result["row_column_adjusted_rand"] = float(adjusted_rand_score(row_labels, col_labels))
        result["embedding_asymmetry_rmse"] = float(np.sqrt(np.mean((row_embedding - col_embedding) ** 2)))
    return result


def compute_disim_importance(disim: dict[str, Any], feature_table: pd.DataFrame | None = None) -> pd.DataFrame:
    row_embedding = disim["row_embedding"]
    col_embedding = disim["col_embedding"]
    col_names = set(col_embedding.index)
    shared = [name for name in row_embedding.index if name in col_names]
    row = row_embedding.loc[shared].to_numpy(dtype=float)
    col = col_embedding.loc[shared].to_numpy(dtype=float)
    row_raw = disim["row_raw_norm"].reindex(shared).to_numpy(dtype=float)
    col_raw = disim["col_raw_norm"].reindex(shared).to_numpy(dtype=float)
    asymmetry = np.linalg.norm(row - col, axis=1)
    leverage = np.maximum(row_raw, col_raw)

    table = pd.DataFrame(
        {
            "gene": shared,
            "disim_asymmetry": asymmetry,
            "row_raw_norm": row_raw,
            "col_raw_norm": col_raw,
            "spectral_leverage": leverage,
        }
    )
    row_assign = disim["assignments"].loc[disim["assignments"]["role"] == "row", ["node", "cluster"]]
    col_assign = disim["assignments"].loc[disim["assignments"]["role"] == "column", ["node", "cluster"]]
    table = table.merge(row_assign.rename(columns={"node": "gene", "cluster": "row_cluster"}), on="gene", how="left")
    table = table.merge(col_assign.rename(columns={"node": "gene", "cluster": "col_cluster"}), on="gene", how="left")
    table["asymmetry_percentile"] = table["disim_asymmetry"].rank(pct=True)
    table["leverage_percentile"] = table["spectral_leverage"].rank(pct=True)
    table["importance_score"] = 0.7 * table["asymmetry_percentile"] + 0.3 * table["leverage_percentile"]
    if feature_table is not None and "gene" in feature_table.columns:
        keep_cols = [c for c in ("gene", "feature_id", "total_counts", "nonzero_entries") if c in feature_table.columns]
        table = table.merge(feature_table[keep_cols], on="gene", how="left")
    return table.sort_values("importance_score", ascending=False).reset_index(drop=True)


def embedding_distance_matrix(embedding: pd.DataFrame) -> pd.DataFrame:
    if len(embedding) < 2:
        return pd.DataFrame(np.zeros((len(embedding), len(embedding))), index=embedding.index, columns=embedding.index)
    d = squareform(pdist(embedding.to_numpy(dtype=float), metric="euclidean"))
    return normalize_matrix(pd.DataFrame(d, index=embedding.index, columns=embedding.index))


def compare_coste_disim(coste: dict[str, Any], disim: dict[str, Any]) -> dict[str, Any]:
    row_emb_dist = embedding_distance_matrix(disim["row_embedding"])
    col_emb_dist = embedding_distance_matrix(disim["col_embedding"])
    return {
        "row_cophenetic_vs_disim_embedding_pearson": _upper_triangle_corr(coste["row_cophenetic"], row_emb_dist),
        "col_cophenetic_vs_disim_embedding_pearson": _upper_triangle_corr(coste["col_cophenetic"], col_emb_dist),
    }


def _upper_triangle_corr(a: pd.DataFrame, b: pd.DataFrame) -> float | None:
    if a.shape != b.shape or a.shape[0] < 3:
        return None
    av = a.to_numpy(dtype=float)
    bv = b.to_numpy(dtype=float)
    idx = np.triu_indices_from(av, k=1)
    x = av[idx]
    y = bv[idx]
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return None
    return float(np.corrcoef(x[ok], y[ok])[0, 1])


def build_matrix(args: argparse.Namespace) -> tuple[pd.DataFrame, dict[str, Any], dict[str, pd.DataFrame]]:
    input_dir = args.input_dir
    if args.mode == "celltype":
        cells = load_cell_table(
            input_dir / "cells.csv.gz",
            input_dir / "WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv",
            max_cells=args.max_cells,
            random_state=args.random_state,
        )
        matrix = compute_celltype_searcher_findee_matrix(cells)
        metadata = {
            "mode": args.mode,
            "analysis_cell_count": int(len(cells)),
            "celltype_count": int(matrix.shape[0]),
        }
        extra = {
            "celltype_counts": cells["celltype"].value_counts().rename_axis("celltype").reset_index(name="n_cells")
        }
        return matrix, metadata, extra
    if args.mode == "gene-bin":
        matrix, selected, metadata = compute_gene_bin_directed_matrix(
            input_dir / "cell_feature_matrix.h5",
            input_dir / "cells.csv.gz",
            max_genes=args.max_genes,
            grid_bins=args.grid_bins,
            target_batch_size=args.target_batch_size,
            use_gpu=not args.cpu_only,
        )
        metadata["mode"] = args.mode
        return matrix, metadata, {"selected_genes": selected}
    raise ValueError(f"Unsupported mode: {args.mode}")


def run(args: argparse.Namespace) -> None:
    if args.workflow == "disim-select-coste":
        run_disim_select_coste(args)
        return

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    benchmark_rows: list[dict[str, Any]] = []

    matrix_result = timed("matrix_build", build_matrix, args)
    matrix, matrix_metadata, extra_tables = matrix_result.payload
    benchmark_rows.append(_timed_row(matrix_result, method="shared", n_items=matrix.shape[0]))

    write_matrix_outputs(matrix, output_dir, args.matrix_output, stem="directed_searcher_findee_matrix")
    for name, table in extra_tables.items():
        table.to_csv(output_dir / f"{name}.csv", index=False)

    coste_result = timed(
        "coste_cophenetic",
        compute_coste_cophenetic,
        matrix,
        method=args.linkage_method,
        metric=args.cophenetic_metric,
    )
    coste = coste_result.payload
    benchmark_rows.append(_timed_row(coste_result, method="COSTE", n_items=matrix.shape[0]))
    coste["row_cophenetic"].to_csv(output_dir / "coste_row_cophenetic.csv")
    coste["col_cophenetic"].to_csv(output_dir / "coste_col_cophenetic.csv")

    disim_result = timed(
        "disim",
        run_disim,
        matrix,
        n_clusters=args.n_clusters,
        n_components=args.n_components,
        affinity_transform=args.affinity_transform,
        zero_diagonal=not args.keep_diagonal,
        random_state=args.random_state,
        return_affinity=True,
    )
    disim = disim_result.payload
    benchmark_rows.append(_timed_row(disim_result, method="PNAS_di_sim", n_items=matrix.shape[0]))
    disim["affinity"].to_csv(output_dir / "disim_affinity_matrix.csv")
    disim["row_embedding"].to_csv(output_dir / "disim_row_embedding.csv")
    disim["col_embedding"].to_csv(output_dir / "disim_col_embedding.csv")
    disim["assignments"].to_csv(output_dir / "disim_assignments.csv", index=False)

    comparison = compare_coste_disim(coste, disim)
    summary = {
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "input_metadata": load_xenium_metadata(args.input_dir),
        "matrix_metadata": matrix_metadata,
        "matrix_stats": matrix_stats(matrix),
        "coste": {
            "row_cophenetic_corr": coste["row_cophenetic_corr"],
            "col_cophenetic_corr": coste["col_cophenetic_corr"],
            "linkage_method": args.linkage_method,
            "cophenetic_metric": args.cophenetic_metric,
        },
        "disim": {
            key: value
            for key, value in disim.items()
            if key
            not in {
                "affinity",
                "row_embedding",
                "col_embedding",
                "assignments",
            }
        },
        "comparison": comparison,
        "pnas_di_sim_reference": {
            "doi": "10.1073/pnas.1525793113",
            "method": "regularized directed graph Laplacian, top singular vectors, row-normalized left/right embeddings, k-means co-clustering",
        },
        "repo_revision": _git_revision(args.repo_root),
    }
    _write_json(output_dir / "benchmark_manifest.json", summary)
    pd.DataFrame(benchmark_rows).to_csv(output_dir / "method_summary.csv", index=False)
    pd.DataFrame([comparison]).to_csv(output_dir / "coste_disim_comparison.csv", index=False)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


def write_matrix_outputs(matrix: pd.DataFrame, output_dir: Path, matrix_output: str, *, stem: str) -> None:
    if matrix_output in {"csv", "both"}:
        matrix.to_csv(output_dir / f"{stem}.csv")
    if matrix_output in {"npy", "both"}:
        np.save(output_dir / f"{stem}.float32.npy", matrix.to_numpy(dtype=np.float32))
        pd.Series(matrix.index, name="row_label").to_csv(output_dir / f"{stem}.rows.csv", index=False)
        pd.Series(matrix.columns, name="column_label").to_csv(output_dir / f"{stem}.columns.csv", index=False)


def run_disim_select_coste(args: argparse.Namespace) -> None:
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    benchmark_rows: list[dict[str, Any]] = []

    if args.mode != "gene-bin":
        raise ValueError("disim-select-coste workflow currently expects --mode gene-bin")
    if args.max_genes != 0:
        _log(f"workflow uses the requested top {args.max_genes} genes; pass --max-genes 0 for the full panel")

    _log("building directed gene matrix")
    matrix_result = timed("matrix_build", build_matrix, args)
    matrix, matrix_metadata, extra_tables = matrix_result.payload
    benchmark_rows.append(_timed_row(matrix_result, method="shared", n_items=matrix.shape[0]))
    write_matrix_outputs(matrix, output_dir, args.matrix_output, stem="directed_searcher_findee_matrix")
    for name, table in extra_tables.items():
        table.to_csv(output_dir / f"{name}.csv", index=False)

    _log("running full di-sim")
    disim_result = timed(
        "full_disim",
        run_disim,
        matrix,
        n_clusters=args.n_clusters,
        n_components=args.n_components,
        affinity_transform=args.affinity_transform,
        zero_diagonal=not args.keep_diagonal,
        random_state=args.random_state,
        return_affinity=False,
    )
    disim = disim_result.payload
    benchmark_rows.append(_timed_row(disim_result, method="PNAS_di_sim_full", n_items=matrix.shape[0]))
    disim["row_embedding"].to_csv(output_dir / "full_disim_row_embedding.csv")
    disim["col_embedding"].to_csv(output_dir / "full_disim_col_embedding.csv")
    disim["assignments"].to_csv(output_dir / "full_disim_assignments.csv", index=False)

    feature_table = extra_tables.get("selected_genes")
    importance = compute_disim_importance(disim, feature_table=feature_table)
    importance.to_csv(output_dir / "full_disim_gene_importance.csv", index=False)

    eligible_importance = importance
    if args.min_selected_total_counts > 0 and "total_counts" in eligible_importance.columns:
        eligible_importance = eligible_importance.loc[
            pd.to_numeric(eligible_importance["total_counts"], errors="coerce").fillna(0)
            >= args.min_selected_total_counts
        ].copy()
    selected_count = min(args.selected_coste_genes, len(eligible_importance))
    selected_genes = eligible_importance.head(selected_count)["gene"].astype(str).tolist()
    selected_importance = eligible_importance.head(selected_count).copy()
    selected_importance.to_csv(output_dir / f"selected_top{selected_count}_genes_for_coste.csv", index=False)
    selected_matrix = matrix.loc[selected_genes, selected_genes].copy()
    write_matrix_outputs(selected_matrix, output_dir, "csv", stem=f"selected_top{selected_count}_directed_matrix")

    _log(f"running COSTE cophenetic on {selected_count} di-sim-selected genes")
    coste_result = timed(
        "selected_coste_cophenetic",
        compute_coste_cophenetic,
        selected_matrix,
        method=args.linkage_method,
        metric=args.cophenetic_metric,
    )
    coste = coste_result.payload
    benchmark_rows.append(_timed_row(coste_result, method="COSTE_selected", n_items=selected_matrix.shape[0]))
    coste["row_cophenetic"].to_csv(output_dir / f"selected_top{selected_count}_coste_row_cophenetic.csv")
    coste["col_cophenetic"].to_csv(output_dir / f"selected_top{selected_count}_coste_col_cophenetic.csv")

    selected_disim = {
        **disim,
        "row_embedding": disim["row_embedding"].loc[selected_genes],
        "col_embedding": disim["col_embedding"].loc[selected_genes],
    }
    comparison = compare_coste_disim(coste, selected_disim)
    pd.DataFrame([comparison]).to_csv(output_dir / f"selected_top{selected_count}_coste_disim_comparison.csv", index=False)

    summary = {
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "workflow": args.workflow,
        "input_metadata": load_xenium_metadata(args.input_dir),
        "matrix_metadata": matrix_metadata,
        "matrix_stats": matrix_stats(matrix),
        "selected_coste_gene_count": int(selected_count),
        "selected_coste_selection_rule": "top importance_score = 0.7 * di-sim row/column asymmetry percentile + 0.3 * spectral leverage percentile",
        "selected_coste_min_total_counts": float(args.min_selected_total_counts),
        "selected_coste_eligible_gene_count": int(len(eligible_importance)),
        "selected_matrix_stats": matrix_stats(selected_matrix),
        "coste_selected": {
            "row_cophenetic_corr": coste["row_cophenetic_corr"],
            "col_cophenetic_corr": coste["col_cophenetic_corr"],
            "linkage_method": args.linkage_method,
            "cophenetic_metric": args.cophenetic_metric,
        },
        "disim_full": {
            key: value
            for key, value in disim.items()
            if key
            not in {
                "affinity",
                "row_embedding",
                "col_embedding",
                "row_raw_norm",
                "col_raw_norm",
                "assignments",
            }
        },
        "selected_comparison": comparison,
        "pnas_di_sim_reference": {
            "doi": "10.1073/pnas.1525793113",
            "method": "regularized directed graph Laplacian, top singular vectors, row-normalized left/right embeddings, k-means co-clustering",
        },
        "repo_revision": _git_revision(args.repo_root),
    }
    _write_json(output_dir / "benchmark_manifest.json", summary)
    pd.DataFrame(benchmark_rows).to_csv(output_dir / "method_summary.csv", index=False)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


def _timed_row(result: TimedResult, *, method: str, n_items: int) -> dict[str, Any]:
    return {
        "step": result.name,
        "method": method,
        "n_items": int(n_items),
        "seconds": float(result.seconds),
        "peak_rss_gb": float(result.peak_rss_gb),
    }


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark PNAS di-sim and COSTE on the Atera WTA breast dataset."
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--workflow",
        choices=["benchmark", "disim-select-coste"],
        default="benchmark",
        help="benchmark runs COSTE and di-sim on the whole matrix; disim-select-coste runs full di-sim then COSTE on important genes.",
    )
    parser.add_argument("--mode", choices=["celltype", "gene-bin"], default="gene-bin")
    parser.add_argument("--max-cells", type=int, default=None, help="Optional cell subsample for celltype mode.")
    parser.add_argument("--max-genes", type=int, default=1024, help="Top expressed genes for gene-bin mode; 0 means all genes.")
    parser.add_argument("--grid-bins", type=positive_int, default=96, help="Spatial bins per axis for gene-bin mode.")
    parser.add_argument("--target-batch-size", type=positive_int, default=128)
    parser.add_argument("--cpu-only", action="store_true", help="Disable torch/cuda for gene-bin matrix construction.")
    parser.add_argument("--matrix-output", choices=["csv", "npy", "both", "none"], default="csv")
    parser.add_argument("--selected-coste-genes", type=positive_int, default=1024)
    parser.add_argument(
        "--min-selected-total-counts",
        type=float,
        default=0.0,
        help="Minimum total Atera counts for genes eligible for selected-gene COSTE.",
    )
    parser.add_argument("--linkage-method", default="average")
    parser.add_argument("--cophenetic-metric", default="euclidean")
    parser.add_argument("--n-components", type=positive_int, default=8)
    parser.add_argument("--n-clusters", type=positive_int, default=8)
    parser.add_argument("--affinity-transform", choices=["exp", "inverse", "linear"], default="exp")
    parser.add_argument("--keep-diagonal", action="store_true")
    parser.add_argument("--random-state", type=int, default=7)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
