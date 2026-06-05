#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
import zarr
from scipy.cluster.hierarchy import cophenet, fcluster, linkage
from scipy.ndimage import gaussian_filter
from scipy.sparse.linalg import svds
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler


DEFAULT_TRANSCRIPTS_ZARR = Path(
    "/data/taobo.hu/pyxenium_lr_benchmark_2026-04/data/source_cache/breast/"
    "WTA_Preview_FFPE_Breast_Cancer_outs/spatialdata.zarr/points/transcripts"
)
DEFAULT_OUTPUT_DIR = Path("/data/taobo.hu/atera_breast_cellfree_transcript_pipeline")
DEFAULT_CELL_BASED_DISIM_DIR = Path(
    "/data/taobo.hu/atera_breast_disim_coste_benchmark/full_disim_18028_select1024_min5000_grid96_k16"
)
DEFAULT_CELL_BASED_HISTOSEG_DIR = Path(
    "/data/taobo.hu/atera_breast_coste1024_histoseg/coste1024_modules16_domains12"
)


@dataclass
class TimedResult:
    step: str
    seconds: float
    peak_rss_gb: float
    payload: Any


class PeakMemoryTracker:
    def __init__(self, interval_seconds: float = 0.5) -> None:
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


def timed(step: str, func, *args, **kwargs) -> TimedResult:
    t0 = time.perf_counter()
    with PeakMemoryTracker() as tracker:
        payload = func(*args, **kwargs)
    return TimedResult(step=step, seconds=time.perf_counter() - t0, peak_rss_gb=tracker.peak_rss / (1024**3), payload=payload)


def _log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def resize_int_array(values: np.ndarray, needed: int) -> np.ndarray:
    if needed < len(values):
        return values
    new_len = max(needed + 1, len(values) * 2)
    out = np.zeros(new_len, dtype=values.dtype)
    out[: len(values)] = values
    return out


def transcript_range(n_rows: int, chunk_rows: int, max_transcripts: int) -> list[tuple[int, int]]:
    limit = n_rows if max_transcripts <= 0 else min(n_rows, max_transcripts)
    return [(start, min(start + chunk_rows, limit)) for start in range(0, limit, chunk_rows)]


def scan_transcript_metadata(
    transcripts_zarr: Path,
    *,
    qv_min: float,
    chunk_rows: int,
    max_transcripts: int,
) -> dict[str, Any]:
    group = zarr.open_group(str(transcripts_zarr), mode="r")
    n_rows = int(group.attrs.get("n_rows", group["x"].shape[0]))
    chunks = transcript_range(n_rows, chunk_rows, max_transcripts)
    counts_by_identity = np.zeros(65536, dtype=np.int64)
    identity_to_name: dict[int, str] = {}
    xmin = ymin = math.inf
    xmax = ymax = -math.inf
    total_kept = 0
    for chunk_idx, (start, stop) in enumerate(chunks, start=1):
        valid = np.asarray(group["valid"][start:stop], dtype=bool)
        qv = np.asarray(group["quality_score"][start:stop], dtype=np.float32)
        mask = valid & (qv >= qv_min)
        if not mask.any():
            continue
        ids = np.asarray(group["gene_identity"][start:stop], dtype=np.int64)[mask]
        max_id = int(ids.max(initial=0))
        counts_by_identity = resize_int_array(counts_by_identity, max_id)
        counts_by_identity[: max_id + 1] += np.bincount(ids, minlength=max_id + 1)
        total_kept += int(len(ids))

        x = np.asarray(group["x"][start:stop], dtype=np.float32)[mask]
        y = np.asarray(group["y"][start:stop], dtype=np.float32)[mask]
        xmin = min(xmin, float(x.min()))
        xmax = max(xmax, float(x.max()))
        ymin = min(ymin, float(y.min()))
        ymax = max(ymax, float(y.max()))

        missing_ids = [identity for identity in np.unique(ids) if int(identity) not in identity_to_name]
        if missing_ids:
            names = np.asarray(group["gene_name"][start:stop])[mask]
            unique_ids, first_idx = np.unique(ids, return_index=True)
            for identity, idx in zip(unique_ids, first_idx):
                identity_int = int(identity)
                if identity_int not in identity_to_name:
                    identity_to_name[identity_int] = str(names[int(idx)])

        if chunk_idx == 1 or chunk_idx % 25 == 0 or chunk_idx == len(chunks):
            _log(f"metadata scan chunk {chunk_idx:,}/{len(chunks):,}: kept {total_kept:,} transcripts")

    identities_all = np.flatnonzero(counts_by_identity > 0).astype(np.int64)
    gene_names_all = [identity_to_name.get(int(identity), f"identity_{int(identity)}") for identity in identities_all]
    gene_mask = np.array([is_gene_expression_name(name) for name in gene_names_all], dtype=bool)
    identities = identities_all[gene_mask]
    gene_names = [name for name, keep in zip(gene_names_all, gene_mask) if bool(keep)]
    gene_counts = counts_by_identity[identities].astype(np.int64)
    return {
        "n_rows": n_rows,
        "n_rows_scanned": chunks[-1][1] if chunks else 0,
        "qv_min": qv_min,
        "n_transcripts_kept": int(total_kept),
        "n_gene_expression_transcripts_kept": int(gene_counts.sum()),
        "n_non_gene_identities_seen": int((~gene_mask).sum()),
        "identities": identities,
        "gene_names": gene_names,
        "counts_by_identity": gene_counts,
        "bounds": {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax},
    }


def is_gene_expression_name(name: str) -> bool:
    value = str(name).strip()
    if not value:
        return False
    banned_prefixes = (
        "NegControl",
        "UnassignedCodeword",
        "DeprecatedCodeword",
        "GenomicControl",
        "Intergenic_Region",
        "Blank",
    )
    if value.startswith(banned_prefixes):
        return False
    if "Codeword" in value and value.startswith(("Neg", "Unassigned", "Deprecated")):
        return False
    return True


def aggregate_transcripts_to_grid(
    transcripts_zarr: Path,
    *,
    identities: np.ndarray,
    bounds: dict[str, float],
    qv_min: float,
    grid_bins: int,
    chunk_rows: int,
    max_transcripts: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    group = zarr.open_group(str(transcripts_zarr), mode="r")
    n_rows = int(group.attrs.get("n_rows", group["x"].shape[0]))
    chunks = transcript_range(n_rows, chunk_rows, max_transcripts)
    n_genes = len(identities)
    identity_to_col = np.full(int(identities.max()) + 1, -1, dtype=np.int32)
    identity_to_col[identities] = np.arange(n_genes, dtype=np.int32)
    counts = np.zeros((grid_bins * grid_bins, n_genes), dtype=np.float32)
    x_edges = np.linspace(bounds["xmin"], bounds["xmax"], grid_bins + 1)
    y_edges = np.linspace(bounds["ymin"], bounds["ymax"], grid_bins + 1)
    total_kept = 0
    for chunk_idx, (start, stop) in enumerate(chunks, start=1):
        valid = np.asarray(group["valid"][start:stop], dtype=bool)
        qv = np.asarray(group["quality_score"][start:stop], dtype=np.float32)
        mask = valid & (qv >= qv_min)
        if not mask.any():
            continue
        raw_ids = np.asarray(group["gene_identity"][start:stop], dtype=np.int64)[mask]
        in_map = raw_ids < len(identity_to_col)
        raw_ids = raw_ids[in_map]
        cols = identity_to_col[raw_ids]
        valid_cols = cols >= 0
        if not valid_cols.any():
            continue
        cols = cols[valid_cols]
        x = np.asarray(group["x"][start:stop], dtype=np.float32)[mask][in_map][valid_cols]
        y = np.asarray(group["y"][start:stop], dtype=np.float32)[mask][in_map][valid_cols]
        xi = np.clip(np.searchsorted(x_edges, x, side="right") - 1, 0, grid_bins - 1)
        yi = np.clip(np.searchsorted(y_edges, y, side="right") - 1, 0, grid_bins - 1)
        bins = yi.astype(np.int64) * grid_bins + xi.astype(np.int64)
        flat = bins * n_genes + cols.astype(np.int64)
        unique_flat, flat_counts = np.unique(flat, return_counts=True)
        counts.ravel()[unique_flat] += flat_counts.astype(np.float32)
        total_kept += int(len(cols))
        if chunk_idx == 1 or chunk_idx % 25 == 0 or chunk_idx == len(chunks):
            _log(f"aggregate chunk {chunk_idx:,}/{len(chunks):,}: aggregated {total_kept:,} transcripts")

    xs = (x_edges[:-1] + x_edges[1:]) / 2.0
    ys = (y_edges[:-1] + y_edges[1:]) / 2.0
    xx, yy = np.meshgrid(xs, ys)
    coords = np.column_stack([xx.ravel(), yy.ravel()]).astype(np.float32)
    bin_totals = counts.sum(axis=1)
    occupied = bin_totals > 0
    return counts[occupied], coords[occupied], occupied.reshape(grid_bins, grid_bins)


def compute_directed_matrix(
    bin_gene: np.ndarray,
    bin_coords: np.ndarray,
    *,
    target_batch_size: int,
    use_gpu: bool,
) -> np.ndarray:
    gene_sums = bin_gene.sum(axis=0).astype(np.float32)
    keep = gene_sums > 0
    if not np.all(keep):
        raise ValueError("bin_gene contains zero-count genes; filter before matrix computation.")
    centers = (bin_gene.T @ bin_coords.astype(np.float32)) / gene_sums[:, None]
    n_genes = bin_gene.shape[1]

    if use_gpu:
        try:
            import torch

            if torch.cuda.is_available():
                device = torch.device("cuda")
                _log(f"building transcript-only directed matrix on CUDA for {n_genes:,} genes")
                w = torch.as_tensor(bin_gene, dtype=torch.float32, device=device)
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
            _log(f"GPU matrix path failed; falling back to CPU: {exc}")

    _log(f"building transcript-only directed matrix on CPU for {n_genes:,} genes")
    out = np.empty((n_genes, n_genes), dtype=np.float32)
    wt = bin_gene.T
    denom = gene_sums[:, None]
    for start in range(0, n_genes, target_batch_size):
        stop = min(start + target_batch_size, n_genes)
        diff = bin_coords[:, None, :] - centers[None, start:stop, :]
        dist = np.sqrt(np.sum(diff * diff, axis=2, dtype=np.float32))
        out[:, start:stop] = (wt @ dist) / denom
        if start == 0 or stop == n_genes or (start // target_batch_size) % 20 == 0:
            _log(f"directed matrix columns {start:,}-{stop:,} / {n_genes:,}")
    return out


def run_disim_array(
    distance_matrix: np.ndarray,
    gene_names: list[str],
    *,
    n_components: int,
    n_clusters: int,
    random_state: int,
) -> dict[str, Any]:
    values = distance_matrix.astype(np.float32, copy=False)
    positive = values[values > 0]
    scale = float(np.median(positive)) if positive.size else 1.0
    affinity = np.exp(-values / max(scale, 1e-12), dtype=np.float32)
    np.fill_diagonal(affinity, 0.0)
    out_degree = affinity.sum(axis=1)
    in_degree = affinity.sum(axis=0)
    tau = float(affinity.sum() / max(affinity.shape[0], 1))
    regularized = affinity / np.sqrt((out_degree + tau)[:, None] * (in_degree + tau)[None, :])
    k = int(max(1, min(n_components, min(regularized.shape) - 1)))
    u, singular_values, vt = svds(regularized, k=k, which="LM")
    order = np.argsort(singular_values)[::-1]
    singular_values = singular_values[order]
    u = u[:, order]
    v = vt[order, :].T
    row_raw_norm = np.linalg.norm(u, axis=1)
    col_raw_norm = np.linalg.norm(v, axis=1)
    row_embedding = row_normalize(u)
    col_embedding = row_normalize(v)
    row_labels = KMeans(n_clusters=min(n_clusters, len(gene_names)), n_init=50, random_state=random_state).fit_predict(row_embedding)
    col_labels = KMeans(n_clusters=min(n_clusters, len(gene_names)), n_init=50, random_state=random_state).fit_predict(col_embedding)
    importance = compute_importance(gene_names, row_embedding, col_embedding, row_raw_norm, col_raw_norm, row_labels, col_labels)
    return {
        "row_embedding": row_embedding,
        "col_embedding": col_embedding,
        "row_labels": row_labels,
        "col_labels": col_labels,
        "importance": importance,
        "singular_values": singular_values.astype(float).tolist(),
        "tau": tau,
        "row_column_adjusted_rand": float(adjusted_rand_score(row_labels, col_labels)),
        "embedding_asymmetry_rmse": float(np.sqrt(np.mean((row_embedding - col_embedding) ** 2))),
    }


def row_normalize(x: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(x, axis=1)
    norms[norms == 0] = 1.0
    return x / norms[:, None]


def compute_importance(
    gene_names: list[str],
    row_embedding: np.ndarray,
    col_embedding: np.ndarray,
    row_raw_norm: np.ndarray,
    col_raw_norm: np.ndarray,
    row_labels: np.ndarray,
    col_labels: np.ndarray,
) -> pd.DataFrame:
    asymmetry = np.linalg.norm(row_embedding - col_embedding, axis=1)
    leverage = np.maximum(row_raw_norm, col_raw_norm)
    table = pd.DataFrame(
        {
            "gene": gene_names,
            "disim_asymmetry": asymmetry,
            "row_raw_norm": row_raw_norm,
            "col_raw_norm": col_raw_norm,
            "spectral_leverage": leverage,
            "row_cluster": row_labels.astype(int),
            "col_cluster": col_labels.astype(int),
        }
    )
    table["asymmetry_percentile"] = table["disim_asymmetry"].rank(pct=True)
    table["leverage_percentile"] = table["spectral_leverage"].rank(pct=True)
    table["importance_score"] = 0.7 * table["asymmetry_percentile"] + 0.3 * table["leverage_percentile"]
    return table.sort_values("importance_score", ascending=False).reset_index(drop=True)


def compute_cophenetic(distance_matrix: np.ndarray, labels: list[str]) -> dict[str, Any]:
    row_dist = pdist(distance_matrix, metric="euclidean")
    col_dist = pdist(distance_matrix.T, metric="euclidean")
    row_linkage = linkage(row_dist, method="average")
    col_linkage = linkage(col_dist, method="average")
    row_corr, row_coph = cophenet(row_linkage, row_dist)
    col_corr, col_coph = cophenet(col_linkage, col_dist)
    row_df = normalize_square(squareform(row_coph), labels)
    col_df = normalize_square(squareform(col_coph), labels)
    return {
        "row_cophenetic": row_df,
        "col_cophenetic": col_df,
        "row_cophenetic_corr": float(row_corr),
        "col_cophenetic_corr": float(col_corr),
    }


def normalize_square(values: np.ndarray, labels: list[str]) -> pd.DataFrame:
    vmin = float(np.nanmin(values))
    vmax = float(np.nanmax(values))
    out = np.zeros_like(values, dtype=float) if vmax <= vmin else (values - vmin) / (vmax - vmin)
    np.fill_diagonal(out, 0.0)
    return pd.DataFrame(out, index=labels, columns=labels)


def cluster_coste_modules(selected: pd.DataFrame, row_coph: pd.DataFrame, col_coph: pd.DataFrame, n_modules: int) -> pd.DataFrame:
    labels = selected["gene"].astype(str).tolist()
    combined = (row_coph.loc[labels, labels].astype(float) + col_coph.loc[labels, labels].astype(float)) / 2.0
    values = combined.to_numpy(dtype=float)
    values = (values + values.T) / 2.0
    np.fill_diagonal(values, 0.0)
    z = linkage(squareform(values, checks=False), method="average")
    module_labels = fcluster(z, t=n_modules, criterion="maxclust")
    out = selected.copy()
    out["coste_module_id"] = module_labels.astype(int)
    out["coste_module_name"] = out["coste_module_id"].map(lambda x: f"CF_COSTE_M{x:02d}")
    return out.sort_values(["coste_module_id", "importance_score"], ascending=[True, False])


def run_transcript_histoseg(
    bin_gene: np.ndarray,
    bin_coords: np.ndarray,
    occupied_grid: np.ndarray,
    gene_names: list[str],
    modules: pd.DataFrame,
    *,
    n_domains: int,
    spatial_weight: float,
    random_state: int,
    output_dir: Path,
) -> dict[str, Any]:
    gene_to_col = {gene: idx for idx, gene in enumerate(gene_names)}
    module_ids = sorted(modules["coste_module_id"].unique().astype(int))
    module_scores = np.zeros((bin_gene.shape[0], len(module_ids)), dtype=np.float32)
    for module_idx, module_id in enumerate(module_ids):
        module_genes = [gene for gene in modules.loc[modules["coste_module_id"] == module_id, "gene"].astype(str) if gene in gene_to_col]
        cols = [gene_to_col[gene] for gene in module_genes]
        if cols:
            module_scores[:, module_idx] = bin_gene[:, cols].sum(axis=1)
    bin_totals = bin_gene.sum(axis=1)
    bin_totals[bin_totals <= 0] = np.nan
    normalized = np.log1p((module_scores / bin_totals[:, None]) * 10000.0)
    normalized = np.nan_to_num(normalized, nan=0.0, posinf=0.0, neginf=0.0)
    scaled = StandardScaler().fit_transform(normalized)
    smoothed = smooth_occupied_scores(scaled, occupied_grid)
    xy = StandardScaler().fit_transform(bin_coords)
    features = np.hstack([smoothed, xy * spatial_weight]).astype(np.float32)
    domain_labels = MiniBatchKMeans(
        n_clusters=n_domains,
        random_state=random_state,
        batch_size=4096,
        n_init=30,
        reassignment_ratio=0.01,
    ).fit_predict(features) + 1
    domain_grid = np.zeros(occupied_grid.shape, dtype=np.int16)
    domain_grid[occupied_grid] = domain_labels.astype(np.int16)
    rows = pd.DataFrame(
        {
            "bin_id": np.arange(len(domain_labels), dtype=int),
            "x": bin_coords[:, 0],
            "y": bin_coords[:, 1],
            "n_transcripts": bin_gene.sum(axis=1),
            "histoseg_structure_id": domain_labels.astype(int),
            "histoseg_structure_name": [f"CellFree-HistoSeg-{x:02d}" for x in domain_labels],
        }
    )
    for module_idx, module_id in enumerate(module_ids):
        rows[f"CF_COSTE_M{module_id:02d}_score"] = smoothed[:, module_idx]
    rows.to_parquet(output_dir / "cellfree_histoseg_domains.parquet", index=False)
    rows.to_csv(output_dir / "cellfree_histoseg_domains.csv", index=False)
    np.save(output_dir / "cellfree_histoseg_domain_grid.npy", domain_grid)
    summary = summarize_cellfree_domains(rows, modules)
    summary.to_csv(output_dir / "cellfree_histoseg_domain_summary.csv", index=False)
    plot_cellfree_domains(rows, output_dir)
    return {
        "domain_rows": rows,
        "domain_grid": domain_grid,
        "summary": summary,
    }


def smooth_occupied_scores(score_values: np.ndarray, occupied_grid: np.ndarray) -> np.ndarray:
    grid_shape = occupied_grid.shape
    out = np.zeros_like(score_values, dtype=np.float32)
    for module_idx in range(score_values.shape[1]):
        grid = np.zeros(grid_shape, dtype=np.float32)
        grid[occupied_grid] = score_values[:, module_idx]
        mask = occupied_grid.astype(np.float32)
        numerator = gaussian_filter(grid * mask, sigma=1.2)
        denominator = gaussian_filter(mask, sigma=1.2)
        denominator[denominator <= 1e-6] = 1.0
        smoothed = numerator / denominator
        out[:, module_idx] = smoothed[occupied_grid]
    return out


def summarize_cellfree_domains(domains: pd.DataFrame, modules: pd.DataFrame) -> pd.DataFrame:
    score_cols = [c for c in domains.columns if c.startswith("CF_COSTE_M") and c.endswith("_score")]
    top_gene_lookup = (
        modules.sort_values("importance_score", ascending=False)
        .groupby("coste_module_name")["gene"]
        .apply(lambda x: "/".join(x.head(5).astype(str)))
        .to_dict()
    )
    rows: list[dict[str, Any]] = []
    for domain_id, group in domains.groupby("histoseg_structure_id"):
        mean_scores = group[score_cols].mean().sort_values(ascending=False)
        top_score_col = str(mean_scores.index[0])
        top_module = top_score_col.replace("_score", "")
        rows.append(
            {
                "histoseg_structure_id": int(domain_id),
                "histoseg_structure_name": str(group["histoseg_structure_name"].iloc[0]),
                "n_bins": int(len(group)),
                "fraction_bins": float(len(group) / len(domains)),
                "n_transcripts": float(group["n_transcripts"].sum()),
                "top_coste_module": top_module,
                "top_coste_module_mean_score": float(mean_scores.iloc[0]),
                "top_coste_module_genes": top_gene_lookup.get(top_module, ""),
            }
        )
    return pd.DataFrame(rows).sort_values("histoseg_structure_id")


def plot_cellfree_domains(domains: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
    scatter = ax.scatter(
        domains["x"],
        domains["y"],
        c=domains["histoseg_structure_id"],
        s=9,
        cmap="tab20",
        linewidths=0,
    )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("Cell-free transcript HistoSeg domains")
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("histoseg domain")
    fig.savefig(output_dir / "cellfree_histoseg_domains.png", dpi=300)
    fig.savefig(output_dir / "cellfree_histoseg_domains.svg")
    plt.close(fig)


def compare_performance(output_dir: Path, timed_rows: list[dict[str, Any]], args: argparse.Namespace) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in timed_rows:
        rows.append({"workflow": "cell_free_transcript", **row})
    cell_based_summary = args.cell_based_disim_dir / "method_summary.csv"
    if cell_based_summary.exists():
        cell_based = pd.read_csv(cell_based_summary)
        for _, row in cell_based.iterrows():
            rows.append(
                {
                    "workflow": "cell_based",
                    "step": str(row["step"]),
                    "method": str(row["method"]),
                    "n_items": int(row["n_items"]),
                    "seconds": float(row["seconds"]),
                    "peak_rss_gb": float(row["peak_rss_gb"]),
                }
            )
    histoseg_manifest = args.cell_based_histoseg_dir / "histoseg_manifest.json"
    if histoseg_manifest.exists():
        manifest = json.loads(histoseg_manifest.read_text(encoding="utf-8"))
        rows.append(
            {
                "workflow": "cell_based",
                "step": "histoseg",
                "method": "cell_based_COSTE_HistoSeg",
                "n_items": int(manifest.get("n_cells", 0)),
                "seconds": float(manifest.get("seconds", np.nan)),
                "peak_rss_gb": np.nan,
            }
        )
    comparison = pd.DataFrame(rows)
    comparison.to_csv(output_dir / "performance_comparison.csv", index=False)
    return comparison


def run(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    timed_rows: list[dict[str, Any]] = []

    scan_result = timed(
        "transcript_metadata_scan",
        scan_transcript_metadata,
        args.transcripts_zarr,
        qv_min=args.qv_min,
        chunk_rows=args.chunk_rows,
        max_transcripts=args.max_transcripts,
    )
    metadata = scan_result.payload
    gene_names = list(metadata["gene_names"])
    identities = np.asarray(metadata["identities"], dtype=np.int64)
    transcript_counts = np.asarray(metadata["counts_by_identity"], dtype=np.int64)
    timed_rows.append(timed_row(scan_result, "cell_free_transcript", len(gene_names)))

    aggregate_result = timed(
        "transcript_grid_aggregation",
        aggregate_transcripts_to_grid,
        args.transcripts_zarr,
        identities=identities,
        bounds=metadata["bounds"],
        qv_min=args.qv_min,
        grid_bins=args.grid_bins,
        chunk_rows=args.chunk_rows,
        max_transcripts=args.max_transcripts,
    )
    bin_gene, bin_coords, occupied_grid = aggregate_result.payload
    keep = bin_gene.sum(axis=0) > 0
    identities = identities[keep]
    transcript_counts = transcript_counts[keep]
    gene_names = [gene for gene, ok in zip(gene_names, keep) if bool(ok)]
    bin_gene = bin_gene[:, keep]
    timed_rows.append(timed_row(aggregate_result, "cell_free_transcript", len(gene_names)))

    feature_summary = pd.DataFrame(
        {
            "gene_identity": identities.astype(int),
            "gene": gene_names,
            "transcript_count": transcript_counts.astype(int),
        }
    ).sort_values("transcript_count", ascending=False)
    feature_summary.to_csv(args.output_dir / "cellfree_transcript_gene_summary.csv", index=False)

    matrix_result = timed(
        "cellfree_directed_matrix",
        compute_directed_matrix,
        bin_gene,
        bin_coords,
        target_batch_size=args.target_batch_size,
        use_gpu=not args.cpu_only,
    )
    distance_matrix = matrix_result.payload
    timed_rows.append(timed_row(matrix_result, "cell_free_transcript", len(gene_names)))
    np.save(args.output_dir / "cellfree_directed_matrix.float32.npy", distance_matrix.astype(np.float32, copy=False))
    pd.Series(gene_names, name="gene").to_csv(args.output_dir / "cellfree_directed_matrix.genes.csv", index=False)

    disim_result = timed(
        "cellfree_full_disim",
        run_disim_array,
        distance_matrix,
        gene_names,
        n_components=args.n_components,
        n_clusters=args.n_clusters,
        random_state=args.random_state,
    )
    disim = disim_result.payload
    timed_rows.append(timed_row(disim_result, "cell_free_transcript", len(gene_names)))

    row_emb = pd.DataFrame(disim["row_embedding"], index=gene_names)
    col_emb = pd.DataFrame(disim["col_embedding"], index=gene_names)
    row_emb.to_csv(args.output_dir / "cellfree_full_disim_row_embedding.csv")
    col_emb.to_csv(args.output_dir / "cellfree_full_disim_col_embedding.csv")
    disim["importance"] = disim["importance"].merge(feature_summary, on="gene", how="left")
    disim["importance"].to_csv(args.output_dir / "cellfree_full_disim_gene_importance.csv", index=False)

    eligible = disim["importance"].loc[disim["importance"]["transcript_count"] >= args.min_selected_transcripts].copy()
    selected = eligible.head(min(args.selected_coste_genes, len(eligible))).copy()
    selected_genes = selected["gene"].astype(str).tolist()
    selected_cols = [gene_names.index(gene) for gene in selected_genes]
    selected_matrix = distance_matrix[np.ix_(selected_cols, selected_cols)]
    selected.to_csv(args.output_dir / f"cellfree_selected_top{len(selected)}_genes_for_coste.csv", index=False)
    np.save(args.output_dir / f"cellfree_selected_top{len(selected)}_directed_matrix.float32.npy", selected_matrix.astype(np.float32, copy=False))

    coste_result = timed("cellfree_selected_coste", compute_cophenetic, selected_matrix, selected_genes)
    coste = coste_result.payload
    timed_rows.append(timed_row(coste_result, "cell_free_transcript", len(selected_genes)))
    coste["row_cophenetic"].to_csv(args.output_dir / f"cellfree_selected_top{len(selected)}_coste_row_cophenetic.csv")
    coste["col_cophenetic"].to_csv(args.output_dir / f"cellfree_selected_top{len(selected)}_coste_col_cophenetic.csv")

    modules = cluster_coste_modules(selected, coste["row_cophenetic"], coste["col_cophenetic"], n_modules=args.n_gene_modules)
    modules.to_csv(args.output_dir / "cellfree_coste_gene_modules.csv", index=False)

    histoseg_result = timed(
        "cellfree_histoseg",
        run_transcript_histoseg,
        bin_gene,
        bin_coords,
        occupied_grid,
        gene_names,
        modules,
        n_domains=args.n_domains,
        spatial_weight=args.spatial_weight,
        random_state=args.random_state,
        output_dir=args.output_dir,
    )
    histoseg = histoseg_result.payload
    timed_rows.append(timed_row(histoseg_result, "cell_free_transcript", len(histoseg["domain_rows"])))

    pd.DataFrame(timed_rows).to_csv(args.output_dir / "method_summary.csv", index=False)
    comparison = compare_performance(args.output_dir, timed_rows, args)

    manifest = {
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "workflow": "cell_free_transcript_disim_coste_histoseg",
        "transcripts_zarr": str(args.transcripts_zarr),
        "qv_min": float(args.qv_min),
        "grid_bins": int(args.grid_bins),
        "max_transcripts": int(args.max_transcripts),
        "n_rows_scanned": int(metadata["n_rows_scanned"]),
        "n_transcripts_kept": int(metadata["n_transcripts_kept"]),
        "n_gene_expression_transcripts_kept": int(metadata["n_gene_expression_transcripts_kept"]),
        "n_non_gene_identities_seen": int(metadata["n_non_gene_identities_seen"]),
        "n_genes": int(len(gene_names)),
        "n_occupied_bins": int(bin_gene.shape[0]),
        "selected_coste_gene_count": int(len(selected_genes)),
        "selected_coste_min_transcripts": int(args.min_selected_transcripts),
        "selected_coste_eligible_gene_count": int(len(eligible)),
        "coste_selected": {
            "row_cophenetic_corr": coste["row_cophenetic_corr"],
            "col_cophenetic_corr": coste["col_cophenetic_corr"],
        },
        "disim_full": {
            "tau": float(disim["tau"]),
            "singular_values": disim["singular_values"],
            "row_column_adjusted_rand": float(disim["row_column_adjusted_rand"]),
            "embedding_asymmetry_rmse": float(disim["embedding_asymmetry_rmse"]),
        },
        "histoseg": {
            "n_domains": int(args.n_domains),
            "n_gene_modules": int(args.n_gene_modules),
            "n_bins": int(len(histoseg["domain_rows"])),
        },
        "outputs": {
            "method_summary": "method_summary.csv",
            "performance_comparison": "performance_comparison.csv",
            "domain_summary": "cellfree_histoseg_domain_summary.csv",
            "domain_grid": "cellfree_histoseg_domain_grid.npy",
            "domain_figure": "cellfree_histoseg_domains.png",
        },
    }
    _write_json(args.output_dir / "cellfree_manifest.json", manifest)
    print(json.dumps({"manifest": manifest, "performance_comparison": comparison.to_dict(orient="records")}, indent=2), flush=True)


def timed_row(result: TimedResult, method: str, n_items: int) -> dict[str, Any]:
    return {
        "step": result.step,
        "method": method,
        "n_items": int(n_items),
        "seconds": float(result.seconds),
        "peak_rss_gb": float(result.peak_rss_gb),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Cell-free transcript-only di-sim -> COSTE -> HistoSeg pipeline.")
    parser.add_argument("--transcripts-zarr", type=Path, default=DEFAULT_TRANSCRIPTS_ZARR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cell-based-disim-dir", type=Path, default=DEFAULT_CELL_BASED_DISIM_DIR)
    parser.add_argument("--cell-based-histoseg-dir", type=Path, default=DEFAULT_CELL_BASED_HISTOSEG_DIR)
    parser.add_argument("--qv-min", type=float, default=20.0)
    parser.add_argument("--grid-bins", type=int, default=96)
    parser.add_argument("--chunk-rows", type=int, default=2_000_000)
    parser.add_argument("--max-transcripts", type=int, default=0, help="0 means all transcripts; positive values are for smoke tests.")
    parser.add_argument("--target-batch-size", type=int, default=128)
    parser.add_argument("--cpu-only", action="store_true")
    parser.add_argument("--n-components", type=int, default=16)
    parser.add_argument("--n-clusters", type=int, default=16)
    parser.add_argument("--selected-coste-genes", type=int, default=1024)
    parser.add_argument("--min-selected-transcripts", type=int, default=5000)
    parser.add_argument("--n-gene-modules", type=int, default=16)
    parser.add_argument("--n-domains", type=int, default=12)
    parser.add_argument("--spatial-weight", type=float, default=0.35)
    parser.add_argument("--random-state", type=int, default=7)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
