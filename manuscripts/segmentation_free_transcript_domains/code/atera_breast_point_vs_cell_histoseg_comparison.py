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
from sklearn.neighbors import NearestNeighbors


DEFAULT_CELL_HISTOSEG_DIR = Path("/data/taobo.hu/atera_breast_coste1024_histoseg/coste1024_modules16_domains12")
DEFAULT_POINT_HISTOSEG_DIR = Path(
    "/data/taobo.hu/atera_breast_transcript_point_histoseg/coste1024_point_knn_train750k_assign_all"
)
DEFAULT_OUTPUT_DIR = Path("/data/taobo.hu/atera_breast_histoseg_comparison/point_vs_cell_coste1024")


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


class DistanceReservoir:
    def __init__(self, max_values: int, seed: int) -> None:
        self.max_values = int(max_values)
        self.rng = np.random.default_rng(seed)
        self.keys = np.empty(0, dtype=np.float32)
        self.values = np.empty(0, dtype=np.float32)

    def add(self, values: np.ndarray) -> None:
        if self.max_values <= 0 or len(values) == 0:
            return
        keys = self.rng.random(len(values), dtype=np.float32)
        self.keys = np.concatenate([self.keys, keys])
        self.values = np.concatenate([self.values, values.astype(np.float32, copy=False)])
        if len(self.values) > self.max_values:
            keep = np.argpartition(self.keys, self.max_values - 1)[: self.max_values]
            keep = keep[np.argsort(self.keys[keep])]
            self.keys = self.keys[keep]
            self.values = self.values[keep]


def timed(step: str, func, *args, **kwargs) -> TimedResult:
    t0 = time.perf_counter()
    with PeakMemoryTracker() as tracker:
        payload = func(*args, **kwargs)
    return TimedResult(step=step, seconds=time.perf_counter() - t0, peak_rss_gb=tracker.peak_rss / (1024**3), payload=payload)


def _log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def load_inputs(cell_dir: Path, point_dir: Path) -> dict[str, Any]:
    cell_domains = pd.read_parquet(
        cell_dir / "histoseg_domains.parquet",
        columns=["cell_id", "x", "y", "celltype", "histoseg_structure_id", "histoseg_structure_name"],
    )
    cell_summary = pd.read_csv(cell_dir / "histoseg_domain_summary.csv")
    point_summary = pd.read_csv(point_dir / "transcript_point_histoseg_domain_summary.csv")
    point_preview = pd.read_parquet(point_dir / "transcript_point_histoseg_preview_points.parquet")
    chunk_dir = point_dir / "transcript_point_domain_chunks"
    chunk_files = sorted(chunk_dir.glob("part_*.parquet"))
    if not chunk_files:
        raise FileNotFoundError(f"No point-domain parquet chunks found in {chunk_dir}")
    return {
        "cell_domains": cell_domains,
        "cell_summary": cell_summary,
        "point_summary": point_summary,
        "point_preview": point_preview,
        "chunk_files": chunk_files,
    }


def compare_domains(inputs: dict[str, Any], *, distance_sample_size: int, random_state: int) -> dict[str, Any]:
    cells = inputs["cell_domains"]
    point_summary = inputs["point_summary"]
    chunk_files = inputs["chunk_files"]

    point_domain_ids = sorted(int(x) for x in point_summary["histoseg_structure_id"].unique())
    cell_domain_ids = sorted(int(x) for x in cells["histoseg_structure_id"].unique())
    max_point_domain = max(point_domain_ids)
    max_cell_domain = max(cell_domain_ids)
    point_domain_count = max_point_domain + 1
    cell_domain_count = max_cell_domain + 1

    celltype_codes, celltype_names = pd.factorize(cells["celltype"].astype(str), sort=True)
    n_celltypes = len(celltype_names)

    cell_coords = cells[["x", "y"]].to_numpy(dtype=np.float32)
    cell_domain = cells["histoseg_structure_id"].to_numpy(dtype=np.int16)
    _log(f"building nearest-cell index for {len(cells):,} cell-domain centroids")
    nbrs = NearestNeighbors(n_neighbors=1, algorithm="kd_tree", leaf_size=40, n_jobs=-1)
    nbrs.fit(cell_coords)

    overlap = np.zeros((point_domain_count, cell_domain_count), dtype=np.int64)
    point_celltype = np.zeros((point_domain_count, n_celltypes), dtype=np.int64)
    distance_reservoir = DistanceReservoir(distance_sample_size, seed=random_state)
    distance_sum = 0.0
    distance_sq_sum = 0.0
    distance_max = 0.0
    n_points = 0

    for idx, path in enumerate(chunk_files, start=1):
        points = pd.read_parquet(path, columns=["x", "y", "histoseg_structure_id"])
        coords = points[["x", "y"]].to_numpy(dtype=np.float32)
        distances, nearest_idx = nbrs.kneighbors(coords, return_distance=True)
        nearest_idx = nearest_idx[:, 0]
        distances = distances[:, 0].astype(np.float32, copy=False)
        point_domain = points["histoseg_structure_id"].to_numpy(dtype=np.int16)
        nearest_cell_domain = cell_domain[nearest_idx]
        flat = point_domain.astype(np.int64) * cell_domain_count + nearest_cell_domain.astype(np.int64)
        overlap += np.bincount(flat, minlength=point_domain_count * cell_domain_count).reshape(
            point_domain_count, cell_domain_count
        )
        nearest_celltype = celltype_codes[nearest_idx]
        flat_celltype = point_domain.astype(np.int64) * n_celltypes + nearest_celltype.astype(np.int64)
        point_celltype += np.bincount(flat_celltype, minlength=point_domain_count * n_celltypes).reshape(
            point_domain_count, n_celltypes
        )
        distance_reservoir.add(distances)
        distance_sum += float(distances.sum(dtype=np.float64))
        distance_sq_sum += float(np.dot(distances.astype(np.float64), distances.astype(np.float64)))
        distance_max = max(distance_max, float(distances.max(initial=0.0)))
        n_points += int(len(points))
        if idx == 1 or idx % 25 == 0 or idx == len(chunk_files):
            _log(f"matched point chunks {idx:,}/{len(chunk_files):,}: {n_points:,} transcript points")

    q_values = [0.5, 0.75, 0.9, 0.95, 0.99]
    distance_quantiles = (
        {f"q{int(q * 100):02d}_um": float(np.quantile(distance_reservoir.values, q)) for q in q_values}
        if len(distance_reservoir.values)
        else {}
    )
    distance_stats = {
        "mean_um": float(distance_sum / max(n_points, 1)),
        "std_um": float(math.sqrt(max(distance_sq_sum / max(n_points, 1) - (distance_sum / max(n_points, 1)) ** 2, 0.0))),
        "max_um": float(distance_max),
        "quantile_sample_size": int(len(distance_reservoir.values)),
        **distance_quantiles,
    }
    return {
        "overlap": overlap,
        "point_celltype": point_celltype,
        "celltype_names": list(map(str, celltype_names)),
        "distance_stats": distance_stats,
        "n_points": int(n_points),
    }


def contingency_metrics(overlap: np.ndarray) -> dict[str, float]:
    table = overlap[1:, 1:].astype(np.float64)
    n = float(table.sum())
    if n <= 0:
        return {}
    row_sum = table.sum(axis=1)
    col_sum = table.sum(axis=0)
    comb = lambda x: x * (x - 1.0) / 2.0
    sum_comb = float(comb(table).sum())
    row_comb = float(comb(row_sum).sum())
    col_comb = float(comb(col_sum).sum())
    total_comb = float(comb(n))
    expected = row_comb * col_comb / total_comb if total_comb > 0 else 0.0
    max_index = 0.5 * (row_comb + col_comb)
    ari = (sum_comb - expected) / (max_index - expected) if max_index > expected else 0.0

    p = table / n
    p_row = row_sum / n
    p_col = col_sum / n
    nz = p > 0
    mi = float((p[nz] * np.log(p[nz] / (p_row[:, None] * p_col[None, :])[nz])).sum())
    h_row = float(-(p_row[p_row > 0] * np.log(p_row[p_row > 0])).sum())
    h_col = float(-(p_col[p_col > 0] * np.log(p_col[p_col > 0])).sum())
    nmi = 2.0 * mi / (h_row + h_col) if h_row + h_col > 0 else 0.0
    homogeneity = mi / h_col if h_col > 0 else 0.0
    completeness = mi / h_row if h_row > 0 else 0.0
    point_to_cell_purity = float(table.max(axis=1).sum() / n)
    cell_to_point_purity = float(table.max(axis=0).sum() / n)
    return {
        "n_transcript_points_compared": float(n),
        "adjusted_rand_index": float(ari),
        "normalized_mutual_info": float(nmi),
        "homogeneity_point_given_cell": float(homogeneity),
        "completeness_cell_given_point": float(completeness),
        "point_to_cell_purity": point_to_cell_purity,
        "cell_to_point_purity": cell_to_point_purity,
    }


def parse_gene_set(value: Any) -> set[str]:
    if pd.isna(value):
        return set()
    return {x.strip() for x in str(value).replace(",", "/").split("/") if x.strip()}


def classify_pathology(celltype: str = "", genes: str = "") -> str:
    text = f"{celltype} {genes}".lower()
    immune = any(
        term in text for term in ["lymph", "myeloid", "immune", "b cell", "plasma", "ccl", "cxcl", "iglc", "ms4a1", "ltb", "mmp9"]
    )
    stroma = any(term in text for term in ["caf", "fibro", "strom", "endothelial", "vascular", "fn1", "mmp11", "dpt", "col", "igf2"])
    tumor = any(term in text for term in ["invasive", "tumor", "11q13", "serpina", "msmb", "kcnq3", "ccnd1"])
    luminal = any(term in text for term in ["dcis", "luminal", "apocrine", "pip", "tff", "hspb8", "clic6", "niban", "tat"])
    blood = any(term in text for term in ["hbb", "hemoglobin"])
    if immune and stroma:
        return "immune-stromal interface"
    if tumor and immune:
        return "tumor/immune interface"
    if tumor and stroma:
        return "tumor/stroma interface"
    if immune:
        return "immune-rich"
    if stroma:
        return "stroma/CAF"
    if tumor:
        return "invasive/tumor-rich"
    if luminal:
        return "luminal/DCIS/apocrine"
    if blood:
        return "blood/RBC-like"
    return "mixed/other"


def build_overlap_outputs(
    overlap: np.ndarray,
    point_celltype: np.ndarray,
    celltype_names: list[str],
    inputs: dict[str, Any],
    output_dir: Path,
) -> dict[str, Any]:
    point_summary = inputs["point_summary"].copy()
    cell_summary = inputs["cell_summary"].copy()

    point_summary["point_pathology_label"] = [
        classify_pathology(genes=row.get("top_coste_module_genes", "")) for _, row in point_summary.iterrows()
    ]
    cell_summary["cell_pathology_label"] = [
        classify_pathology(row.get("dominant_celltype", ""), row.get("top_coste_module_genes", ""))
        for _, row in cell_summary.iterrows()
    ]

    point_names = {
        int(row.histoseg_structure_id): str(row.histoseg_structure_name) for _, row in point_summary.iterrows()
    }
    cell_names = {
        int(row.histoseg_structure_id): str(row.histoseg_structure_name) for _, row in cell_summary.iterrows()
    }
    point_ids = sorted(point_names)
    cell_ids = sorted(cell_names)
    table = overlap[np.ix_(point_ids, cell_ids)]
    row_sum = table.sum(axis=1)
    col_sum = table.sum(axis=0)
    total = int(table.sum())
    counts_df = pd.DataFrame(
        table,
        index=[point_names[i] for i in point_ids],
        columns=[cell_names[j] for j in cell_ids],
    )
    row_fraction = counts_df.div(np.maximum(row_sum, 1), axis=0)
    col_fraction = counts_df.div(np.maximum(col_sum, 1), axis=1)
    counts_df.to_csv(output_dir / "point_vs_cell_overlap_counts.csv")
    row_fraction.to_csv(output_dir / "point_vs_cell_overlap_point_fraction.csv")
    col_fraction.to_csv(output_dir / "point_vs_cell_overlap_cell_fraction.csv")

    point_by_id = point_summary.set_index("histoseg_structure_id")
    cell_by_id = cell_summary.set_index("histoseg_structure_id")
    best_rows: list[dict[str, Any]] = []
    for row_idx, point_id in enumerate(point_ids):
        counts = table[row_idx]
        best_col_idx = int(np.argmax(counts))
        cell_id = cell_ids[best_col_idx]
        count = int(counts[best_col_idx])
        point_total = int(row_sum[row_idx])
        cell_total = int(col_sum[best_col_idx])
        point_row = point_by_id.loc[point_id]
        cell_row = cell_by_id.loc[cell_id]
        top_celltype_code = int(np.argmax(point_celltype[point_id])) if point_celltype[point_id].sum() else -1
        top_celltype = celltype_names[top_celltype_code] if top_celltype_code >= 0 else ""
        top_celltype_fraction = (
            float(point_celltype[point_id, top_celltype_code] / max(point_celltype[point_id].sum(), 1))
            if top_celltype_code >= 0
            else 0.0
        )
        point_genes = str(point_row.get("top_coste_module_genes", ""))
        cell_genes = str(cell_row.get("top_coste_module_genes", ""))
        gene_union = parse_gene_set(point_genes) | parse_gene_set(cell_genes)
        gene_jaccard = len(parse_gene_set(point_genes) & parse_gene_set(cell_genes)) / max(len(gene_union), 1)
        point_label = str(point_row.get("point_pathology_label", "mixed/other"))
        cell_label = str(cell_row.get("cell_pathology_label", "mixed/other"))
        label_relation = "same" if point_label == cell_label else ("related" if "tumor" in point_label + cell_label else "different")
        best_rows.append(
            {
                "point_domain_id": int(point_id),
                "point_domain": point_names[point_id],
                "point_transcripts": point_total,
                "point_fraction_total": float(point_total / max(total, 1)),
                "point_top_module": str(point_row.get("top_coste_module", "")),
                "point_top_genes": point_genes,
                "point_pathology_label": point_label,
                "best_cell_domain_id": int(cell_id),
                "best_cell_domain": cell_names[cell_id],
                "overlap_transcripts": count,
                "point_to_cell_fraction": float(count / max(point_total, 1)),
                "cell_domain_coverage": float(count / max(cell_total, 1)),
                "jaccard_transcript_weighted": float(count / max(point_total + cell_total - count, 1)),
                "matched_cell_dominant_celltype": str(cell_row.get("dominant_celltype", "")),
                "matched_cell_dominant_celltype_fraction": float(cell_row.get("dominant_celltype_fraction", np.nan)),
                "cell_pathology_label": cell_label,
                "nearest_celltype_within_point_domain": top_celltype,
                "nearest_celltype_fraction_within_point_domain": top_celltype_fraction,
                "matched_cell_top_module": str(cell_row.get("top_coste_module", "")),
                "matched_cell_top_genes": cell_genes,
                "top_gene_jaccard": float(gene_jaccard),
                "pathology_label_relation": label_relation,
            }
        )
    best = pd.DataFrame(best_rows).sort_values(["point_to_cell_fraction", "jaccard_transcript_weighted"], ascending=False)
    best.to_csv(output_dir / "point_vs_cell_best_domain_matches.csv", index=False)

    metrics = contingency_metrics(overlap)
    metrics["weighted_mean_best_jaccard"] = float(
        np.average(best["jaccard_transcript_weighted"], weights=best["point_transcripts"])
    )
    metrics["weighted_mean_point_to_cell_fraction"] = float(
        np.average(best["point_to_cell_fraction"], weights=best["point_transcripts"])
    )
    pd.DataFrame([metrics]).to_csv(output_dir / "point_vs_cell_global_metrics.csv", index=False)
    _write_json(output_dir / "point_vs_cell_global_metrics.json", metrics)
    return {
        "counts": counts_df,
        "row_fraction": row_fraction,
        "col_fraction": col_fraction,
        "best_matches": best,
        "metrics": metrics,
        "point_summary": point_summary,
        "cell_summary": cell_summary,
    }


def match_preview_points(inputs: dict[str, Any], output_dir: Path) -> pd.DataFrame:
    cells = inputs["cell_domains"]
    preview = inputs["point_preview"].copy()
    nbrs = NearestNeighbors(n_neighbors=1, algorithm="kd_tree", leaf_size=40, n_jobs=-1)
    nbrs.fit(cells[["x", "y"]].to_numpy(dtype=np.float32))
    distances, nearest_idx = nbrs.kneighbors(preview[["x", "y"]].to_numpy(dtype=np.float32), return_distance=True)
    nearest = cells.iloc[nearest_idx[:, 0]].reset_index(drop=True)
    preview["nearest_cell_domain_id"] = nearest["histoseg_structure_id"].to_numpy(dtype=int)
    preview["nearest_cell_domain_name"] = nearest["histoseg_structure_name"].astype(str).to_numpy()
    preview["nearest_celltype"] = nearest["celltype"].astype(str).to_numpy()
    preview["nearest_cell_distance_um"] = distances[:, 0].astype(np.float32)
    preview.to_parquet(output_dir / "point_preview_nearest_cell_domains.parquet", index=False)
    preview.to_csv(output_dir / "point_preview_nearest_cell_domains.csv", index=False)
    return preview


def plot_heatmap(row_fraction: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    image = ax.imshow(row_fraction.to_numpy(dtype=float), vmin=0, vmax=max(0.5, float(row_fraction.max().max())), cmap="viridis")
    ax.set_xticks(np.arange(row_fraction.shape[1]))
    ax.set_xticklabels([c.replace("COSTE-HistoSeg-", "C") for c in row_fraction.columns], rotation=45, ha="right")
    ax.set_yticks(np.arange(row_fraction.shape[0]))
    ax.set_yticklabels([i.replace("TranscriptPoint-HistoSeg-", "P") for i in row_fraction.index])
    ax.set_xlabel("cell-based COSTE-HistoSeg domain")
    ax.set_ylabel("transcript-point HistoSeg domain")
    ax.set_title("Transcript-point to cell-domain overlap (row fraction)")
    for i in range(row_fraction.shape[0]):
        for j in range(row_fraction.shape[1]):
            value = float(row_fraction.iloc[i, j])
            if value >= 0.05:
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="white" if value > 0.25 else "black", fontsize=9)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("fraction of point-domain transcripts")
    fig.savefig(output_dir / "point_vs_cell_overlap_heatmap.png", dpi=300)
    fig.savefig(output_dir / "point_vs_cell_overlap_heatmap.svg")
    plt.close(fig)


def plot_spatial_comparison(inputs: dict[str, Any], preview: pd.DataFrame, output_dir: Path) -> None:
    cells = inputs["cell_domains"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5), constrained_layout=True)
    panels = [
        (
            axes[0],
            cells["x"],
            cells["y"],
            cells["histoseg_structure_id"],
            "Cell-based COSTE-HistoSeg\n170,057 cells",
            1.5,
            0.8,
        ),
        (
            axes[1],
            preview["x"],
            preview["y"],
            preview["histoseg_structure_id"],
            "Transcript-point HistoSeg\n500,000 point preview",
            0.35,
            0.75,
        ),
        (
            axes[2],
            preview["x"],
            preview["y"],
            preview["nearest_cell_domain_id"],
            "Point preview colored by\nnearest cell-based domain",
            0.35,
            0.75,
        ),
    ]
    # All three panels share an explicit [1, 12] normalization so the single
    # "domain id" colour key is consistent across them. Without a fixed
    # vmin/vmax each scatter would auto-normalize to its own data range, and a
    # panel that happened not to span the full 1..12 (e.g. a nearest-cell-domain
    # that never reaches a given id) would map the same colour to a different
    # id, silently breaking the shared legend. NOTE: panel 2 is numbered by
    # transcript-POINT domain id while panels 1 and 3 are numbered by CELL
    # domain id (two independent clusterings); the shared colour scale is
    # intentional for spatial comparison, and the panel titles disambiguate.
    for ax, x, y, color, title, size, alpha in panels:
        scatter = ax.scatter(x, y, c=color, s=size, cmap="tab20", vmin=1, vmax=12, linewidths=0, alpha=alpha, rasterized=True)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
        ax.set_title(title)
    cbar = fig.colorbar(scatter, ax=axes, fraction=0.02, pad=0.01)
    cbar.set_label("domain id")
    fig.savefig(output_dir / "point_vs_cell_spatial_comparison.png", dpi=300)
    fig.savefig(output_dir / "point_vs_cell_spatial_comparison.svg")
    plt.close(fig)


def write_report(outputs: dict[str, Any], distance_stats: dict[str, Any], output_dir: Path) -> None:
    metrics = outputs["metrics"]
    best = outputs["best_matches"].copy()
    top_matches = best.sort_values("point_fraction_total", ascending=False).head(8)
    lines = [
        "# Transcript-Point vs Cell-Based COSTE-HistoSeg Comparison",
        "",
        "Comparison unit: every selected-gene transcript point was matched to the nearest cell-based COSTE-HistoSeg cell centroid. The overlap is transcript-weighted and does not use a spatial grid.",
        "",
        "## Global overlap",
        "",
        f"- Transcript points compared: {int(metrics['n_transcript_points_compared']):,}",
        f"- Adjusted Rand index: {metrics['adjusted_rand_index']:.4f}",
        f"- Normalized mutual information: {metrics['normalized_mutual_info']:.4f}",
        f"- Point-to-cell purity: {metrics['point_to_cell_purity']:.4f}",
        f"- Cell-to-point purity: {metrics['cell_to_point_purity']:.4f}",
        f"- Weighted mean best Jaccard: {metrics['weighted_mean_best_jaccard']:.4f}",
        "",
        "Nearest-cell distance for transcript points:",
        "",
        f"- Mean: {distance_stats['mean_um']:.2f} um",
        f"- Median sample: {distance_stats.get('q50_um', float('nan')):.2f} um",
        f"- 95th percentile sample: {distance_stats.get('q95_um', float('nan')):.2f} um",
        "",
        "## Largest point-domain matches",
        "",
        "| point domain | point fraction | best cell domain | point-to-cell fraction | Jaccard | nearest dominant cell type | point genes | cell genes | interpretation |",
        "| --- | ---: | --- | ---: | ---: | --- | --- | --- | --- |",
    ]
    for _, row in top_matches.iterrows():
        lines.append(
            "| {point_domain} | {point_fraction_total:.3f} | {best_cell_domain} | {point_to_cell_fraction:.3f} | "
            "{jaccard_transcript_weighted:.3f} | {nearest_celltype_within_point_domain} | {point_top_genes} | "
            "{matched_cell_top_genes} | {point_pathology_label} -> {cell_pathology_label} ({pathology_label_relation}) |".format(
                **row.to_dict()
            )
        )
    lines.extend(
        [
            "",
            "## Pathology interpretation",
            "",
            "The point-level map is cell-free in computation, but its domains can be interpreted by projecting them onto the nearest cell-based molecular HistoSeg domains. High-purity matches mean a transcript-point domain mostly falls inside one cell-level molecular region. Low ARI/NMI with moderate purity means the two maps agree on broad regions but split and merge them differently. The pathology labels below are heuristic molecular labels based on dominant nearest cell type plus top COSTE genes; they are not a formal H&E diagnosis.",
            "",
            "The largest point domains are dominated by COSTE signatures such as SERPINA6/MSMB/KCNQ3/SERPINA1 and RERGL/BMPER/IGF2/DPT/CCL22. These signatures map mainly to tumor-rich, DCIS/luminal, stromal/CAF, immune-rich, and immune-stromal interface cell-level regions depending on their nearest cell-domain context. This supports the view that the point-level HistoSeg is detecting molecular micro-regions rather than simply reproducing the cell-level segmentation.",
            "",
        ]
    )
    (output_dir / "point_vs_cell_pathology_interpretation.md").write_text("\n".join(lines), encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    timed_rows: list[dict[str, Any]] = []
    load_result = timed("load_inputs", load_inputs, args.cell_histoseg_dir, args.point_histoseg_dir)
    inputs = load_result.payload
    timed_rows.append(timed_row(load_result, 0))

    compare_result = timed(
        "match_all_point_domains_to_nearest_cell_domain",
        compare_domains,
        inputs,
        distance_sample_size=args.distance_sample_size,
        random_state=args.random_state,
    )
    comparison = compare_result.payload
    timed_rows.append(timed_row(compare_result, comparison["n_points"]))

    output_result = timed(
        "write_overlap_outputs",
        build_overlap_outputs,
        comparison["overlap"],
        comparison["point_celltype"],
        comparison["celltype_names"],
        inputs,
        args.output_dir,
    )
    outputs = output_result.payload
    timed_rows.append(timed_row(output_result, int(comparison["n_points"])))

    preview_result = timed("match_preview_points", match_preview_points, inputs, args.output_dir)
    preview = preview_result.payload
    timed_rows.append(timed_row(preview_result, len(preview)))

    plot_heatmap(outputs["row_fraction"], args.output_dir)
    plot_spatial_comparison(inputs, preview, args.output_dir)
    write_report(outputs, comparison["distance_stats"], args.output_dir)

    method_summary = pd.DataFrame(timed_rows)
    method_summary.to_csv(args.output_dir / "method_summary.csv", index=False)
    manifest = {
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "workflow": "transcript_point_vs_cell_based_coste_histoseg_comparison",
        "cell_histoseg_dir": str(args.cell_histoseg_dir),
        "point_histoseg_dir": str(args.point_histoseg_dir),
        "comparison_unit": "selected_gene_transcript_point matched to nearest cell-based HistoSeg cell",
        "uses_spatial_grid": False,
        "n_cell_domains": int(inputs["cell_domains"]["histoseg_structure_id"].nunique()),
        "n_point_domains": int(inputs["point_summary"]["histoseg_structure_id"].nunique()),
        "n_cells": int(len(inputs["cell_domains"])),
        "n_point_chunks": int(len(inputs["chunk_files"])),
        "n_transcript_points_compared": int(comparison["n_points"]),
        "distance_stats": comparison["distance_stats"],
        "global_metrics": outputs["metrics"],
        "outputs": {
            "overlap_counts": "point_vs_cell_overlap_counts.csv",
            "overlap_point_fraction": "point_vs_cell_overlap_point_fraction.csv",
            "overlap_cell_fraction": "point_vs_cell_overlap_cell_fraction.csv",
            "best_matches": "point_vs_cell_best_domain_matches.csv",
            "global_metrics": "point_vs_cell_global_metrics.json",
            "pathology_report": "point_vs_cell_pathology_interpretation.md",
            "overlap_heatmap": "point_vs_cell_overlap_heatmap.png",
            "spatial_comparison": "point_vs_cell_spatial_comparison.png",
        },
    }
    _write_json(args.output_dir / "point_vs_cell_comparison_manifest.json", manifest)
    print(json.dumps({"manifest": manifest, "method_summary": method_summary.to_dict(orient="records")}, indent=2), flush=True)


def timed_row(result: TimedResult, n_items: int) -> dict[str, Any]:
    return {
        "step": result.step,
        "n_items": int(n_items),
        "seconds": float(result.seconds),
        "peak_rss_gb": float(result.peak_rss_gb),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare transcript-point HistoSeg domains with cell-based COSTE-HistoSeg.")
    parser.add_argument("--cell-histoseg-dir", type=Path, default=DEFAULT_CELL_HISTOSEG_DIR)
    parser.add_argument("--point-histoseg-dir", type=Path, default=DEFAULT_POINT_HISTOSEG_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--distance-sample-size", type=int, default=1_000_000)
    parser.add_argument("--random-state", type=int, default=7)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
