#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.ndimage import median_filter
from scipy.sparse import csc_matrix, csr_matrix
from scipy.spatial.distance import squareform
from sklearn.cluster import MiniBatchKMeans
from sklearn.neighbors import NearestNeighbors, kneighbors_graph
from sklearn.preprocessing import StandardScaler


DEFAULT_INPUT_DIR = Path("/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast")
DEFAULT_COSTE_DIR = Path(
    "/data/taobo.hu/atera_breast_disim_coste_benchmark/"
    "full_disim_18028_select1024_min5000_grid96_k16"
)
DEFAULT_OUTPUT_DIR = Path("/data/taobo.hu/atera_breast_coste1024_histoseg")


def _log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def _decode(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def load_cells(input_dir: Path) -> pd.DataFrame:
    cells = pd.read_csv(
        input_dir / "cells.csv.gz",
        usecols=["cell_id", "x_centroid", "y_centroid", "transcript_counts", "total_counts", "cell_area"],
    )
    groups = pd.read_csv(input_dir / "WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv")
    groups = groups.rename(columns={"group": "celltype", "Clusters": "celltype", "Barcode": "cell_id"})
    cells = cells.merge(groups[["cell_id", "celltype", "color"]], on="cell_id", how="left")
    cells = cells.rename(columns={"x_centroid": "x", "y_centroid": "y"})
    cells["celltype"] = cells["celltype"].fillna("Unassigned").astype(str)
    return cells


def load_selected_coste(coste_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = pd.read_csv(coste_dir / "selected_top1024_genes_for_coste.csv")
    row_coph = pd.read_csv(coste_dir / "selected_top1024_coste_row_cophenetic.csv", index_col=0)
    col_path = coste_dir / "selected_top1024_coste_col_cophenetic.csv"
    if col_path.exists():
        col_coph = pd.read_csv(col_path, index_col=0)
        common = [g for g in selected["gene"].astype(str) if g in row_coph.index and g in col_coph.index]
        combined = (row_coph.loc[common, common].astype(float) + col_coph.loc[common, common].astype(float)) / 2.0
    else:
        common = [g for g in selected["gene"].astype(str) if g in row_coph.index]
        combined = row_coph.loc[common, common].astype(float)
    selected = selected.loc[selected["gene"].astype(str).isin(common)].copy()
    selected["gene"] = selected["gene"].astype(str)
    selected = selected.set_index("gene").loc[common].reset_index()
    return selected, combined


def cluster_coste_genes(selected: pd.DataFrame, coste_distance: pd.DataFrame, n_modules: int) -> pd.DataFrame:
    values = coste_distance.to_numpy(dtype=float)
    values = np.nan_to_num(values, nan=np.nanmedian(values), posinf=np.nanmax(values), neginf=np.nanmin(values))
    values = (values + values.T) / 2.0
    np.fill_diagonal(values, 0.0)
    condensed = squareform(values, checks=False)
    z = linkage(condensed, method="average")
    labels = fcluster(z, t=n_modules, criterion="maxclust")
    modules = selected.copy()
    modules["coste_module_id"] = labels.astype(int)
    modules["coste_module_name"] = modules["coste_module_id"].map(lambda x: f"COSTE_M{x:02d}")
    return modules.sort_values(["coste_module_id", "importance_score"], ascending=[True, False])


def load_selected_expression(input_dir: Path, genes: list[str]) -> tuple[csr_matrix, list[str], list[str]]:
    matrix_h5 = input_dir / "cell_feature_matrix.h5"
    with h5py.File(matrix_h5, "r") as handle:
        shape = tuple(int(x) for x in handle["matrix/shape"][:])
        feature_names = _decode(handle["matrix/features/name"][:])
        feature_type = _decode(handle["matrix/features/feature_type"][:])
        data = np.asarray(handle["matrix/data"][:], dtype=np.float32)
        indices = np.asarray(handle["matrix/indices"][:], dtype=np.int64)
        indptr = np.asarray(handle["matrix/indptr"][:], dtype=np.int64)
        barcodes = _decode(handle["matrix/barcodes"][:])

    gene_to_index: dict[str, int] = {}
    for i, (name, kind) in enumerate(zip(feature_names, feature_type)):
        if kind == "Gene Expression" and name not in gene_to_index:
            gene_to_index[name] = i
    missing = [gene for gene in genes if gene not in gene_to_index]
    if missing:
        raise ValueError(f"Missing selected genes from H5 matrix: {missing[:10]} ({len(missing)} total)")
    feature_idx = np.array([gene_to_index[gene] for gene in genes], dtype=np.int64)
    feature_by_cell = csc_matrix((data, indices, indptr), shape=shape)
    cell_by_gene = feature_by_cell[feature_idx, :].T.tocsr()
    return cell_by_gene, genes, barcodes


def compute_module_scores(
    cell_by_gene: csr_matrix,
    genes: list[str],
    modules: pd.DataFrame,
    cells: pd.DataFrame,
) -> tuple[pd.DataFrame, np.ndarray]:
    gene_to_col = {gene: idx for idx, gene in enumerate(genes)}
    module_ids = sorted(modules["coste_module_id"].unique().astype(int))
    module_to_col = {module_id: idx for idx, module_id in enumerate(module_ids)}
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    for _, row in modules.iterrows():
        gene = str(row["gene"])
        if gene not in gene_to_col:
            continue
        rows.append(gene_to_col[gene])
        cols.append(module_to_col[int(row["coste_module_id"])])
        vals.append(1.0)
    indicator = csr_matrix((vals, (rows, cols)), shape=(len(genes), len(module_ids)), dtype=np.float32)
    module_counts = (cell_by_gene @ indicator).toarray().astype(np.float32)
    totals = pd.to_numeric(cells["total_counts"], errors="coerce").fillna(0).to_numpy(dtype=np.float32)
    totals[totals <= 0] = np.nan
    normalized = np.log1p((module_counts / totals[:, None]) * 10000.0)
    normalized = np.nan_to_num(normalized, nan=0.0, posinf=0.0, neginf=0.0)
    score_values = StandardScaler().fit_transform(normalized)
    score_cols = [f"COSTE_M{module_id:02d}_score" for module_id in module_ids]
    scores = pd.DataFrame(score_values, columns=score_cols)
    scores.insert(0, "cell_id", cells["cell_id"].to_numpy())
    return scores, module_counts


def smooth_scores(coords: np.ndarray, score_values: np.ndarray, n_neighbors: int) -> np.ndarray:
    graph = kneighbors_graph(coords, n_neighbors=n_neighbors, mode="connectivity", include_self=True, n_jobs=-1)
    graph = graph.tocsr().astype(np.float32)
    row_sum = np.asarray(graph.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    graph = graph.multiply(1.0 / row_sum[:, None])
    return graph @ score_values


def cluster_domains(
    cells: pd.DataFrame,
    module_scores: pd.DataFrame,
    *,
    n_domains: int,
    n_neighbors: int,
    spatial_weight: float,
    random_state: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    coords = cells[["x", "y"]].to_numpy(dtype=np.float32)
    score_values = module_scores.drop(columns=["cell_id"]).to_numpy(dtype=np.float32)
    smoothed = smooth_scores(coords, score_values, n_neighbors=n_neighbors)
    xy = StandardScaler().fit_transform(coords)
    features = np.hstack([smoothed, xy * spatial_weight]).astype(np.float32)
    labels = MiniBatchKMeans(
        n_clusters=n_domains,
        random_state=random_state,
        batch_size=8192,
        n_init=20,
        reassignment_ratio=0.01,
    ).fit_predict(features)
    domains = cells.copy()
    domains["histoseg_structure_id"] = labels.astype(int) + 1
    domains["histoseg_structure_name"] = domains["histoseg_structure_id"].map(lambda x: f"COSTE-HistoSeg-{x:02d}")
    domains["histoseg_signed_distance_um"] = approximate_boundary_distance(coords, labels)
    for idx, col in enumerate(module_scores.columns[1:]):
        domains[col] = smoothed[:, idx]
    return domains, smoothed


def approximate_boundary_distance(coords: np.ndarray, labels: np.ndarray) -> np.ndarray:
    distances = np.zeros(len(labels), dtype=np.float32)
    for label in np.unique(labels):
        in_domain = labels == label
        outside = ~in_domain
        if outside.sum() == 0:
            distances[in_domain] = np.nan
            continue
        nn = NearestNeighbors(n_neighbors=1, algorithm="auto").fit(coords[outside])
        dist, _ = nn.kneighbors(coords[in_domain])
        distances[in_domain] = dist[:, 0].astype(np.float32)
    return distances


def summarize_domains(domains: pd.DataFrame, modules: pd.DataFrame) -> pd.DataFrame:
    score_cols = [c for c in domains.columns if c.startswith("COSTE_M") and c.endswith("_score")]
    rows: list[dict[str, Any]] = []
    top_gene_lookup = (
        modules.sort_values("importance_score", ascending=False)
        .groupby("coste_module_name")["gene"]
        .apply(lambda x: "/".join(x.head(5).astype(str)))
        .to_dict()
    )
    for domain_id, group in domains.groupby("histoseg_structure_id"):
        mean_scores = group[score_cols].mean().sort_values(ascending=False)
        top_module_score_col = str(mean_scores.index[0])
        top_module = top_module_score_col.replace("_score", "")
        celltype_counts = group["celltype"].value_counts()
        rows.append(
            {
                "histoseg_structure_id": int(domain_id),
                "histoseg_structure_name": str(group["histoseg_structure_name"].iloc[0]),
                "n_cells": int(len(group)),
                "fraction_cells": float(len(group) / len(domains)),
                "median_boundary_distance_um": float(group["histoseg_signed_distance_um"].median()),
                "top_coste_module": top_module,
                "top_coste_module_mean_score": float(mean_scores.iloc[0]),
                "top_coste_module_genes": top_gene_lookup.get(top_module, ""),
                "dominant_celltype": str(celltype_counts.index[0]) if len(celltype_counts) else "",
                "dominant_celltype_fraction": float(celltype_counts.iloc[0] / len(group)) if len(celltype_counts) else 0.0,
            }
        )
    return pd.DataFrame(rows).sort_values("histoseg_structure_id")


def plot_domains(domains: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
    scatter = ax.scatter(
        domains["x"],
        domains["y"],
        c=domains["histoseg_structure_id"],
        s=1,
        cmap="tab20",
        linewidths=0,
        rasterized=True,
    )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("COSTE-selected HistoSeg domains")
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("histoseg domain")
    fig.savefig(output_dir / "coste_histoseg_domains.png", dpi=300)
    fig.savefig(output_dir / "coste_histoseg_domains.svg")
    plt.close(fig)


def write_grid_preview(domains: pd.DataFrame, output_dir: Path, bins: int = 256) -> None:
    x = domains["x"].to_numpy(dtype=float)
    y = domains["y"].to_numpy(dtype=float)
    labels = domains["histoseg_structure_id"].to_numpy(dtype=int)
    x_edges = np.linspace(x.min(), x.max(), bins + 1)
    y_edges = np.linspace(y.min(), y.max(), bins + 1)
    xi = np.clip(np.searchsorted(x_edges, x, side="right") - 1, 0, bins - 1)
    yi = np.clip(np.searchsorted(y_edges, y, side="right") - 1, 0, bins - 1)
    grid = np.zeros((bins, bins), dtype=np.int16)
    counts = np.zeros((bins, bins), dtype=np.int16)
    for xx, yy, label in zip(xi, yi, labels):
        if counts[yy, xx] == 0:
            grid[yy, xx] = label
        else:
            grid[yy, xx] = label if label == grid[yy, xx] else grid[yy, xx]
        counts[yy, xx] += 1
    grid = median_filter(grid, size=3)
    np.save(output_dir / "coste_histoseg_domain_grid.npy", grid)


def run(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    _log("loading COSTE-selected genes and cophenetic matrix")
    selected, coste_distance = load_selected_coste(args.coste_dir)
    modules = cluster_coste_genes(selected, coste_distance, n_modules=args.n_gene_modules)
    modules.to_csv(args.output_dir / "coste_gene_modules.csv", index=False)

    _log("loading Xenium cells and selected gene expression")
    cells = load_cells(args.input_dir)
    cell_by_gene, genes, barcodes = load_selected_expression(args.input_dir, modules["gene"].astype(str).tolist())
    cells = cells.set_index("cell_id").reindex(barcodes).reset_index()
    if cells[["x", "y"]].isna().any().any():
        raise ValueError("H5 barcodes and cells.csv.gz do not align")

    _log("computing COSTE module scores per cell")
    module_scores, _module_counts = compute_module_scores(cell_by_gene, genes, modules, cells)
    module_scores.to_parquet(args.output_dir / "coste_module_scores.parquet", index=False)

    _log("clustering spatial domains")
    domains, _smoothed = cluster_domains(
        cells,
        module_scores,
        n_domains=args.n_domains,
        n_neighbors=args.n_neighbors,
        spatial_weight=args.spatial_weight,
        random_state=args.random_state,
    )
    domains.to_parquet(args.output_dir / "histoseg_domains.parquet", index=False)
    domains[
        [
            "cell_id",
            "x",
            "y",
            "celltype",
            "histoseg_structure_id",
            "histoseg_structure_name",
            "histoseg_signed_distance_um",
        ]
    ].to_csv(args.output_dir / "histoseg_domains.csv", index=False)

    summary = summarize_domains(domains, modules)
    summary.to_csv(args.output_dir / "histoseg_domain_summary.csv", index=False)
    plot_domains(domains, args.output_dir)
    write_grid_preview(domains, args.output_dir)

    manifest = {
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "input_dir": str(args.input_dir),
        "coste_dir": str(args.coste_dir),
        "n_selected_genes": int(len(selected)),
        "n_gene_modules": int(args.n_gene_modules),
        "n_domains": int(args.n_domains),
        "n_cells": int(len(domains)),
        "n_neighbors": int(args.n_neighbors),
        "spatial_weight": float(args.spatial_weight),
        "seconds": float(time.perf_counter() - t0),
        "outputs": {
            "domains": "histoseg_domains.parquet",
            "domain_summary": "histoseg_domain_summary.csv",
            "gene_modules": "coste_gene_modules.csv",
            "module_scores": "coste_module_scores.parquet",
            "domain_figure": "coste_histoseg_domains.png",
        },
    }
    (args.output_dir / "histoseg_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build HistoSeg-style domains from Atera breast COSTE-selected genes.")
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--coste-dir", type=Path, default=DEFAULT_COSTE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-gene-modules", type=int, default=16)
    parser.add_argument("--n-domains", type=int, default=12)
    parser.add_argument("--n-neighbors", type=int, default=16)
    parser.add_argument("--spatial-weight", type=float, default=0.35)
    parser.add_argument("--random-state", type=int, default=7)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
