from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_SOURCE_DIR = Path("/data/taobo.hu/atera_urna_segmentation_bias_20260604")
DEFAULT_CELLS = Path("/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast/cells.csv.gz")
DEFAULT_CELL_GROUPS = Path(
    "/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast/WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv"
)

DEFAULT_GENES = ["CCND1", "C1QA", "JCHAIN", "CDH5", "ERBB2"]
TARGET_OVERRIDES = {
    "CCND1": "11q13 Invasive Tumor Cells",
    "C1QA": "Macrophages",
    "JCHAIN": "Plasma Cells",
    "CDH5": "Endothelial Cells",
    "ERBB2": "CXCL14+ Fibroblasts",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract real spatial windows for Atera uRNA/COSTE examples.")
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--cells", type=Path, default=DEFAULT_CELLS)
    parser.add_argument("--cell-groups", type=Path, default=DEFAULT_CELL_GROUPS)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--genes", nargs="*", default=DEFAULT_GENES)
    parser.add_argument("--window-um", type=float, default=850.0)
    parser.add_argument("--grid-um", type=float, default=120.0)
    parser.add_argument("--max-urna-points", type=int, default=7000)
    parser.add_argument("--max-other-cells", type=int, default=3500)
    parser.add_argument("--seed", type=int, default=17)
    return parser.parse_args()


def load_cells(cells_path: Path, cell_groups_path: Path) -> pd.DataFrame:
    cells = pd.read_csv(cells_path)
    groups = pd.read_csv(cell_groups_path)
    groups = groups.rename(columns={"Barcode": "cell_id", "Clusters": "celltype", "group": "celltype"})
    cells = cells.rename(columns={"x_centroid": "x", "y_centroid": "y"})
    merged = cells[["cell_id", "x", "y"]].merge(groups[["cell_id", "celltype"]], on="cell_id", how="left")
    merged = merged.dropna(subset=["x", "y", "celltype"]).copy()
    merged["celltype"] = merged["celltype"].astype(str)
    merged = merged.loc[merged["celltype"] != "Unassigned"].reset_index(drop=True)
    return merged


def load_urna_coords(source_dir: Path, genes: list[str]) -> pd.DataFrame:
    coords_path = source_dir / "selected_urna_coords_530genes.parquet"
    if not coords_path.exists():
        coords_path = source_dir / "selected_urna_coords_39genes.parquet"
    coords = pd.read_parquet(coords_path)
    coords = coords.rename(columns={"feature_name": "gene"})
    coords = coords.loc[coords["gene"].isin(genes), ["gene", "x", "y"]].copy()
    return coords


def choose_window(
    gene_coords: pd.DataFrame,
    target_cells: pd.DataFrame,
    window_um: float,
    grid_um: float,
) -> tuple[float, float]:
    if gene_coords.empty:
        raise ValueError("No uRNA coordinates available for gene")
    combined = pd.concat(
        [
            gene_coords[["x", "y"]].assign(weight=1.0),
            target_cells[["x", "y"]].assign(weight=1.8),
        ],
        ignore_index=True,
    )
    min_x, max_x = combined["x"].min(), combined["x"].max()
    min_y, max_y = combined["y"].min(), combined["y"].max()
    x_bins = np.arange(min_x, max_x + grid_um, grid_um)
    y_bins = np.arange(min_y, max_y + grid_um, grid_um)

    gene_hist, x_edges, y_edges = np.histogram2d(gene_coords["x"], gene_coords["y"], bins=[x_bins, y_bins])
    target_hist, _, _ = np.histogram2d(target_cells["x"], target_cells["y"], bins=[x_bins, y_bins])
    if target_hist.max() == 0:
        target_hist = np.ones_like(gene_hist)

    gene_score = gene_hist / max(gene_hist.max(), 1)
    target_score = target_hist / max(target_hist.max(), 1)
    score = gene_score + target_score + 2.0 * np.sqrt(gene_score * target_score)
    best_i, best_j = np.unravel_index(np.nanargmax(score), score.shape)
    center_x = float((x_edges[best_i] + x_edges[best_i + 1]) / 2)
    center_y = float((y_edges[best_j] + y_edges[best_j + 1]) / 2)

    # Refine once by maximizing the exact number of uRNAs and target cells in nearby candidate centers.
    candidates = []
    for di in range(-2, 3):
        for dj in range(-2, 3):
            ii = min(max(best_i + di, 0), len(x_edges) - 2)
            jj = min(max(best_j + dj, 0), len(y_edges) - 2)
            candidates.append(((x_edges[ii] + x_edges[ii + 1]) / 2, (y_edges[jj] + y_edges[jj + 1]) / 2))

    best = (center_x, center_y)
    best_score = -np.inf
    half = window_um / 2
    for cx, cy in candidates:
        g = in_window(gene_coords, cx, cy, half).sum()
        t = in_window(target_cells, cx, cy, half).sum()
        exact_score = g + 20 * t + 2 * np.sqrt(max(g, 0) * max(t, 0))
        if exact_score > best_score:
            best_score = exact_score
            best = (float(cx), float(cy))
    return best


def in_window(df: pd.DataFrame, cx: float, cy: float, half: float) -> pd.Series:
    return (df["x"] >= cx - half) & (df["x"] <= cx + half) & (df["y"] >= cy - half) & (df["y"] <= cy + half)


def downsample(df: pd.DataFrame, max_rows: int, rng: np.random.Generator) -> pd.DataFrame:
    if len(df) <= max_rows:
        return df
    return df.iloc[rng.choice(len(df), size=max_rows, replace=False)].copy()


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cells = load_cells(args.cells, args.cell_groups)
    coords = load_urna_coords(args.source_dir, args.genes)
    summary = pd.read_csv(args.source_dir / "urna_segmentation_bias_gene_summary.csv")

    point_rows = []
    meta_rows = []
    half = args.window_um / 2

    for gene in args.genes:
        gene_coords = coords.loc[coords["gene"] == gene].copy()
        if gene_coords.empty:
            continue
        summary_row = summary.loc[summary["gene"] == gene].iloc[0].to_dict()
        target = TARGET_OVERRIDES.get(gene, summary_row["urna_best_cluster"])
        target_cells = cells.loc[cells["celltype"] == target].copy()
        cx, cy = choose_window(gene_coords, target_cells, args.window_um, args.grid_um)

        gene_window = gene_coords.loc[in_window(gene_coords, cx, cy, half)].copy()
        cell_window = cells.loc[in_window(cells, cx, cy, half)].copy()
        target_window = cell_window.loc[cell_window["celltype"] == target].copy()
        other_window = cell_window.loc[cell_window["celltype"] != target].copy()

        gene_window = downsample(gene_window, args.max_urna_points, rng)
        other_window = downsample(other_window, args.max_other_cells, rng)

        gene_window["example_gene"] = gene
        gene_window["point_type"] = "uRNA"
        gene_window["celltype"] = ""
        gene_window["target_celltype"] = target
        gene_window["center_x"] = cx
        gene_window["center_y"] = cy
        point_rows.append(gene_window[["example_gene", "point_type", "celltype", "target_celltype", "x", "y", "center_x", "center_y"]])

        target_window["example_gene"] = gene
        target_window["point_type"] = "target_cell"
        target_window["target_celltype"] = target
        target_window["center_x"] = cx
        target_window["center_y"] = cy
        point_rows.append(
            target_window[["example_gene", "point_type", "celltype", "target_celltype", "x", "y", "center_x", "center_y"]]
        )

        other_window["example_gene"] = gene
        other_window["point_type"] = "other_cell"
        other_window["target_celltype"] = target
        other_window["center_x"] = cx
        other_window["center_y"] = cy
        point_rows.append(
            other_window[["example_gene", "point_type", "celltype", "target_celltype", "x", "y", "center_x", "center_y"]]
        )

        meta_rows.append(
            {
                "gene": gene,
                "target_celltype": target,
                "center_x": cx,
                "center_y": cy,
                "window_um": args.window_um,
                "window_urna_points_plotted": int(len(gene_window)),
                "window_target_cells": int(len(target_window)),
                "window_other_cells_plotted": int(len(other_window)),
                "total_gene_urna": int(summary_row["unassigned_count"]),
                "urna_fraction": float(summary_row["unassigned_fraction"]),
                "urna_best_cluster": summary_row["urna_best_cluster"],
                "urna_best_sss": float(summary_row["urna_best_sss"]),
                "full_best_cluster": summary_row["full_best_cluster"],
                "sss_spearman_vs_full": float(summary_row["sss_spearman_vs_full"]),
                "best_cluster_match": bool(summary_row["best_cluster_match"]),
            }
        )

    points = pd.concat(point_rows, ignore_index=True)
    metadata = pd.DataFrame(meta_rows)
    points.to_csv(args.output_dir / "spatial_example_points.csv", index=False)
    metadata.to_csv(args.output_dir / "spatial_example_metadata.csv", index=False)
    (args.output_dir / "spatial_example_config.json").write_text(
        json.dumps(
            {
                "source_dir": str(args.source_dir),
                "cells": str(args.cells),
                "cell_groups": str(args.cell_groups),
                "genes": args.genes,
                "window_um": args.window_um,
                "grid_um": args.grid_um,
                "max_urna_points": args.max_urna_points,
                "max_other_cells": args.max_other_cells,
                "seed": args.seed,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
