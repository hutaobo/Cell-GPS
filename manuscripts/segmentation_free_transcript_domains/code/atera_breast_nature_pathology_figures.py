#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import FancyArrowPatch, Rectangle
from skimage.transform import warp


DEFAULT_HE_IMAGE = Path(
    "/data/taobo.hu/pyxenium_lazyslide_breast_wta_20260507/data/"
    "WTA_Preview_FFPE_Breast_Cancer_he_image.ome.tif"
)
DEFAULT_HE_ALIGNMENT = Path(
    "/data/taobo.hu/pyxenium_lazyslide_breast_wta_20260507/data/"
    "WTA_Preview_FFPE_Breast_Cancer_he_alignment.csv"
)
DEFAULT_HE_KEYPOINTS = Path(
    "/data/taobo.hu/pyxenium_lazyslide_breast_wta_20260507/data/"
    "WTA_Preview_FFPE_Breast_Cancer_keypoints.csv"
)
DEFAULT_CELL_HISTOSEG = Path("/data/taobo.hu/atera_breast_coste1024_histoseg/coste1024_modules16_domains12")
DEFAULT_POINT_HISTOSEG = Path(
    "/data/taobo.hu/atera_breast_transcript_point_histoseg/coste1024_point_knn_train750k_assign_all"
)
DEFAULT_COMPARISON = Path("/data/taobo.hu/atera_breast_histoseg_comparison/point_vs_cell_coste1024")
DEFAULT_OUTPUT_DIR = Path("/data/taobo.hu/atera_breast_histoseg_nature_figures/pathology_explanations")

PIXEL_SIZE_UM = 0.2125
OVERVIEW_LEVEL = 5
OVERVIEW_DOWNSAMPLE = 32.0
ROI_LEVEL = 4
ROI_DOWNSAMPLE = 16.0


@dataclass
class ROI:
    xmin: float
    xmax: float
    ymin: float
    ymax: float


STORIES: dict[str, dict[str, Any]] = {
    "tumor_rich": {
        "title": "Tumor-rich molecular signal",
        "question": "Does transcript-point HistoSeg recover the tumor-rich signal?",
        "claim": "SERPINA6, MSMB, KCNQ3 and SERPINA1-rich transcript domains recover invasive tumor regions.",
        "answer": "Yes. The strongest point domain maps onto tumor-rich cell-domain context, while broader tumor signals split into transcript-level subdomains.",
        "answer_short": "P10 sits in the C10 tumor-rich region; full-data P10-C10 = 0.90.",
        "color": "#C2185B",
        "accent": "#7A003C",
        "point_domains": [1, 4, 7, 9, 10, 12],
        "local_examples": [
            {
                "label": "P10 over C10",
                "point_domain": 10,
                "cell_domain": 10,
                "width": 900,
                "height": 850,
                "genes": "SERPINA6, MSMB, KCNQ3, SERPINA1",
            }
        ],
        "analysis": [
            "The strongest one-to-one match is P10 to C10: 90.3% of P10 transcripts fall nearest to C10.",
            "The largest tumor-rich point domain P01 shares SERPINA6, KCNQ3, MSMB and SERPINA1 with C02, but spreads across several tumor-associated cell domains.",
            "Interpretation: the point-level map captures tumor-rich molecular signal and splits broad invasive regions into transcript-level subdomains.",
        ],
    },
    "immune_stroma": {
        "title": "Immune/stroma interface",
        "question": "Does transcript-point HistoSeg recover the immune-stromal interface?",
        "claim": "CCL22, DPT, IGF2, BMPER and RERGL-rich point domains identify immune-stromal interface regions.",
        "answer": "Yes. The point domains concentrate at immune-stromal interface regions and share immune-stromal marker logic with their cell-domain context.",
        "answer_short": "P06 localizes in the C08 immune-stromal interface; full-data P06-C08 = 0.63.",
        "color": "#008C8C",
        "accent": "#004F57",
        "point_domains": [2, 3, 5, 6, 8, 11],
        "local_examples": [
            {
                "label": "P06 over C08",
                "point_domain": 6,
                "cell_domain": 8,
                "width": 1450,
                "height": 1150,
                "genes": "CCL22, DPT, IGF2, BMPER, RERGL",
            }
        ],
        "analysis": [
            "P06 maps to C08 with 62.7% point-to-cell overlap and a 0.511 transcript-weighted Jaccard.",
            "Both sides carry immune-stromal marker genes: point P06 has RERGL, BMPER, IGF2, DPT and CCL22, while C08 has CCL22, SAA2, MS4A1, LTB and IGLC7.",
            "Interpretation: the point-level domain is not just immune or stromal; it marks an interface microenvironment.",
        ],
    },
    "mixed_luminal_apocrine": {
        "title": "Mixed luminal/apocrine microenvironment",
        "question": "Why do some immune-stromal-like point domains sit inside luminal/DCIS or apocrine contexts?",
        "claim": "Immune/stroma-like transcript domains appear inside luminal/DCIS and apocrine cell-domain contexts.",
        "answer": "They are best read as local microenvironment mixing inside epithelial pathological contexts, not as one pure luminal, apocrine, immune, or stromal zone.",
        "answer_short": "P02/P03 sit in C04/C03 epithelial contexts but keep immune-stromal marker logic, supporting local mixing.",
        "color": "#E68613",
        "accent": "#8F4A00",
        "point_domains": [2, 3, 11],
        "local_examples": [
            {
                "label": "P02 over luminal C04",
                "point_domain": 2,
                "cell_domain": 4,
                "width": 1200,
                "height": 950,
                "genes": "RERGL, BMPER, IGF2, DPT, CCL22",
            },
            {
                "label": "P03 over apocrine C03",
                "point_domain": 3,
                "cell_domain": 3,
                "width": 700,
                "height": 650,
                "genes": "RERGL, BMPER, IGF2, DPT, CCL22",
            },
        ],
        "analysis": [
            "P02 and P11 lie mainly in luminal/DCIS C04, while P03 almost entirely overlays apocrine C03.",
            "Their point-level marker logic is immune-stromal interface rather than pure luminal or pure apocrine.",
            "Interpretation: these regions likely contain local microenvironment mixing inside epithelial pathological contexts.",
        ],
    },
}


def set_theme() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.55,
            "ytick.major.width": 0.55,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def load_image_level(path: Path, level: int) -> np.ndarray:
    with tifffile.TiffFile(path) as tif:
        image = tif.series[0].levels[level].asarray()
    if image.ndim == 3 and image.shape[0] == 3:
        image = np.moveaxis(image, 0, -1)
    return image


def he_level_scales(path: Path) -> list[float]:
    with tifffile.TiffFile(path) as tif:
        levels = tif.series[0].levels
        shapes = []
        for level in levels:
            shape = level.shape
            axes = level.axes
            y_idx = axes.index("Y")
            x_idx = axes.index("X")
            shapes.append((int(shape[x_idx]), int(shape[y_idx])))
    return [shapes[0][0] / width for width, _height in shapes]


def load_alignment(path: Path) -> np.ndarray:
    matrix = np.loadtxt(path, delimiter=",")
    if matrix.shape != (3, 3):
        raise ValueError(f"Expected 3x3 HE alignment matrix, got {matrix.shape} from {path}")
    return matrix


def alignment_qc(path: Path, matrix: np.ndarray) -> dict[str, float]:
    if not path.exists():
        return {}
    keypoints = pd.read_csv(path)
    src = np.c_[keypoints[["alignmentX", "alignmentY"]].to_numpy(dtype=float), np.ones(len(keypoints))]
    pred = src @ matrix.T
    target = keypoints[["fixedX", "fixedY"]].to_numpy(dtype=float)
    errors = np.sqrt(((pred[:, :2] - target) ** 2).sum(axis=1))
    return {
        "n_keypoints": int(len(keypoints)),
        "mean_error_px": float(errors.mean()),
        "median_error_px": float(np.median(errors)),
        "max_error_px": float(errors.max()),
        "mean_error_um": float(errors.mean() * PIXEL_SIZE_UM),
        "max_error_um": float(errors.max() * PIXEL_SIZE_UM),
    }


def he_roi_patch(
    he_rgb: np.ndarray,
    he_level_scale: float,
    alignment_matrix: np.ndarray,
    roi: ROI,
    fixed_downsample: float,
    order: int = 1,
    background_value: float = 255.0,
) -> np.ndarray:
    x0 = roi.xmin / PIXEL_SIZE_UM
    x1 = roi.xmax / PIXEL_SIZE_UM
    y0 = roi.ymin / PIXEL_SIZE_UM
    y1 = roi.ymax / PIXEL_SIZE_UM
    width = max(1, int(math.ceil((x1 - x0) / fixed_downsample)))
    height = max(1, int(math.ceil((y1 - y0) / fixed_downsample)))
    level_scale_matrix = np.array(
        [[he_level_scale, 0.0, 0.0], [0.0, he_level_scale, 0.0], [0.0, 0.0, 1.0]],
        dtype=float,
    )
    transform_matrix = alignment_matrix @ level_scale_matrix
    inverse_matrix = np.linalg.inv(transform_matrix)

    def inverse_map(coords: np.ndarray) -> np.ndarray:
        flat = np.asarray(coords).reshape(-1, 2)
        fixed_x = x0 + flat[:, 0] * fixed_downsample
        fixed_y = y0 + flat[:, 1] * fixed_downsample
        hom = np.column_stack([fixed_x, fixed_y, np.ones(len(flat))])
        src = hom @ inverse_matrix.T
        return src[:, :2].reshape(np.asarray(coords).shape)

    channels = []
    for idx in range(he_rgb.shape[2]):
        warped = warp(
            he_rgb[..., idx],
            inverse_map=inverse_map,
            output_shape=(height, width),
            order=order,
            mode="constant",
            cval=background_value,
            preserve_range=True,
        )
        channels.append(warped.astype(np.uint8))
    return np.stack(channels, axis=-1)


def load_data(args: argparse.Namespace) -> dict[str, Any]:
    cells = pd.read_parquet(
        args.cell_histoseg_dir / "histoseg_domains.parquet",
        columns=["x", "y", "celltype", "histoseg_structure_id", "histoseg_structure_name"],
    )
    points = pd.read_parquet(args.comparison_dir / "point_preview_nearest_cell_domains.parquet")
    best = pd.read_csv(args.comparison_dir / "point_vs_cell_best_domain_matches.csv")
    overlap = pd.read_csv(args.comparison_dir / "point_vs_cell_overlap_point_fraction.csv", index_col=0)
    metrics = json.loads((args.comparison_dir / "point_vs_cell_global_metrics.json").read_text(encoding="utf-8"))
    overview = load_image_level(args.he_image, OVERVIEW_LEVEL)
    roi = load_image_level(args.he_image, ROI_LEVEL)
    scales = he_level_scales(args.he_image)
    alignment_matrix = load_alignment(args.he_alignment)
    return {
        "cells": cells,
        "points": points,
        "best": best,
        "overlap": overlap,
        "metrics": metrics,
        "overview_image": overview,
        "roi_image": roi,
        "he_scales": scales,
        "alignment_matrix": alignment_matrix,
        "alignment_qc": alignment_qc(args.he_keypoints, alignment_matrix),
        "he_patch_cache": {},
    }


def add_panel_label(ax: plt.Axes, label: str, color: str = "black") -> None:
    ax.text(
        -0.07,
        1.08,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.0,
        fontweight="bold",
        color=color,
    )


def add_scale_bar(ax: plt.Axes, roi: ROI, length_um: float, color: str = "black") -> None:
    width = roi.xmax - roi.xmin
    height = roi.ymax - roi.ymin
    x0 = roi.xmin + 0.06 * width
    x1 = x0 + length_um
    y = roi.ymax - 0.07 * height
    ax.plot([x0, x1], [y, y], color=color, lw=1.6, solid_capstyle="butt")
    ax.text((x0 + x1) / 2, y - 0.035 * height, f"{int(length_um)} µm", color=color, ha="center", va="bottom", fontsize=5.4)


def style_he_axis(ax: plt.Axes, roi: ROI, title: str = "", scale_um: float | None = None) -> None:
    ax.set_xlim(roi.xmin, roi.xmax)
    ax.set_ylim(roi.ymax, roi.ymin)
    ax.set_aspect("equal")
    ax.set_anchor("N")
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=6.5, pad=1.5)
    if scale_um is not None:
        add_scale_bar(ax, roi, scale_um, color="black")


def clamp_roi(center_x: float, center_y: float, width: float, height: float, bounds: ROI) -> ROI:
    xmin = max(bounds.xmin, center_x - width / 2)
    xmax = min(bounds.xmax, center_x + width / 2)
    ymin = max(bounds.ymin, center_y - height / 2)
    ymax = min(bounds.ymax, center_y + height / 2)
    if xmax - xmin < width and xmin <= bounds.xmin:
        xmax = min(bounds.xmax, xmin + width)
    if xmax - xmin < width and xmax >= bounds.xmax:
        xmin = max(bounds.xmin, xmax - width)
    if ymax - ymin < height and ymin <= bounds.ymin:
        ymax = min(bounds.ymax, ymin + height)
    if ymax - ymin < height and ymax >= bounds.ymax:
        ymin = max(bounds.ymin, ymax - height)
    return ROI(xmin, xmax, ymin, ymax)


def example_roi(points: pd.DataFrame, example: dict[str, Any], bounds: ROI) -> ROI:
    subset = points.loc[
        (points["histoseg_structure_id"] == int(example["point_domain"]))
        & (points["nearest_cell_domain_id"] == int(example["cell_domain"]))
    ]
    if len(subset) < 50:
        subset = points.loc[points["histoseg_structure_id"] == int(example["point_domain"])]
    center_x, center_y = dense_center(subset, float(example["width"]), float(example["height"]))
    return clamp_roi(center_x, center_y, float(example["width"]), float(example["height"]), bounds)


def dense_center(points: pd.DataFrame, width: float, height: float) -> tuple[float, float]:
    if len(points) < 100:
        return float(points["x"].median()), float(points["y"].median())
    bin_size = max(80.0, min(width, height) / 4.0)
    bins = pd.DataFrame(
        {
            "xbin": np.floor(points["x"].to_numpy(dtype=float) / bin_size).astype(np.int64),
            "ybin": np.floor(points["y"].to_numpy(dtype=float) / bin_size).astype(np.int64),
        },
        index=points.index,
    )
    counts = bins.value_counts(sort=True)
    xbin, ybin = counts.index[0]
    local = points.loc[(bins["xbin"] == xbin) & (bins["ybin"] == ybin)]
    if len(local) < 20:
        local = points
    return float(local["x"].median()), float(local["y"].median())


def cached_he_patch(data: dict[str, Any], roi: ROI, level: int, downsample: float) -> np.ndarray:
    key = (
        level,
        downsample,
        round(roi.xmin, 3),
        round(roi.xmax, 3),
        round(roi.ymin, 3),
        round(roi.ymax, 3),
    )
    cache = data["he_patch_cache"]
    if key not in cache:
        image = data["overview_image"] if level == OVERVIEW_LEVEL else data["roi_image"]
        cache[key] = he_roi_patch(
            image,
            he_level_scale=float(data["he_scales"][level]),
            alignment_matrix=data["alignment_matrix"],
            roi=roi,
            fixed_downsample=downsample,
        )
    return cache[key]


def draw_he(ax: plt.Axes, data: dict[str, Any], roi: ROI, level: int, downsample: float) -> None:
    patch = cached_he_patch(data, roi, level, downsample)
    ax.imshow(patch, extent=[roi.xmin, roi.xmax, roi.ymax, roi.ymin], origin="upper", interpolation="bilinear")
    style_he_axis(ax, roi)


def plot_global_overlay(
    ax: plt.Axes,
    data: dict[str, Any],
    story: dict[str, Any],
    rois: list[ROI],
    tissue_bounds: ROI,
) -> None:
    draw_he(ax, data, tissue_bounds, OVERVIEW_LEVEL, OVERVIEW_DOWNSAMPLE)
    points = data["points"]
    cells = data["cells"]
    point_domains = set(story["point_domains"])
    cell_domains = set(int(x) for x in data["best"].loc[data["best"]["point_domain_id"].isin(point_domains), "best_cell_domain_id"])
    point_subset = points.loc[points["histoseg_structure_id"].isin(point_domains)]
    if len(point_subset) > 55_000:
        point_subset = point_subset.sample(55_000, random_state=13)
    cell_subset = cells.loc[cells["histoseg_structure_id"].isin(cell_domains)]
    if len(cell_subset) > 22_000:
        cell_subset = cell_subset.sample(22_000, random_state=17)
    ax.scatter(cell_subset["x"], cell_subset["y"], s=0.9, c="#111111", alpha=0.18, linewidths=0, rasterized=True)
    ax.scatter(point_subset["x"], point_subset["y"], s=0.26, c=story["color"], alpha=0.46, linewidths=0, rasterized=True)
    for idx, roi in enumerate(rois, start=1):
        ax.add_patch(Rectangle((roi.xmin, roi.ymin), roi.xmax - roi.xmin, roi.ymax - roi.ymin, fill=False, ec=story["accent"], lw=0.85))
        ax.text(
            roi.xmin + 15,
            roi.ymin + 40,
            f"zoom {idx}",
            color=story["accent"],
            fontsize=5.5,
            weight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.4},
        )
    add_scale_bar(ax, tissue_bounds, 1000, color="black")
    ax.set_title("Global H&E overlap", fontsize=6.5, pad=1.5)


def plot_local_overlay(
    ax: plt.Axes,
    data: dict[str, Any],
    story: dict[str, Any],
    example: dict[str, Any],
    roi: ROI,
    *,
    show_title: bool = True,
) -> None:
    draw_he(ax, data, roi, ROI_LEVEL, ROI_DOWNSAMPLE)
    points = data["points"]
    cells = data["cells"]
    point_subset = points.loc[
        (points["histoseg_structure_id"] == int(example["point_domain"]))
        & (points["x"].between(roi.xmin, roi.xmax))
        & (points["y"].between(roi.ymin, roi.ymax))
    ]
    cell_subset = cells.loc[
        (cells["histoseg_structure_id"] == int(example["cell_domain"]))
        & (cells["x"].between(roi.xmin, roi.xmax))
        & (cells["y"].between(roi.ymin, roi.ymax))
    ]
    all_cells = cells.loc[
        (cells["x"].between(roi.xmin, roi.xmax))
        & (cells["y"].between(roi.ymin, roi.ymax))
    ]
    if len(all_cells):
        ax.scatter(all_cells["x"], all_cells["y"], s=0.9, c="#4D4D4D", alpha=0.12, linewidths=0, rasterized=True)
    ax.scatter(cell_subset["x"], cell_subset["y"], s=3.8, facecolors="none", edgecolors="#1A1A1A", linewidths=0.22, alpha=0.42, rasterized=True)
    ax.scatter(point_subset["x"], point_subset["y"], s=0.58, c=story["color"], alpha=0.72, linewidths=0, rasterized=True)
    add_scale_bar(ax, roi, 200 if roi.xmax - roi.xmin > 900 else 100, color="black")
    if show_title:
        ax.set_title(example["label"], fontsize=6.4, pad=1.5)


def plot_overlap_stats(
    ax: plt.Axes,
    data: dict[str, Any],
    story: dict[str, Any],
    *,
    compact: bool = False,
    show_title: bool = True,
    show_xlabel: bool = True,
    show_legend: bool = True,
) -> None:
    best = data["best"].loc[data["best"]["point_domain_id"].isin(story["point_domains"])].copy()
    best = best.sort_values("point_fraction_total", ascending=True)
    labels = [f"P{int(x):02d}-C{int(y):02d}" for x, y in zip(best["point_domain_id"], best["best_cell_domain_id"])]
    y = np.arange(len(best))
    ax.barh(y - 0.17, best["point_to_cell_fraction"], height=0.28, color=story["color"], alpha=0.92, label="point-to-cell fraction")
    ax.barh(y + 0.17, best["jaccard_transcript_weighted"], height=0.28, color="#BDBDBD", alpha=0.95, label="Jaccard")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5.6 if compact else 6)
    ax.set_xlim(0, 1.0)
    if show_xlabel:
        ax.set_xlabel("Overlap metric", fontsize=6.2 if compact else 7)
    else:
        ax.set_xlabel("")
        ax.set_xticklabels([])
    if show_title:
        ax.set_title("Full-data overlap statistics", fontsize=6.5, pad=1.5)
    if show_legend:
        if compact:
            ax.legend(loc="upper right", fontsize=5.2, handlelength=0.9)
        else:
            ax.legend(
                loc="lower center",
                bbox_to_anchor=(0.5, -0.25),
                ncol=2,
                fontsize=5.5,
                handlelength=0.9,
                columnspacing=0.9,
            )
    ax.grid(axis="x", lw=0.3, color="#D9D9D9")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="x", labelsize=5.5 if compact else 6)


def plot_marker_table(ax: plt.Axes, data: dict[str, Any], story: dict[str, Any]) -> None:
    ax.axis("off")
    best = data["best"].loc[data["best"]["point_domain_id"].isin(story["point_domains"])].copy()
    best = best.sort_values("point_to_cell_fraction", ascending=False).head(4)
    metrics = data["metrics"]
    lines = [
        story["claim"],
        "",
        f"Global comparison: ARI {metrics['adjusted_rand_index']:.2f}, NMI {metrics['normalized_mutual_info']:.2f},",
        f"point-to-cell purity {metrics['point_to_cell_purity']:.2f}.",
        "",
        "Key marker evidence:",
    ]
    for _, row in best.iterrows():
        genes = str(row["point_top_genes"]).replace("/", ", ")
        if len(genes) > 45:
            genes = genes[:42] + "..."
        lines.append(
            f"P{int(row['point_domain_id']):02d}-C{int(row['best_cell_domain_id']):02d}: "
            f"{row['point_to_cell_fraction']:.2f} overlap; {genes}"
        )
    lines.append("")
    lines.extend(story["analysis"])
    y = 0.98
    for idx, line in enumerate(lines):
        weight = "bold" if idx == 0 or line == "Key marker evidence:" else "normal"
        size = 6.8 if idx == 0 else 6.2
        ax.text(0.0, y, line, ha="left", va="top", fontsize=size, fontweight=weight, wrap=True)
        y -= 0.092 if line else 0.055


def plot_global_he_context(
    ax: plt.Axes,
    data: dict[str, Any],
    story: dict[str, Any],
    rois: list[ROI],
    tissue_bounds: ROI,
) -> None:
    draw_he(ax, data, tissue_bounds, OVERVIEW_LEVEL, OVERVIEW_DOWNSAMPLE)
    for idx, roi in enumerate(rois, start=1):
        ax.add_patch(
            Rectangle(
                (roi.xmin, roi.ymin),
                roi.xmax - roi.xmin,
                roi.ymax - roi.ymin,
                fill=False,
                ec=story["accent"],
                lw=0.9,
            )
        )
        ax.text(
            roi.xmin + 15,
            roi.ymin + 40,
            f"zoom {idx}",
            color=story["accent"],
            fontsize=5.5,
            weight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.4},
        )
    add_scale_bar(ax, tissue_bounds, 1000, color="black")
    ax.set_title("1. H&E context", fontsize=6.5, pad=1.5)


def domain_label(examples: list[dict[str, Any]], key: str, prefix: str) -> str:
    values = sorted({int(example[key]) for example in examples})
    return "/".join(f"{prefix}{value:02d}" for value in values)


def plot_example_point_global(
    ax: plt.Axes,
    data: dict[str, Any],
    story: dict[str, Any],
    examples: list[dict[str, Any]],
    rois: list[ROI],
    tissue_bounds: ROI,
) -> None:
    draw_he(ax, data, tissue_bounds, OVERVIEW_LEVEL, OVERVIEW_DOWNSAMPLE)
    points = data["points"]
    point_domains = {int(example["point_domain"]) for example in examples}
    point_subset = points.loc[points["histoseg_structure_id"].isin(point_domains)]
    if len(point_subset) > 45_000:
        point_subset = point_subset.sample(45_000, random_state=31)
    ax.scatter(point_subset["x"], point_subset["y"], s=0.34, c=story["color"], alpha=0.58, linewidths=0, rasterized=True)
    for idx, roi in enumerate(rois, start=1):
        ax.add_patch(Rectangle((roi.xmin, roi.ymin), roi.xmax - roi.xmin, roi.ymax - roi.ymin, fill=False, ec=story["accent"], lw=0.9))
        ax.text(
            roi.xmin + 15,
            roi.ymin + 40,
            f"zoom {idx}",
            color=story["accent"],
            fontsize=5.5,
            weight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.4},
        )
    add_scale_bar(ax, tissue_bounds, 1000, color="black")
    ax.set_title(f"1. {domain_label(examples, 'point_domain', 'P')} point global", fontsize=6.3, pad=1.5)


def plot_example_cell_global(
    ax: plt.Axes,
    data: dict[str, Any],
    story: dict[str, Any],
    examples: list[dict[str, Any]],
    rois: list[ROI],
    tissue_bounds: ROI,
) -> None:
    draw_he(ax, data, tissue_bounds, OVERVIEW_LEVEL, OVERVIEW_DOWNSAMPLE)
    cells = data["cells"]
    cell_domains = {int(example["cell_domain"]) for example in examples}
    cell_subset = cells.loc[cells["histoseg_structure_id"].isin(cell_domains)]
    if len(cell_subset) > 26_000:
        cell_subset = cell_subset.sample(26_000, random_state=37)
    ax.scatter(cell_subset["x"], cell_subset["y"], s=1.25, c=story["accent"], alpha=0.42, linewidths=0, rasterized=True)
    for idx, roi in enumerate(rois, start=1):
        ax.add_patch(Rectangle((roi.xmin, roi.ymin), roi.xmax - roi.xmin, roi.ymax - roi.ymin, fill=False, ec=story["accent"], lw=0.9))
        ax.text(
            roi.xmin + 15,
            roi.ymin + 40,
            f"zoom {idx}",
            color=story["accent"],
            fontsize=5.5,
            weight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.4},
        )
    add_scale_bar(ax, tissue_bounds, 1000, color="black")
    ax.set_title(f"2. {domain_label(examples, 'cell_domain', 'C')} cell global", fontsize=6.3, pad=1.5)


def plot_story_caption(ax: plt.Axes, story: dict[str, Any], step: str, title: str, body: str) -> None:
    ax.axis("off")
    ax.add_patch(
        Rectangle(
            (0.0, 0.0),
            1.0,
            1.0,
            transform=ax.transAxes,
            fc="#F7F7F7",
            ec="#D0D0D0",
            lw=0.35,
            clip_on=False,
        )
    )
    ax.text(
        0.035,
        0.82,
        step,
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=6.2,
        fontweight="bold",
        color="white",
        bbox={"boxstyle": "round,pad=0.18", "facecolor": story["color"], "edgecolor": "none"},
    )
    ax.text(0.23, 0.86, title, transform=ax.transAxes, ha="left", va="top", fontsize=5.7, fontweight="bold")
    ax.text(
        0.035,
        0.55,
        textwrap.fill(body, width=24),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.0,
        linespacing=1.15,
    )


def plot_story_answer(ax: plt.Axes, data: dict[str, Any], story: dict[str, Any]) -> None:
    ax.axis("off")
    best = data["best"].loc[data["best"]["point_domain_id"].isin(story["point_domains"])].copy()
    best = best.sort_values("point_to_cell_fraction", ascending=False).head(3)
    support = "; ".join(
        [
            f"P{int(row.point_domain_id):02d}-C{int(row.best_cell_domain_id):02d} {row.point_to_cell_fraction:.2f}"
            for row in best.itertuples(index=False)
        ]
    )
    ax.add_patch(
        Rectangle(
            (0.0, 0.0),
            1.0,
            1.0,
            transform=ax.transAxes,
            fc="#FFFFFF",
            ec=story["accent"],
            lw=0.8,
            clip_on=False,
        )
    )
    ax.text(0.04, 0.94, "5. Pathology answer", transform=ax.transAxes, ha="left", va="top", fontsize=6.4, fontweight="bold")
    ax.text(
        0.04,
        0.78,
        textwrap.fill(story.get("answer_short", story["answer"]), width=25),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.1,
        linespacing=1.18,
    )
    ax.text(0.04, 0.36, "Full-data support", transform=ax.transAxes, ha="left", va="top", fontsize=5.7, fontweight="bold")
    ax.text(
        0.04,
        0.25,
        textwrap.fill(support, width=24),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.0,
        linespacing=1.15,
    )


def add_local_example_label(ax: plt.Axes, roi: ROI, label: str) -> None:
    width = roi.xmax - roi.xmin
    height = roi.ymax - roi.ymin
    ax.text(
        roi.xmin + 0.035 * width,
        roi.ymin + 0.07 * height,
        label,
        ha="left",
        va="top",
        fontsize=5.3,
        color="black",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.6},
    )


def plot_story_support_stats(ax: plt.Axes, data: dict[str, Any], story: dict[str, Any]) -> None:
    best = data["best"].loc[data["best"]["point_domain_id"].isin(story["point_domains"])].copy()
    best = best.sort_values("point_to_cell_fraction", ascending=False).head(3)
    best = best.sort_values("point_to_cell_fraction", ascending=True)
    labels = [f"P{int(x):02d}-C{int(y):02d}" for x, y in zip(best["point_domain_id"], best["best_cell_domain_id"])]
    y = np.arange(len(best))
    ax.barh(y - 0.15, best["point_to_cell_fraction"], height=0.26, color=story["color"], alpha=0.92)
    ax.barh(y + 0.15, best["jaccard_transcript_weighted"], height=0.26, color="#BDBDBD", alpha=0.95)
    for yi, label, row in zip(y, labels, best.itertuples(index=False)):
        ax.text(0.02, yi, label, va="center", ha="left", fontsize=5.3)
        ax.text(min(float(row.point_to_cell_fraction) + 0.03, 0.96), yi - 0.15, f"{row.point_to_cell_fraction:.2f}", va="center", fontsize=4.8)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Overlap metric", fontsize=5.9)
    ax.set_title("4. Top full-data support", fontsize=6.4, pad=1.5)
    ax.grid(axis="x", lw=0.3, color="#D9D9D9")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="x", labelsize=5.2)
    ax.tick_params(axis="y", length=0)


def add_flow_arrows(fig: plt.Figure, axes: list[plt.Axes], color: str) -> None:
    fig.canvas.draw()
    for left, right in zip(axes[:-1], axes[1:]):
        left_box = left.get_position()
        right_box = right.get_position()
        x0 = left_box.x1 + 0.006
        x1 = right_box.x0 - 0.006
        if x1 <= x0:
            continue
        y = (left_box.y0 + left_box.y1) / 2.0
        fig.add_artist(
            FancyArrowPatch(
                (x0, y),
                (x1, y),
                transform=fig.transFigure,
                arrowstyle="-|>",
                mutation_scale=8,
                lw=0.7,
                color=color,
                alpha=0.85,
                clip_on=False,
            )
        )


def save_pub(fig: plt.Figure, output_dir: Path, stem: str, *, exact_size: bool = False) -> list[str]:
    paths: list[str] = []
    for ext, kwargs in [
        ("png", {"dpi": 450}),
        ("svg", {}),
        ("pdf", {}),
        ("tiff", {"dpi": 600}),
    ]:
        path = output_dir / f"{stem}.{ext}"
        if exact_size:
            fig.savefig(path, **kwargs)
        else:
            fig.savefig(path, bbox_inches="tight", **kwargs)
        paths.append(str(path))
    plt.close(fig)
    return paths


def make_story_plate(data: dict[str, Any], story_key: str, story: dict[str, Any], output_dir: Path, tissue_bounds: ROI) -> dict[str, Any]:
    rois = [example_roi(data["points"], example, tissue_bounds) for example in story["local_examples"]]
    fig = plt.figure(figsize=(7.2, 5.25))
    grid = GridSpec(2, 3, figure=fig, width_ratios=[1.2, 1.12, 1.02], height_ratios=[1.02, 1.0], hspace=0.28, wspace=0.34)

    ax_a = fig.add_subplot(grid[:, 0])
    plot_global_overlay(ax_a, data, story, rois, tissue_bounds)
    add_panel_label(ax_a, "a")

    local_grid = GridSpecFromSubplotSpec(len(rois), 1, subplot_spec=grid[0, 1], hspace=0.15)
    local_axes = []
    for idx, (example, roi) in enumerate(zip(story["local_examples"], rois)):
        ax = fig.add_subplot(local_grid[idx, 0])
        plot_local_overlay(ax, data, story, example, roi)
        local_axes.append(ax)
    add_panel_label(local_axes[0], "b")

    ax_c = fig.add_subplot(grid[1, 1])
    plot_overlap_stats(ax_c, data, story)
    add_panel_label(ax_c, "c")

    ax_d = fig.add_subplot(grid[:, 2])
    plot_marker_table(ax_d, data, story)
    add_panel_label(ax_d, "d")

    fig.suptitle(story["title"], x=0.02, y=0.995, ha="left", fontsize=8.2, fontweight="bold")
    paths = save_pub(fig, output_dir, f"nature_pathology_{story_key}_plate")
    return {
        "story": story_key,
        "title": story["title"],
        "claim": story["claim"],
        "local_rois": [roi.__dict__ for roi in rois],
        "outputs": paths,
    }


def make_story_flow(data: dict[str, Any], story_key: str, story: dict[str, Any], output_dir: Path, tissue_bounds: ROI) -> dict[str, Any]:
    rois = [example_roi(data["points"], example, tissue_bounds) for example in story["local_examples"]]
    examples = story["local_examples"]
    fig = plt.figure(figsize=(183 / 25.4, 92 / 25.4))
    grid = GridSpec(
        2,
        5,
        figure=fig,
        height_ratios=[0.62, 0.38],
        width_ratios=[1.05, 1.05, 1.05, 0.96, 1.12],
        hspace=0.1,
        wspace=0.36,
    )
    fig.subplots_adjust(left=0.035, right=0.985, top=0.84, bottom=0.075)

    ax_point_global = fig.add_subplot(grid[0, 0])
    plot_example_point_global(ax_point_global, data, story, examples, rois, tissue_bounds)

    ax_cell_global = fig.add_subplot(grid[0, 1])
    plot_example_cell_global(ax_cell_global, data, story, examples, rois, tissue_bounds)

    local_grid = GridSpecFromSubplotSpec(len(rois), 1, subplot_spec=grid[0, 2], hspace=0.16)
    local_axes: list[plt.Axes] = []
    for idx, (example, roi) in enumerate(zip(examples, rois)):
        ax = fig.add_subplot(local_grid[idx, 0])
        plot_local_overlay(ax, data, story, example, roi, show_title=False)
        if idx == 0:
            ax.set_title("3. Local overlap", fontsize=6.2, pad=1.4)
        add_local_example_label(ax, roi, example["label"])
        local_axes.append(ax)

    ax_stats = fig.add_subplot(grid[0, 3])
    plot_story_support_stats(ax_stats, data, story)

    ax_answer = fig.add_subplot(grid[0, 4])
    plot_story_answer(ax_answer, data, story)

    genes = story["local_examples"][0]["genes"]
    point_domains = domain_label(examples, "point_domain", "P")
    cell_domains = domain_label(examples, "cell_domain", "C")
    captions = [
        ("1", "Point global", f"Show {point_domains} alone on aligned H&E."),
        ("2", "Cell global", f"Show matched {cell_domains} alone on the same H&E."),
        (
            "3",
            "Inspect overlap",
            f"Zoom into the densest {point_domains}/{cell_domains} overlap." if len(rois) == 1 else "Inspect each point/cell pair locally.",
        ),
        ("4", "Check all points", "Full 60.9M-point comparison; color is point fraction, gray is Jaccard."),
        ("5", "Answer question", "Combine marker logic, H&E context, and statistics."),
    ]
    for col, (step, title, body) in enumerate(captions):
        ax_caption = fig.add_subplot(grid[1, col])
        plot_story_caption(ax_caption, story, step, title, body)

    fig.suptitle(story["question"], x=0.02, y=0.975, ha="left", fontsize=7.0, fontweight="bold")
    add_flow_arrows(fig, [ax_point_global, ax_cell_global, local_axes[0], ax_stats, ax_answer], story["accent"])
    paths = save_pub(fig, output_dir, f"nature_pathology_story_{story_key}_steps", exact_size=True)
    return {
        "story": story_key,
        "title": story["title"],
        "question": story["question"],
        "claim": story["claim"],
        "local_rois": [roi.__dict__ for roi in rois],
        "outputs": paths,
    }


def make_triptych(data: dict[str, Any], output_dir: Path, tissue_bounds: ROI) -> dict[str, Any]:
    fig = plt.figure(figsize=(7.2, 7.0))
    grid = GridSpec(4, 3, figure=fig, height_ratios=[0.95, 0.95, 0.95, 0.62], hspace=0.34, wspace=0.3)
    row_info: list[dict[str, Any]] = []
    for row_idx, (story_key, story) in enumerate(STORIES.items()):
        example = story["local_examples"][0]
        roi = example_roi(data["points"], example, tissue_bounds)
        row_info.append({"story": story_key, "roi": roi.__dict__})
        ax_global = fig.add_subplot(grid[row_idx, 0])
        plot_global_overlay(ax_global, data, story, [roi], tissue_bounds)
        if row_idx == 0:
            add_panel_label(ax_global, "a")
        ax_global.set_title(story["title"], fontsize=6.5, pad=1.5)
        ax_local = fig.add_subplot(grid[row_idx, 1])
        plot_local_overlay(ax_local, data, story, example, roi)
        if row_idx == 0:
            add_panel_label(ax_local, "b")
        ax_stats = fig.add_subplot(grid[row_idx, 2])
        plot_overlap_stats(
            ax_stats,
            data,
            story,
            compact=True,
            show_title=row_idx == 0,
            show_xlabel=row_idx == 2,
            show_legend=row_idx == 0,
        )
        if row_idx == 0:
            add_panel_label(ax_stats, "c")

    ax_text = fig.add_subplot(grid[3, :])
    ax_text.axis("off")
    notes = [
        "Interpretation summary",
        "1. Tumor-rich domains recover SERPINA6, MSMB, KCNQ3 and SERPINA1-positive invasive tumor regions.",
        "2. Immune/stroma domains recover CCL22, DPT, IGF2, BMPER and RERGL interface regions.",
        "3. Immune/stroma-like domains inside luminal/DCIS or apocrine contexts suggest local microenvironment mixing.",
        "All overlap statistics use 60,906,532 selected-gene transcript points; HE panels use a 500,000-point visual preview.",
    ]
    y = 0.98
    for idx, note in enumerate(notes):
        ax_text.text(0.0, y, note, fontsize=6.0 if idx else 6.5, fontweight="bold" if idx == 0 else "normal", va="top")
        y -= 0.17
    add_panel_label(ax_text, "d")
    paths = save_pub(fig, output_dir, "nature_pathology_triptych_summary")
    return {"outputs": paths, "rows": row_info}


def write_contract(output_dir: Path) -> None:
    contract = {
        "core_conclusion": (
            "Transcript-point HistoSeg recovers tumor-rich, immune-stromal interface, and mixed epithelial "
            "microenvironment domains that are spatially related to but finer than cell-based COSTE-HistoSeg."
        ),
        "figure_archetype": "image plate + quant",
        "backend": "Python/matplotlib",
        "target": "Nature-style imaging plus quantitative validation",
        "story_flow_export": {
            "width_mm": 183,
            "height_mm": 92,
            "png_dpi": 450,
            "tiff_dpi": 600,
            "font_range_pt": "5-7 for story-flow text",
            "layout": "double-column width, exact-size export without tight bounding-box resizing",
        },
        "panel_logic": {
            "story_flow": "turns each pathology question into a five-step point-domain/global cell-domain/global visual argument",
            "global_HE_overlap": "shows tissue-wide localization of the point-domain signature on histology",
            "local_HE_overlap": "shows representative microscopic colocalization with cell-domain context",
            "statistics": "uses full 60.9M transcript-point overlap metrics",
            "analysis": "connects overlap, dominant cell type, and COSTE marker genes to pathology interpretation",
        },
        "image_integrity": (
            "HE image is shown as an aligned spatial context; scatter overlays are transparent and use "
            "a fixed preview sample for readability. Quantitative bars use full transcript-point counts."
        ),
    }
    (output_dir / "nature_pathology_figure_contract.json").write_text(json.dumps(contract, indent=2), encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    set_theme()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = load_data(args)
    cells = data["cells"]
    tissue_bounds = ROI(
        xmin=0.0,
        xmax=float(max(cells["x"].max(), data["points"]["x"].max())) + 150,
        ymin=0.0,
        ymax=float(max(cells["y"].max(), data["points"]["y"].max())) + 150,
    )
    outputs: list[dict[str, Any]] = []
    story_flows: list[dict[str, Any]] = []
    for story_key, story in STORIES.items():
        outputs.append(make_story_plate(data, story_key, story, args.output_dir, tissue_bounds))
        story_flows.append(make_story_flow(data, story_key, story, args.output_dir, tissue_bounds))
    triptych = make_triptych(data, args.output_dir, tissue_bounds)
    write_contract(args.output_dir)

    source_rows: list[dict[str, Any]] = []
    for story_key, story in STORIES.items():
        story_best = data["best"].loc[data["best"]["point_domain_id"].isin(story["point_domains"])].copy()
        story_best.insert(0, "story", story_key)
        source_rows.extend(story_best.to_dict(orient="records"))
    pd.DataFrame(source_rows).to_csv(args.output_dir / "nature_pathology_source_domain_statistics.csv", index=False)
    manifest = {
        "created_at": pd.Timestamp.utcnow().isoformat(),
        "workflow": "nature_style_pathology_explanation_figures",
        "he_image": str(args.he_image),
        "he_alignment": str(args.he_alignment),
        "he_keypoints": str(args.he_keypoints),
        "alignment_qc": data["alignment_qc"],
        "cell_histoseg_dir": str(args.cell_histoseg_dir),
        "point_histoseg_dir": str(args.point_histoseg_dir),
        "comparison_dir": str(args.comparison_dir),
        "uses_spatial_grid": False,
        "statistics_unit": "all selected-gene transcript points",
        "visual_overlay_unit": "500000 transcript-point preview plus cell centroids",
        "stories": outputs,
        "story_flows": story_flows,
        "triptych": triptych,
        "outputs": {
            "contract": "nature_pathology_figure_contract.json",
            "source_statistics": "nature_pathology_source_domain_statistics.csv",
        },
    }
    (args.output_dir / "nature_pathology_figure_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Nature-style pathology explanation figures for transcript-point HistoSeg.")
    parser.add_argument("--he-image", type=Path, default=DEFAULT_HE_IMAGE)
    parser.add_argument("--he-alignment", type=Path, default=DEFAULT_HE_ALIGNMENT)
    parser.add_argument("--he-keypoints", type=Path, default=DEFAULT_HE_KEYPOINTS)
    parser.add_argument("--cell-histoseg-dir", type=Path, default=DEFAULT_CELL_HISTOSEG)
    parser.add_argument("--point-histoseg-dir", type=Path, default=DEFAULT_POINT_HISTOSEG)
    parser.add_argument("--comparison-dir", type=Path, default=DEFAULT_COMPARISON)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
