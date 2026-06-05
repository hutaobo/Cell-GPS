from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D


HERE = Path(__file__).resolve().parent
RESULTS = Path(__file__).resolve().parents[2] / "benchmarking" / "atera_wta_breast_urna_segmentation_bias_results"
SPATIAL = HERE / "spatial_examples"
FIGURES = HERE / "figures"
TABLES = HERE / "tables"

EXAMPLE_GENES = ["CCND1", "C1QA", "JCHAIN", "CDH5", "ERBB2"]
REPRESENTATIVE_CELLTYPES = [
    "11q13 Invasive Tumor Cells",
    "Basal-like Structured DCIS Cells",
    "CXCL14+ Fibroblasts",
    "Macrophages",
    "Plasma Cells",
    "Endothelial Cells",
]
CELLTYPE_LABELS = {
    "11q13 Invasive Tumor Cells": "11q13\nTumor",
    "Basal-like Structured DCIS Cells": "Basal-like\nDCIS",
    "CXCL14+ Fibroblasts": "CXCL14+\nFib.",
    "Macrophages": "Macro.",
    "Plasma Cells": "Plasma",
    "Endothelial Cells": "Endoth.",
}
GENE_LABELS = {
    "CCND1": "CCND1\n11q13 tumor",
    "C1QA": "C1QA\nmacrophage",
    "JCHAIN": "JCHAIN\nplasma",
    "CDH5": "CDH5\nendothelial",
    "ERBB2": "ERBB2\ndiscordant",
}
TARGET_COLORS = {
    "CCND1": "#b23a48",
    "C1QA": "#4c956c",
    "JCHAIN": "#6f5aa7",
    "CDH5": "#178f8f",
    "ERBB2": "#c48a22",
}
INK = "#1f2933"
MUTED = "#6b7886"
LIGHT = "#d8dde5"


def load_inputs():
    points = pd.read_csv(SPATIAL / "spatial_example_points.csv")
    metadata = pd.read_csv(SPATIAL / "spatial_example_metadata.csv")
    summary = pd.read_csv(RESULTS / "urna_segmentation_bias_gene_summary.csv")
    sss = pd.read_csv(RESULTS / "urna_only_t_and_c_result.csv", index_col=0)
    marker_summary = pd.read_csv(TABLES / "marker_group_summary.csv")
    return points, metadata, summary, sss, marker_summary


def setup_matplotlib():
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel_label(ax, label: str, x: float = -0.12, y: float = 1.08):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=9, fontweight="bold", color=INK, va="top", ha="left")


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#b9c2cd")
    ax.spines["bottom"].set_color("#b9c2cd")
    ax.tick_params(colors=INK, labelsize=6, length=2.5, pad=1)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def plot_spatial(ax, points: pd.DataFrame, metadata: pd.Series, gene: str):
    sub = points.loc[points["example_gene"] == gene].copy()
    u = sub.loc[sub["point_type"] == "uRNA"]
    target = sub.loc[sub["point_type"] == "target_cell"]
    other = sub.loc[sub["point_type"] == "other_cell"]
    color = TARGET_COLORS[gene]
    cx = float(metadata["center_x"])
    cy = float(metadata["center_y"])
    half = float(metadata["window_um"]) / 2

    ax.scatter(other["x"], other["y"], s=3.0, color="#c7cdd4", alpha=0.45, linewidths=0, rasterized=True)
    ax.scatter(target["x"], target["y"], s=10.0, facecolors="none", edgecolors=color, alpha=0.80, linewidths=0.55, rasterized=True)
    ax.scatter(u["x"], u["y"], s=2.2, color="#111827", alpha=0.33, linewidths=0, rasterized=True)

    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#a6b3c1")
        spine.set_linewidth(0.7)

    ax.text(0.03, 0.96, gene, transform=ax.transAxes, fontsize=8, fontweight="bold", color=INK, ha="left", va="top")
    ax.text(
        0.03,
        0.86,
        short_cluster(str(metadata["target_celltype"])),
        transform=ax.transAxes,
        fontsize=6.2,
        color=color,
        ha="left",
        va="top",
        linespacing=1.0,
    )

    scale_um = 200
    x0 = cx + half - 250
    y0 = cy + half - 65
    ax.plot([x0, x0 + scale_um], [y0, y0], color=INK, lw=1.1)
    ax.text(x0 + scale_um / 2, y0 - 24, "200 um", fontsize=5.5, ha="center", va="top", color=INK)


def plot_sss_heatmap(ax, sss: pd.DataFrame, metadata: pd.DataFrame):
    heat = 1 - sss.loc[EXAMPLE_GENES, REPRESENTATIVE_CELLTYPES].astype(float).clip(0, 1)
    image = ax.imshow(heat.to_numpy(), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(REPRESENTATIVE_CELLTYPES)))
    ax.set_xticklabels([CELLTYPE_LABELS[c] for c in REPRESENTATIVE_CELLTYPES], fontsize=6)
    ax.set_yticks(range(len(EXAMPLE_GENES)))
    ax.set_yticklabels([GENE_LABELS[g] for g in EXAMPLE_GENES], fontsize=6)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for i, gene in enumerate(EXAMPLE_GENES):
        target = metadata.set_index("gene").loc[gene, "target_celltype"]
        if target in REPRESENTATIVE_CELLTYPES:
            j = REPRESENTATIVE_CELLTYPES.index(target)
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#ffffff", lw=1.4))
    ax.set_title("uRNA-only COSTE proximity (1 - SSS; higher is closer)", fontsize=8, color=INK, pad=5)
    return image


def plot_concordance_context(ax, summary: pd.DataFrame, metadata: pd.DataFrame):
    vals = summary["sss_spearman_vs_full"].dropna()
    ax.hist(vals, bins=np.linspace(-0.35, 1.02, 32), color="#d7dde5", edgecolor="white", linewidth=0.4)
    ax.axvline(vals.median(), color=INK, lw=1.0)
    ymax = ax.get_ylim()[1]
    meta = metadata.set_index("gene")
    matched_rhos = []
    for gene in EXAMPLE_GENES:
        rho = float(meta.loc[gene, "sss_spearman_vs_full"])
        ax.plot([rho, rho], [0, ymax * 0.18], color=TARGET_COLORS[gene], lw=1.4)
        if gene == "ERBB2":
            ax.text(rho, ymax * 0.24, gene, fontsize=5.8, color=TARGET_COLORS[gene], ha="center")
        else:
            matched_rhos.append(rho)
    if matched_rhos:
        ax.text(0.98, ymax * 0.24, "matched\nexamples", fontsize=5.8, color=MUTED, ha="right", va="top", linespacing=0.9)
    ax.set_xlabel("rho: uRNA-only vs all-transcript SSS")
    ax.set_ylabel("genes")
    ax.set_title("Global concordance context", fontsize=8, color=INK, pad=5)
    style_axis(ax)


def plot_marker_summary(ax, marker_summary: pd.DataFrame):
    order = ["11q13", "Macrophage", "Plasma", "Vascular", "Tumor/DCIS interface"]
    df = marker_summary.set_index("paper_marker_group").loc[order].reset_index()
    x = np.arange(len(df))
    ax.bar(x - 0.18, df["median_spearman_vs_full"], width=0.34, color="#2166ac", label="median rho")
    ax.bar(x + 0.18, df["best_cluster_match_rate"], width=0.34, color="#4c956c", label="match rate")
    ax.set_ylim(0, 1.08)
    ax.set_xticks(x)
    ax.set_xticklabels(["11q13", "Macro.", "Plasma", "Vascular", "Tumor/\nDCIS"], fontsize=6)
    ax.set_ylabel("fraction or rho")
    ax.set_title("Curated marker validation", fontsize=8, color=INK, pad=5)
    ax.legend(frameon=False, fontsize=5.8, loc="lower left", handlelength=1.0)
    style_axis(ax)


def plot_example_table(ax, metadata: pd.DataFrame):
    ax.set_axis_off()
    rows = []
    for gene in EXAMPLE_GENES:
        row = metadata.set_index("gene").loc[gene]
        full = str(row["full_best_cluster"])
        urna = str(row["urna_best_cluster"])
        rows.append(
            [
                gene,
                f"{float(row['urna_best_sss']):.3f}",
                f"{float(row['sss_spearman_vs_full']):.3f}",
                "yes" if bool(row["best_cluster_match"]) else "no",
                short_cluster(urna),
                short_cluster(full),
            ]
        )
    headers = ["gene", "uRNA\nSSS", "rho", "match", "uRNA best", "full best"]
    table = ax.table(cellText=rows, colLabels=headers, loc="center", cellLoc="left", colLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(5.8)
    table.scale(1, 1.35)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#d0d7df")
        cell.set_linewidth(0.35)
        if r == 0:
            cell.set_text_props(weight="bold", color=INK)
            cell.set_facecolor("#f1f4f7")
        elif c == 0:
            gene = rows[r - 1][0]
            cell.set_text_props(weight="bold", color=TARGET_COLORS[gene])


def short_cluster(name: str) -> str:
    replacements = {
        "11q13 Invasive Tumor Cells": "11q13 tumor",
        "Basal-like Structured DCIS Cells": "Basal-like DCIS",
        "CXCL14+ Fibroblasts": "CXCL14+ fib.",
        "Macrophages": "macrophage",
        "Plasma Cells": "plasma",
        "Endothelial Cells": "endothelial",
    }
    return replacements.get(name, name[:18])


def make_figure():
    setup_matplotlib()
    FIGURES.mkdir(parents=True, exist_ok=True)
    points, metadata, summary, sss, marker_summary = load_inputs()
    metadata = metadata.set_index("gene").loc[EXAMPLE_GENES].reset_index()

    fig = plt.figure(figsize=(7.15, 8.0), dpi=600)
    gs = GridSpec(
        4,
        5,
        figure=fig,
        height_ratios=[1.18, 1.22, 1.25, 1.00],
        hspace=0.55,
        wspace=0.32,
        top=0.93,
        bottom=0.04,
        left=0.06,
        right=0.985,
    )

    for idx, gene in enumerate(EXAMPLE_GENES):
        ax = fig.add_subplot(gs[0, idx])
        plot_spatial(ax, points, metadata.set_index("gene").loc[gene], gene)
        panel_label(ax, chr(ord("a") + idx), x=-0.08, y=1.12)

    ax_heat = fig.add_subplot(gs[1:3, 0:3])
    plot_sss_heatmap(ax_heat, sss, metadata)
    panel_label(ax_heat, "f", x=-0.08, y=1.04)

    ax_hist = fig.add_subplot(gs[1, 3:5])
    plot_concordance_context(ax_hist, summary, metadata)
    panel_label(ax_hist, "g", x=-0.10, y=1.12)

    ax_marker = fig.add_subplot(gs[2, 3:5])
    plot_marker_summary(ax_marker, marker_summary)
    panel_label(ax_marker, "h", x=-0.10, y=1.12)

    ax_table = fig.add_subplot(gs[3, :])
    plot_example_table(ax_table, metadata)
    panel_label(ax_table, "i", x=-0.012, y=1.08)

    legend_items = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#111827", markeredgewidth=0, markersize=4, alpha=0.55, label="uRNA"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#c7cdd4", markeredgewidth=0, markersize=4, alpha=0.70, label="other cells"),
        Line2D([0], [0], marker="o", color=INK, markerfacecolor="none", markersize=4, lw=0, label="target cells"),
    ]
    fig.legend(handles=legend_items, loc="upper center", bbox_to_anchor=(0.50, 0.992), frameon=False, ncol=3, fontsize=6.5)

    fig.savefig(FIGURES / "figure_5_spatial_evidence_multipanel.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "figure_5_spatial_evidence_multipanel.png", bbox_inches="tight", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    make_figure()
