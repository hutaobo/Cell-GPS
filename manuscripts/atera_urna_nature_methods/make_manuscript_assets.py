from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "benchmarking" / "atera_wta_breast_urna_segmentation_bias_results"
OUT = Path(__file__).resolve().parent
FIGURES = OUT / "figures"
TABLES = OUT / "tables"


MARKER_GENES = {
    "11q13": ["CCND1", "ELOVL2", "KCNJ3", "FGF19"],
    "Basal-like DCIS": ["KRT23", "DSC3", "SOSTDC1", "KLK5", "KLK7", "ITGB6", "MMP7"],
    "Macrophage": ["C3", "C1QA", "CSF1R", "CD163", "SIGLEC1"],
    "Plasma": ["IGHA1", "IGHM", "IGHA2", "JCHAIN"],
    "Rare epithelial": ["PIP", "HSPB8", "CLIC6", "TAT"],
    "Tumor/DCIS interface": ["FOXA1", "CTTN", "ANO1", "ESR1", "FADD", "PGR", "ERBB2", "FGF3"],
    "Vascular": ["EPAS1", "CDH5", "MMRN2"],
}

REPRESENTATIVE_CLUSTERS = [
    "11q13 Invasive Tumor Cells",
    "Basal-like Structured DCIS Cells",
    "CXCL14+ Fibroblasts",
    "Macrophages",
    "Plasma Cells",
    "Luminal-like Amorphous DCIS Cells",
    "Apocrine Cells",
    "Pericytes",
    "Endothelial Cells",
]

COMPACT_CLUSTERS = [
    "11q13 Invasive Tumor Cells",
    "Basal-like Structured DCIS Cells",
    "CXCL14+ Fibroblasts",
    "Macrophages",
    "Plasma Cells",
    "Endothelial Cells",
]

CLUSTER_LABELS = {
    "11q13 Invasive Tumor Cells": "11q13",
    "Basal-like Structured DCIS Cells": "Basal",
    "CXCL14+ Fibroblasts": "CXCL14+",
    "Macrophages": "Macro.",
    "Plasma Cells": "Plasma",
    "Luminal-like Amorphous DCIS Cells": "Luminal",
    "Apocrine Cells": "Apocrine",
    "Pericytes": "Pericyte",
    "Endothelial Cells": "Endoth.",
    "CAFs, DCIS Associated": "CAFs/DCIS",
}

SPATIAL_EXAMPLE_GENES = ["CCND1", "C1QA", "JCHAIN", "CDH5", "ERBB2"]
CONCORDANT_EXAMPLE_GENES = ["CCND1", "C1QA", "JCHAIN", "CDH5"]

GENE_COLORS = {
    "CCND1": "#b23a48",
    "C1QA": "#4c956c",
    "JCHAIN": "#6f5aa7",
    "CDH5": "#178f8f",
    "ERBB2": "#c48a22",
}

EXAMPLE_TARGETS = {
    "CCND1": "11q13 Invasive Tumor Cells",
    "C1QA": "Macrophages",
    "JCHAIN": "Plasma Cells",
    "CDH5": "Endothelial Cells",
    "ERBB2": "CXCL14+ Fibroblasts",
}

MARKER_GROUP_ORDER = [
    "11q13",
    "Basal-like DCIS",
    "Macrophage",
    "Plasma",
    "Rare epithelial",
    "Tumor/DCIS interface",
    "Vascular",
]

MARKER_GROUP_LABELS = {
    "11q13": "11q13 tumor",
    "Basal-like DCIS": "Basal-like\nDCIS",
    "Macrophage": "Macrophage",
    "Plasma": "Plasma",
    "Rare epithelial": "Rare\nepithelial",
    "Tumor/DCIS interface": "Tumor/DCIS\ninterface",
    "Vascular": "Vascular",
}

PALETTE = {
    "ink": "#1f2933",
    "muted": "#657786",
    "red": "#b23a48",
    "blue": "#2166ac",
    "teal": "#178f8f",
    "green": "#4c956c",
    "gold": "#c48a22",
    "gray": "#d7dde5",
    "light": "#f4f6f8",
}


def setup():
    FIGURES.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
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


def load_inputs():
    with open(SOURCE / "run_metadata.json", "r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    counts = pd.read_csv(SOURCE / "transcript_gene_counts.csv")
    summary = pd.read_csv(SOURCE / "urna_segmentation_bias_gene_summary.csv")
    urna_sss = pd.read_csv(SOURCE / "urna_only_t_and_c_result.csv", index_col=0)
    return metadata, counts, summary, urna_sss


def marker_group_for_gene(gene):
    for group, genes in MARKER_GENES.items():
        if gene in genes:
            return group
    return ""


def add_marker_groups(summary):
    summary = summary.copy()
    summary["curated_marker_group"] = summary["gene"].map(marker_group_for_gene)
    summary.loc[summary["marker_group"].isna(), "marker_group"] = ""
    summary["paper_marker_group"] = summary["curated_marker_group"]
    return summary


def save_tables(metadata, counts, summary):
    selected = summary.copy()
    marker_rows = selected[selected["paper_marker_group"] != ""].copy()

    overall = {
        "total_parquet_rows": metadata["total_rows"],
        "gene_rows": metadata["gene_rows"],
        "qv_ge_20_gene_rows": metadata["qv_gene_rows"],
        "qv_ge_20_unassigned_gene_rows": metadata["unassigned_rows"],
        "unassigned_fraction_of_qv_gene_rows": metadata["unassigned_rows"] / metadata["qv_gene_rows"],
        "genes_in_panel": metadata["n_genes"],
        "genes_with_qv_ge_20_uRNA": metadata["n_genes_with_unassigned"],
        "genes_with_qv_ge_20_uRNA_fraction": metadata["n_genes_with_unassigned"] / metadata["n_genes"],
        "genes_analyzed_by_urna_only_COSTE": metadata["n_urna_sss_genes"],
        "median_spearman_vs_full": selected["sss_spearman_vs_full"].median(),
        "mean_spearman_vs_full": selected["sss_spearman_vs_full"].mean(),
        "q25_spearman_vs_full": selected["sss_spearman_vs_full"].quantile(0.25),
        "q75_spearman_vs_full": selected["sss_spearman_vs_full"].quantile(0.75),
        "best_cluster_exact_match_rate": selected["best_cluster_match"].mean(),
        "genes_spearman_lt_0_5": int((selected["sss_spearman_vs_full"] < 0.5).sum()),
        "median_gene_uRNA_fraction_all_panel_genes": counts["unassigned_fraction"].median(),
        "median_gene_uRNA_fraction_analyzed_genes": selected["unassigned_fraction"].median(),
    }
    overall_table = pd.DataFrame(
        [{"metric": key, "value": value} for key, value in overall.items()]
    )
    overall_table.to_csv(TABLES / "summary_statistics.csv", index=False)

    cluster_counts = (
        selected.groupby("urna_best_cluster", dropna=False)
        .size()
        .reset_index(name="n_genes")
        .sort_values("n_genes", ascending=False)
    )
    cluster_counts.to_csv(TABLES / "celltype_attribution_counts.csv", index=False)

    by_marker = (
        marker_rows.groupby("paper_marker_group")
        .agg(
            n_genes=("gene", "size"),
            median_unassigned_count=("unassigned_count", "median"),
            median_unassigned_fraction=("unassigned_fraction", "median"),
            median_spearman_vs_full=("sss_spearman_vs_full", "median"),
            best_cluster_match_rate=("best_cluster_match", "mean"),
        )
        .reset_index()
        .sort_values("paper_marker_group")
    )
    by_marker.to_csv(TABLES / "marker_group_summary.csv", index=False)

    marker_rows = marker_rows[
        [
            "gene",
            "paper_marker_group",
            "unassigned_count",
            "unassigned_fraction",
            "urna_best_cluster",
            "urna_best_sss",
            "full_best_cluster",
            "full_best_sss",
            "best_cluster_match",
            "sss_spearman_vs_full",
        ]
    ].sort_values(["paper_marker_group", "gene"])
    marker_rows.to_csv(TABLES / "marker_validation_table.csv", index=False)

    discordant = selected.sort_values("sss_spearman_vs_full", ascending=True).head(30)
    discordant[
        [
            "gene",
            "unassigned_count",
            "unassigned_fraction",
            "urna_best_cluster",
            "full_best_cluster",
            "best_cluster_match",
            "sss_spearman_vs_full",
        ]
    ].to_csv(TABLES / "lowest_concordance_genes.csv", index=False)

    top_urna = counts.sort_values("unassigned_count", ascending=False).head(50)
    top_urna.to_csv(TABLES / "top_50_uRNA_genes.csv", index=False)

    summary_lines = [
        "# Paper-level Atera uRNA/COSTE summary",
        "",
        f"High-quality gene transcripts (QV >= 20): {metadata['qv_gene_rows']:,}.",
        f"High-quality unassigned gene transcripts: {metadata['unassigned_rows']:,} ({overall['unassigned_fraction_of_qv_gene_rows']:.1%}).",
        f"Genes with at least one high-quality uRNA: {metadata['n_genes_with_unassigned']:,} of {metadata['n_genes']:,} ({overall['genes_with_qv_ge_20_uRNA_fraction']:.2%}).",
        f"Genes analyzed by uRNA-only COSTE: {metadata['n_urna_sss_genes']:,}.",
        f"Median Spearman concordance against all-transcript SSS: {overall['median_spearman_vs_full']:.3f} (IQR {overall['q25_spearman_vs_full']:.3f}-{overall['q75_spearman_vs_full']:.3f}).",
        f"Best-cluster exact match rate: {overall['best_cluster_exact_match_rate']:.3f}.",
        f"Genes with Spearman concordance < 0.5: {overall['genes_spearman_lt_0_5']:,}.",
        "",
        "Top uRNA-only COSTE attributions:",
    ]
    for _, row in cluster_counts.head(8).iterrows():
        summary_lines.append(f"- {row['urna_best_cluster']}: {int(row['n_genes'])} genes")
    (OUT / "paper_level_summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    return overall, cluster_counts, by_marker, marker_rows


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#b7c0cc")
    ax.spines["bottom"].set_color("#b7c0cc")
    ax.tick_params(colors=PALETTE["ink"], labelsize=8)
    ax.title.set_color(PALETTE["ink"])
    ax.yaxis.label.set_color(PALETTE["ink"])
    ax.xaxis.label.set_color(PALETTE["ink"])


def panel_label(ax, label, x=-0.13, y=1.08):
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        color=PALETTE["ink"],
        ha="left",
        va="top",
    )


def fmt_millions(value):
    return f"{value / 1e6:.1f}M"


def marker_rows(summary):
    return summary.loc[summary["paper_marker_group"] != ""].copy()


def marker_group_summary(summary):
    rows = marker_rows(summary)
    grouped = (
        rows.groupby("paper_marker_group")
        .agg(
            n_genes=("gene", "size"),
            median_unassigned_count=("unassigned_count", "median"),
            median_unassigned_fraction=("unassigned_fraction", "median"),
            median_spearman_vs_full=("sss_spearman_vs_full", "median"),
            best_cluster_match_rate=("best_cluster_match", "mean"),
        )
        .reindex(MARKER_GROUP_ORDER)
        .reset_index()
    )
    return grouped.dropna(subset=["n_genes"])


def median_marker_heat(summary, urna_sss, groups, clusters):
    rows = marker_rows(summary).set_index("gene")
    heat = []
    for group in groups:
        genes = [gene for gene, row in rows.iterrows() if row["paper_marker_group"] == group and gene in urna_sss.index]
        heat.append((1 - urna_sss.loc[genes, clusters].astype(float).clip(0, 1)).median(axis=0).to_numpy())
    return np.vstack(heat)


def figure_graphical_abstract(metadata, overall):
    fig, ax = plt.subplots(figsize=(12, 4), dpi=220)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    boxes = [
        (0.04, 0.22, 0.20, 0.56, "Atera WTA breast section", "18,028 genes\n624.1M high-quality gene transcripts"),
        (0.30, 0.22, 0.20, 0.56, "uRNA extraction", "122.7M high-quality unassigned RNAs\n99.97% genes represented"),
        (0.56, 0.22, 0.20, 0.56, "uRNA-only COSTE", "single transcripts treated as nodes\nSpatial Separation Score profiles"),
        (0.82, 0.22, 0.14, 0.56, "Orthogonal check", f"median rho={overall['median_spearman_vs_full']:.3f}\nmatch rate={overall['best_cluster_exact_match_rate']:.3f}"),
    ]
    for x, y, w, h, title, text in boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor="#ffffff", edgecolor="#8aa0b8", linewidth=1.4)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h - 0.12, title, ha="center", va="center", fontsize=10, fontweight="bold", color=PALETTE["ink"])
        ax.text(x + w / 2, y + 0.22, text, ha="center", va="center", fontsize=8.3, color=PALETTE["muted"], linespacing=1.4)
    for x in [0.255, 0.515, 0.775]:
        ax.annotate("", xy=(x + 0.035, 0.5), xytext=(x, 0.5), arrowprops=dict(arrowstyle="-|>", lw=1.6, color=PALETTE["teal"]))

    rng = np.random.default_rng(7)
    for _ in range(85):
        px = rng.uniform(0.065, 0.215)
        py = rng.uniform(0.30, 0.60)
        color = PALETTE["red"] if rng.random() < 0.20 else PALETTE["blue"]
        ax.scatter(px, py, s=rng.uniform(3, 10), color=color, alpha=0.65)
    ax.text(0.04, 0.08, "Concept: unassigned RNAs are retained as spatial evidence instead of discarded as segmentation waste.", fontsize=9, color=PALETTE["ink"])
    fig.tight_layout()
    fig.savefig(FIGURES / "graphical_abstract.png", bbox_inches="tight")
    fig.savefig(FIGURES / "graphical_abstract.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_overview(metadata, counts, summary, cluster_counts):
    fig = plt.figure(figsize=(7.15, 5.6), dpi=600)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[0.95, 1.05], hspace=0.46, wspace=0.34)

    ax = fig.add_subplot(gs[0, 0])
    labels = ["all rows", "gene rows", "QV>=20\ngene", "QV>=20\nuRNA"]
    values = [metadata["total_rows"], metadata["gene_rows"], metadata["qv_gene_rows"], metadata["unassigned_rows"]]
    colors = [PALETTE["gray"], "#9eb3c7", PALETTE["blue"], PALETTE["red"]]
    bars = ax.bar(labels, np.array(values) / 1e6, color=colors, width=0.72)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 18, fmt_millions(value), ha="center", va="bottom", fontsize=6.2, color=PALETTE["ink"])
    ax.set_ylabel("transcripts (millions)")
    ax.set_title("Atera WTA transcript scale", fontsize=8, pad=4)
    ax.tick_params(axis="x", rotation=0, labelsize=6.4)
    style_axes(ax)
    panel_label(ax, "a", x=-0.16)

    ax = fig.add_subplot(gs[0, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    stages = [
        ("18,028", "assayed panel\ngenes"),
        ("18,023", "genes with\n>=1 uRNA"),
        ("top 500 high-uRNA\n+ curated markers", "selection rule"),
        ("530", "analyzed genes\nfor COSTE"),
    ]
    y_positions = [0.74, 0.55, 0.36, 0.17]
    widths = [0.78, 0.74, 0.70, 0.48]
    for idx, ((number, label), y0, width) in enumerate(zip(stages, y_positions, widths)):
        x0 = 0.50 - width / 2
        rect = plt.Rectangle((x0, y0), width, 0.125, facecolor="#ffffff", edgecolor="#9eb3c7", linewidth=0.9)
        ax.add_patch(rect)
        ax.text(x0 + 0.045, y0 + 0.064, number, ha="left", va="center", fontsize=7.1, fontweight="bold", color=PALETTE["ink"], linespacing=0.9)
        ax.text(x0 + width - 0.035, y0 + 0.064, label, ha="right", va="center", fontsize=5.4, color=PALETTE["muted"], linespacing=0.95)
        if idx < len(stages) - 1:
            ax.annotate("", xy=(0.50, y_positions[idx + 1] + 0.14), xytext=(0.50, y0 - 0.012), arrowprops=dict(arrowstyle="-|>", lw=0.9, color=PALETTE["teal"]))
    ax.text(0.5, 0.03, "All-panel coverage is shown for 18,028 genes; COSTE/concordance is limited to the 530 analyzed genes.", ha="center", va="center", fontsize=5.6, color=PALETTE["ink"])
    ax.set_title("Gene universe and analysis set", fontsize=8, pad=4)
    panel_label(ax, "b", x=-0.12)

    ax = fig.add_subplot(gs[1, 0])
    bins = np.linspace(0, counts["unassigned_fraction"].quantile(0.995), 48)
    ax.hist(counts["unassigned_fraction"], bins=bins, color="#cfd6df", edgecolor="white", linewidth=0.25, label="all panel genes")
    ax.hist(summary["unassigned_fraction"], bins=bins, histtype="step", color=PALETTE["teal"], linewidth=1.3, label="530 analyzed genes")
    median_all = counts["unassigned_fraction"].median()
    ax.axvline(median_all, color=PALETTE["ink"], linewidth=1.0)
    ax.text(median_all + 0.006, ax.get_ylim()[1] * 0.92, "median\n0.160", fontsize=5.9, color=PALETTE["ink"], va="top")
    ax.set_xlabel("uRNA fraction per gene")
    ax.set_ylabel("genes")
    ax.set_title("uRNA fraction across the full panel", fontsize=8, pad=4)
    ax.legend(frameon=False, fontsize=5.8, loc="upper right")
    style_axes(ax)
    panel_label(ax, "c", x=-0.16)

    ax = fig.add_subplot(gs[1, 1])
    top = counts.sort_values("unassigned_count", ascending=False).head(12).sort_values("unassigned_count")
    colors = [GENE_COLORS.get(gene, PALETTE["gold"]) for gene in top["gene"]]
    ax.barh(top["gene"], top["unassigned_count"] / 1e3, color=colors)
    ax.set_xlabel("uRNA count (thousands)")
    ax.set_title("High-coverage uRNA genes include markers", fontsize=8, pad=4)
    style_axes(ax)
    panel_label(ax, "d", x=-0.16)

    fig.tight_layout()
    fig.savefig(FIGURES / "legacy_dataset_and_attribution.png", bbox_inches="tight")
    fig.savefig(FIGURES / "legacy_dataset_and_attribution.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_concordance(summary, overall):
    fig = plt.figure(figsize=(7.15, 5.6), dpi=600)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 0.95], hspace=0.45, wspace=0.34)
    df = summary.copy()

    ax = fig.add_subplot(gs[0, 0])
    counts = (
        df.groupby("urna_best_cluster", dropna=False)
        .size()
        .reset_index(name="n_genes")
        .sort_values("n_genes", ascending=False)
        .head(10)
        .sort_values("n_genes")
    )
    ax.barh(counts["urna_best_cluster"], counts["n_genes"], color=PALETTE["green"])
    ax.set_xlabel("genes")
    ax.set_title("uRNA-only best COSTE cluster", fontsize=8, pad=4)
    style_axes(ax)
    panel_label(ax, "a", x=-0.16)

    ax = fig.add_subplot(gs[0, 1])
    match_colors = np.where(df["best_cluster_match"], PALETTE["blue"], PALETTE["red"])
    ax.scatter(
        df["unassigned_count"],
        df["sss_spearman_vs_full"],
        s=20,
        c=match_colors,
        alpha=0.70,
        linewidths=0,
    )
    ax.set_xscale("log")
    ax.axhline(overall["median_spearman_vs_full"], color=PALETTE["ink"], linewidth=1)
    ax.set_xlabel("uRNA count per gene")
    ax.set_ylabel("Spearman rho vs all-transcript SSS")
    ax.set_title("Coverage and concordance", fontsize=8, pad=4)
    style_axes(ax)
    panel_label(ax, "b", x=-0.16)

    ax = fig.add_subplot(gs[1, 0])
    ax.hist(df["sss_spearman_vs_full"], bins=np.linspace(-0.1, 1.05, 35), color=PALETTE["blue"], alpha=0.85)
    ax.axvline(overall["median_spearman_vs_full"], color=PALETTE["ink"], linewidth=1.3, label="median")
    ax.text(overall["median_spearman_vs_full"] + 0.02, ax.get_ylim()[1] * 0.88, f"median rho\n{overall['median_spearman_vs_full']:.3f}", fontsize=6, color=PALETTE["ink"], va="top")
    ax.set_xlabel("Spearman rho")
    ax.set_ylabel("genes")
    ax.set_title("Per-gene SSS concordance", fontsize=8, pad=4)
    style_axes(ax)
    panel_label(ax, "c", x=-0.16)

    ax = fig.add_subplot(gs[1, 1])
    metrics = pd.Series(
        {
            "best-cluster\nmatch": overall["best_cluster_exact_match_rate"],
            "rho >= 0.5": (df["sss_spearman_vs_full"] >= 0.5).mean(),
            "rho >= 0.8": (df["sss_spearman_vs_full"] >= 0.8).mean(),
            "rho = 1": (df["sss_spearman_vs_full"] >= 0.999).mean(),
        }
    )
    ax.bar(metrics.index, metrics.values, color=[PALETTE["blue"], PALETTE["green"], PALETTE["gold"], PALETTE["teal"]])
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of genes")
    ax.set_title("Concordance thresholds", fontsize=8, pad=4)
    style_axes(ax)
    panel_label(ax, "d", x=-0.16)

    fig.tight_layout()
    fig.savefig(FIGURES / "legacy_concordance.png", bbox_inches="tight")
    fig.savefig(FIGURES / "legacy_concordance.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_marker_validation(summary, urna_sss):
    groups = [group for group in MARKER_GROUP_ORDER if group in set(marker_rows(summary)["paper_marker_group"])]
    heat = median_marker_heat(summary, urna_sss, groups, REPRESENTATIVE_CLUSTERS)
    grouped = marker_group_summary(summary).set_index("paper_marker_group").loc[groups].reset_index()

    fig = plt.figure(figsize=(7.15, 4.9), dpi=600)
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.35, 1.0], height_ratios=[1.0, 0.85], hspace=0.42, wspace=0.42)

    ax = fig.add_subplot(gs[:, 0])
    image = ax.imshow(heat, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(REPRESENTATIVE_CLUSTERS)))
    ax.set_xticklabels([CLUSTER_LABELS[c] for c in REPRESENTATIVE_CLUSTERS], rotation=45, ha="right", fontsize=6.0)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([MARKER_GROUP_LABELS[g] for g in groups], fontsize=6.2)
    ax.tick_params(length=0)
    ax.set_title("Median uRNA-only COSTE proximity by marker group", fontsize=8, pad=5)
    expected = {
        "11q13": ["11q13 Invasive Tumor Cells"],
        "Basal-like DCIS": ["Basal-like Structured DCIS Cells"],
        "Macrophage": ["Macrophages"],
        "Plasma": ["Plasma Cells"],
        "Rare epithelial": ["Luminal-like Amorphous DCIS Cells", "Apocrine Cells"],
        "Tumor/DCIS interface": ["Basal-like Structured DCIS Cells", "11q13 Invasive Tumor Cells", "CXCL14+ Fibroblasts"],
        "Vascular": ["Pericytes", "Endothelial Cells"],
    }
    for i, group in enumerate(groups):
        for cluster in expected.get(group, []):
            if cluster in REPRESENTATIVE_CLUSTERS:
                j = REPRESENTATIVE_CLUSTERS.index(cluster)
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#000000", lw=1.2))
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(image, ax=ax, fraction=0.028, pad=0.025)
    cbar.set_label("1 - SSS", fontsize=6, labelpad=1.0)
    cbar.ax.tick_params(labelsize=5.6, length=2)
    panel_label(ax, "a", x=-0.20, y=1.05)

    ax = fig.add_subplot(gs[0, 1])
    y = np.arange(len(grouped))
    colors = [PALETTE["red"] if group == "Tumor/DCIS interface" else PALETTE["blue"] for group in grouped["paper_marker_group"]]
    ax.barh(y - 0.15, grouped["median_spearman_vs_full"], height=0.28, color=colors, label="median rho")
    ax.barh(y + 0.15, grouped["best_cluster_match_rate"], height=0.28, color=PALETTE["green"], label="match rate")
    ax.set_yticks(y)
    ax.set_yticklabels([MARKER_GROUP_LABELS[g].replace("\n", " ") for g in grouped["paper_marker_group"]], fontsize=5.8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("fraction or rho")
    ax.set_title("Marker group recovery", fontsize=8, pad=5)
    ax.legend(frameon=False, fontsize=5.8, loc="upper left", bbox_to_anchor=(0.0, -0.18), ncol=2, handlelength=1.1)
    style_axes(ax)
    panel_label(ax, "b", x=-0.20, y=1.10)

    ax = fig.add_subplot(gs[1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    callouts = [
        ("11q13", "CCND1", "tumor"),
        ("Macrophage", "C1QA", "immune"),
        ("Plasma", "JCHAIN", "immune"),
        ("Vascular", "CDH5", "vascular"),
    ]
    for idx, (group, gene, compartment) in enumerate(callouts):
        x = 0.04 + idx * 0.235
        rect = plt.Rectangle((x, 0.24), 0.20, 0.50, facecolor="#ffffff", edgecolor=GENE_COLORS[gene], linewidth=1.0)
        ax.add_patch(rect)
        ax.text(x + 0.10, 0.58, gene, ha="center", va="center", fontsize=7.2, fontweight="bold", color=GENE_COLORS[gene])
        ax.text(x + 0.10, 0.40, compartment, ha="center", va="center", fontsize=5.8, color=PALETTE["muted"])
    ax.text(0.5, 0.07, "These concordant marker groups provide interpretable positive controls for spatial inspection.", ha="center", va="center", fontsize=6.0, color=PALETTE["ink"])
    ax.set_title("Positive-control groups selected for spatial examples", fontsize=7.4, pad=5)
    panel_label(ax, "c", x=-0.10, y=1.30)

    fig.tight_layout()
    fig.savefig(FIGURES / "legacy_marker_validation.png", bbox_inches="tight")
    fig.savefig(FIGURES / "legacy_marker_validation.pdf", bbox_inches="tight")
    plt.close(fig)

    marker_gene_order = [gene for genes in MARKER_GENES.values() for gene in genes if gene in urna_sss.index]
    marker_summary = summary[summary["gene"].isin(marker_gene_order)][
        ["gene", "paper_marker_group", "urna_best_cluster", "full_best_cluster", "best_cluster_match", "sss_spearman_vs_full"]
    ].copy()
    marker_summary.to_csv(TABLES / "legacy_marker_callouts.csv", index=False)


def figure_gene_selection_rationale(metadata, counts, summary, urna_sss, overall):
    rows = summary.set_index("gene").loc[SPATIAL_EXAMPLE_GENES].reset_index()
    heat = 1 - urna_sss.loc[SPATIAL_EXAMPLE_GENES, COMPACT_CLUSTERS].astype(float).clip(0, 1)

    fig = plt.figure(figsize=(7.15, 5.6), dpi=600)
    gs = GridSpec(3, 3, figure=fig, height_ratios=[0.72, 1.15, 1.05], width_ratios=[1.0, 1.0, 1.05], hspace=0.58, wspace=0.52)

    ax = fig.add_subplot(gs[0, :])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    stages = [
        ("18,028", "assayed panel\ngenes"),
        ("18,023", "with >=1\nuRNA"),
        ("top 500\n+ markers", "selection\nrule"),
        ("530", "COSTE\nanalyzed genes"),
        ("marker controls", "expected\ncompartments"),
        ("4 + 1", "spatial examples"),
    ]
    x_positions = [0.010, 0.172, 0.334, 0.500, 0.664, 0.835]
    widths = [0.118, 0.118, 0.124, 0.120, 0.128, 0.118]
    for idx, ((number, label), x, width) in enumerate(zip(stages, x_positions, widths)):
        face = "#f5f7fa" if idx < 4 else "#ffffff"
        edge = PALETTE["teal"] if idx == 4 else "#9eb3c7"
        rect = plt.Rectangle((x, 0.24), width, 0.50, facecolor=face, edgecolor=edge, linewidth=1.0)
        ax.add_patch(rect)
        num_size = 7.4
        if idx in {3, 4}:
            num_size = 6.1
        ax.text(x + width / 2, 0.58, number, ha="center", va="center", fontsize=num_size, fontweight="bold", color=PALETTE["ink"])
        ax.text(x + width / 2, 0.39, label, ha="center", va="center", fontsize=5.4, color=PALETTE["muted"], linespacing=1.0)
        if idx < len(stages) - 1:
            ax.annotate("", xy=(x_positions[idx + 1] - 0.012, 0.50), xytext=(x + width + 0.010, 0.50), arrowprops=dict(arrowstyle="-|>", lw=0.9, color=PALETTE["teal"]))
    ax.text(0.5, 0.08, "Four examples are concordant positive controls; ERBB2 is retained as a boundary-sensitive discordant control.", ha="center", va="center", fontsize=6.2, color=PALETTE["ink"])
    ax.set_title("Selection path for the spatial examples", fontsize=8, pad=5)
    panel_label(ax, "a", x=-0.035, y=1.10)

    ax = fig.add_subplot(gs[1, 0])
    group_df = marker_group_summary(summary)
    focus_groups = ["11q13", "Macrophage", "Plasma", "Vascular", "Tumor/DCIS interface"]
    focus = group_df.set_index("paper_marker_group").loc[focus_groups].reset_index()
    y = np.arange(len(focus))
    ax.barh(y, focus["median_unassigned_count"] / 1e3, color=[GENE_COLORS.get("ERBB2", PALETTE["gold"]) if g == "Tumor/DCIS interface" else PALETTE["blue"] for g in focus["paper_marker_group"]])
    ax.set_yticks(y)
    ax.set_yticklabels([MARKER_GROUP_LABELS[g].replace("\n", " ") for g in focus["paper_marker_group"]], fontsize=5.7)
    ax.invert_yaxis()
    ax.set_xlabel("median uRNA count (thousands)")
    ax.set_title("Marker groups have uRNA support", fontsize=8, pad=5)
    style_axes(ax)
    panel_label(ax, "b", x=-0.32, y=1.10)

    ax = fig.add_subplot(gs[1, 1:])
    image = ax.imshow(heat.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(COMPACT_CLUSTERS)))
    ax.set_xticklabels([CLUSTER_LABELS[c] for c in COMPACT_CLUSTERS], fontsize=6)
    ax.set_yticks(range(len(SPATIAL_EXAMPLE_GENES)))
    ax.set_yticklabels(SPATIAL_EXAMPLE_GENES, fontsize=6.2)
    ax.tick_params(length=0)
    for i, gene in enumerate(SPATIAL_EXAMPLE_GENES):
        cluster = EXAMPLE_TARGETS[gene]
        if cluster in COMPACT_CLUSTERS:
            j = COMPACT_CLUSTERS.index(cluster)
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#000000", lw=1.2))
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(image, ax=ax, fraction=0.026, pad=0.015)
    cbar.set_label("1 - SSS", fontsize=6, labelpad=1.0)
    cbar.ax.tick_params(labelsize=5.6, length=2)
    ax.set_title("uRNA-only proximity identifies expected targets", fontsize=8, pad=5)
    panel_label(ax, "c", x=-0.12, y=1.10)

    ax = fig.add_subplot(gs[2, :])
    ax.set_axis_off()
    table_rows = []
    for _, row in rows.iterrows():
        gene = row["gene"]
        table_rows.append(
            [
                gene,
                f"{row['unassigned_count'] / 1e3:.1f}k",
                f"{row['unassigned_fraction']:.3f}",
                CLUSTER_LABELS.get(EXAMPLE_TARGETS[gene], EXAMPLE_TARGETS[gene]),
                f"{1 - float(row['urna_best_sss']):.3f}",
                f"{float(row['sss_spearman_vs_full']):.3f}",
                "match" if bool(row["best_cluster_match"]) else "discordant",
            ]
        )
    headers = ["gene", "uRNA\ncount", "uRNA\nfraction", "target", "target\nproximity", "rho", "call"]
    table = ax.table(cellText=table_rows, colLabels=headers, loc="center", cellLoc="left", colLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(5.8)
    table.scale(1.0, 1.25)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#d0d7df")
        cell.set_linewidth(0.35)
        if r == 0:
            cell.set_text_props(weight="bold", color=PALETTE["ink"])
            cell.set_facecolor("#f1f4f7")
        elif c == 0:
            gene = table_rows[r - 1][0]
            cell.set_text_props(weight="bold", color=GENE_COLORS[gene])
        elif c == 6:
            gene = table_rows[r - 1][0]
            cell.set_text_props(color=GENE_COLORS[gene])
    ax.set_title("Per-gene evidence carried into the spatial validation figure", fontsize=8, pad=8)
    panel_label(ax, "d", x=-0.05, y=1.08)

    fig.tight_layout()
    fig.savefig(FIGURES / "legacy_gene_selection_rationale.png", bbox_inches="tight")
    fig.savefig(FIGURES / "legacy_gene_selection_rationale.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_combined_lead_in(metadata, counts, summary, overall, cluster_counts):
    df = summary.copy()
    marker_summary = marker_group_summary(summary).set_index("paper_marker_group")

    fig = plt.figure(figsize=(7.15, 7.12), dpi=600)
    gs = GridSpec(
        4,
        6,
        figure=fig,
        height_ratios=[1.0, 1.08, 0.94, 0.86],
        hspace=0.36,
        wspace=0.72,
    )

    ax = fig.add_subplot(gs[0, 0:2])
    labels = ["all\nrows", "gene\nrows", "QV>=20\ngene", "QV>=20\nuRNA"]
    values = [
        metadata["total_rows"],
        metadata["gene_rows"],
        metadata["qv_gene_rows"],
        metadata["unassigned_rows"],
    ]
    colors = ["#cfd6df", "#9eb3c7", PALETTE["blue"], PALETTE["red"]]
    bars = ax.bar(np.arange(len(values)), np.array(values) / 1e6, color=colors, width=0.68)
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 18,
            fmt_millions(value),
            ha="center",
            va="bottom",
            fontsize=5.4,
            color=PALETTE["ink"],
        )
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=5.8)
    ax.set_ylabel("transcripts (millions)")
    ax.set_title("Molecular scale", fontsize=7.6, pad=4)
    ax.text(
        0.98,
        0.78,
        "19.7%\nhigh-quality\nuRNAs",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.5,
        color=PALETTE["red"],
        linespacing=1.05,
    )
    style_axes(ax)
    panel_label(ax, "a", x=-0.20, y=1.12)

    ax = fig.add_subplot(gs[0, 2:4])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    stages = [
        ("18,028", "assayed panel genes"),
        ("18,023", "genes with >=1 uRNA"),
        ("top 500\n+ markers", "selection rule"),
        ("530", "COSTE/concordance"),
    ]
    y_positions = [0.78, 0.57, 0.36, 0.15]
    widths = [0.90, 0.84, 0.78, 0.70]
    for idx, ((number, label), y0, width) in enumerate(zip(stages, y_positions, widths)):
        x0 = 0.50 - width / 2
        rect = plt.Rectangle(
            (x0, y0),
            width,
            0.125,
            facecolor="#ffffff" if idx == 3 else "#f5f7fa",
            edgecolor=PALETTE["teal"] if idx == 3 else "#a9b8c8",
            linewidth=0.9,
        )
        ax.add_patch(rect)
        ax.text(
            x0 + 0.035,
            y0 + 0.063,
            number,
            ha="left",
            va="center",
            fontsize=5.8 if idx == 2 else 6.8,
            fontweight="bold",
            color=PALETTE["ink"],
            linespacing=0.85,
        )
        ax.text(
            x0 + width - 0.030,
            y0 + 0.063,
            label,
            ha="right",
            va="center",
            fontsize=4.8,
            color=PALETTE["muted"],
            linespacing=0.95,
        )
        if idx < len(stages) - 1:
            ax.annotate(
                "",
                xy=(0.50, y_positions[idx + 1] + 0.140),
                xytext=(0.50, y0 - 0.012),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=PALETTE["teal"]),
            )
    ax.text(
        0.5,
        0.015,
        "Only these 530 genes enter COSTE/concordance.",
        ha="center",
        va="bottom",
        fontsize=5.2,
        color=PALETTE["ink"],
    )
    ax.set_title("Analysis universe", fontsize=7.6, pad=4)
    panel_label(ax, "b", x=-0.14, y=1.12)

    ax = fig.add_subplot(gs[0, 4:6])
    max_fraction = counts["unassigned_fraction"].quantile(0.995)
    bins = np.linspace(0, max_fraction, 44)
    ax.hist(
        counts["unassigned_fraction"].clip(upper=max_fraction),
        bins=bins,
        color="#cfd6df",
        edgecolor="white",
        linewidth=0.25,
        label="18,028 panel",
    )
    ax.hist(
        df["unassigned_fraction"].clip(upper=max_fraction),
        bins=bins,
        histtype="step",
        color=PALETTE["teal"],
        linewidth=1.15,
        label="530 analyzed",
    )
    median_all = counts["unassigned_fraction"].median()
    ax.axvline(median_all, color=PALETTE["ink"], linewidth=0.9)
    ax.text(
        median_all + 0.010,
        ax.get_ylim()[1] * 0.88,
        "median\n0.160",
        fontsize=5.4,
        color=PALETTE["ink"],
        va="top",
    )
    ax.set_xlabel("uRNA fraction per gene")
    ax.set_ylabel("genes")
    ax.set_title("Full-panel uRNA coverage", fontsize=7.6, pad=4)
    ax.legend(frameon=False, fontsize=5.0, loc="upper right", handlelength=1.0)
    style_axes(ax)
    panel_label(ax, "c", x=-0.20, y=1.12)

    ax = fig.add_subplot(gs[1, 0:3])
    top_clusters = cluster_counts.head(8).copy()
    top_clusters["label"] = top_clusters["urna_best_cluster"].map(
        lambda value: CLUSTER_LABELS.get(value, value.replace(" Cells", ""))
    )
    top_clusters = top_clusters.sort_values("n_genes")
    bar_colors = [
        PALETTE["green"] if "Macrophage" in cluster else PALETTE["blue"]
        for cluster in top_clusters["urna_best_cluster"]
    ]
    ax.barh(top_clusters["label"], top_clusters["n_genes"], color=bar_colors)
    for y, value in enumerate(top_clusters["n_genes"]):
        ax.text(value + 4, y, str(int(value)), va="center", fontsize=5.6, color=PALETTE["ink"])
    ax.set_xlim(0, max(top_clusters["n_genes"]) * 1.22)
    ax.set_xlabel("genes among 530")
    ax.set_title("uRNA-only COSTE best cell-group calls", fontsize=7.6, pad=4)
    style_axes(ax)
    panel_label(ax, "d", x=-0.15, y=1.10)

    ax = fig.add_subplot(gs[1, 3:6])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    metrics = [
        ("median rho", overall["median_spearman_vs_full"], PALETTE["blue"]),
        ("best-group match", overall["best_cluster_exact_match_rate"], PALETTE["green"]),
        ("rho >= 0.5", (df["sss_spearman_vs_full"] >= 0.5).mean(), PALETTE["gold"]),
        ("rho >= 0.8", (df["sss_spearman_vs_full"] >= 0.8).mean(), PALETTE["teal"]),
    ]
    for idx, (label, value, color) in enumerate(metrics):
        row = 3 - idx
        y0 = 0.12 + row * 0.19
        ax.text(0.02, y0 + 0.045, label, ha="left", va="center", fontsize=5.8, color=PALETTE["ink"])
        ax.add_patch(plt.Rectangle((0.36, y0 + 0.018), 0.50, 0.052, facecolor="#edf1f5", edgecolor="none"))
        ax.add_patch(plt.Rectangle((0.36, y0 + 0.018), 0.50 * value, 0.052, facecolor=color, edgecolor="none"))
        ax.text(0.89, y0 + 0.045, f"{value:.3f}", ha="left", va="center", fontsize=5.8, color=PALETTE["ink"])
    ax.text(
        0.02,
        0.035,
        "Middle-ground agreement: uRNAs recover many compartments\nwithout copying the cell-assigned signal.",
        ha="left",
        va="bottom",
        fontsize=5.25,
        color=PALETTE["muted"],
        linespacing=1.1,
    )
    ax.set_title("Global concordance, summarized", fontsize=7.6, pad=4)
    panel_label(ax, "e", x=-0.10, y=1.10)

    ax = fig.add_subplot(gs[2, :])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    marker_flow = [
        ("11q13", "11q13 tumor", "CCND1"),
        ("Basal-like DCIS", "Basal DCIS", "KRT23"),
        ("Macrophage", "Macrophage", "C1QA"),
        ("Plasma", "Plasma", "JCHAIN"),
        ("Rare epithelial", "Luminal / apocrine", "PIP"),
        ("Vascular", "Endoth. / pericyte", "CDH5"),
        ("Tumor/DCIS interface", "mixed boundary", "ERBB2"),
    ]
    x_edges = np.linspace(0.015, 0.985, len(marker_flow) + 1)
    for idx, (group, target, exemplar) in enumerate(marker_flow):
        x0 = x_edges[idx] + 0.006
        width = x_edges[idx + 1] - x_edges[idx] - 0.012
        edge = PALETTE["gold"] if group == "Tumor/DCIS interface" else "#9eb3c7"
        face = "#fff8ed" if group == "Tumor/DCIS interface" else "#ffffff"
        rect = plt.Rectangle((x0, 0.15), width, 0.60, facecolor=face, edgecolor=edge, linewidth=0.85)
        ax.add_patch(rect)
        if group in marker_summary.index:
            match = marker_summary.loc[group, "best_cluster_match_rate"]
            rho = marker_summary.loc[group, "median_spearman_vs_full"]
        else:
            match = np.nan
            rho = np.nan
        status = "mixed" if group == "Tumor/DCIS interface" else "recovered"
        ax.text(
            x0 + width / 2,
            0.66,
            MARKER_GROUP_LABELS[group].replace("\n", " "),
            ha="center",
            va="center",
            fontsize=5.4,
            fontweight="bold",
            color=PALETTE["ink"],
        )
        ax.text(
            x0 + width / 2,
            0.50,
            target,
            ha="center",
            va="center",
            fontsize=5.0,
            color=PALETTE["muted"],
        )
        ax.text(
            x0 + width / 2,
            0.35,
            f"{status}\nmatch {match:.2g}; rho {rho:.2g}",
            ha="center",
            va="center",
            fontsize=4.75,
            color=PALETTE["gold"] if group == "Tumor/DCIS interface" else PALETTE["teal"],
            linespacing=1.0,
        )
        ax.text(
            x0 + width / 2,
            0.075,
            exemplar,
            ha="center",
            va="center",
            fontsize=5.1,
            color=GENE_COLORS.get(exemplar, PALETTE["ink"]),
            fontweight="bold" if exemplar in GENE_COLORS else "normal",
        )
    ax.text(
        0.015,
        0.93,
        "Curated markers define interpretable positive controls before local spatial inspection",
        ha="left",
        va="center",
        fontsize=6.0,
        color=PALETTE["ink"],
    )
    ax.set_title("Marker-control logic", fontsize=7.6, pad=4)
    panel_label(ax, "f", x=-0.02, y=1.08)

    ax = fig.add_subplot(gs[3, :])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    example_cards = [
        ("CCND1", "11q13 tumor", "matched positive"),
        ("C1QA", "macrophage", "matched positive"),
        ("JCHAIN", "plasma cell", "matched positive"),
        ("CDH5", "endothelial", "matched positive"),
        ("ERBB2", "CXCL14+ fibroblast", "discordant control"),
    ]
    widths = [0.165, 0.165, 0.165, 0.165, 0.185]
    x0 = 0.025
    for idx, (gene, target, status) in enumerate(example_cards):
        width = widths[idx]
        color = GENE_COLORS[gene]
        face = "#fff8ed" if gene == "ERBB2" else "#ffffff"
        ax.add_patch(plt.Rectangle((x0, 0.26), width, 0.54, facecolor=face, edgecolor=color, linewidth=1.0))
        ax.text(x0 + width / 2, 0.66, gene, ha="center", va="center", fontsize=7.1, fontweight="bold", color=color)
        ax.text(x0 + width / 2, 0.50, target, ha="center", va="center", fontsize=5.4, color=PALETTE["ink"])
        ax.text(x0 + width / 2, 0.36, status, ha="center", va="center", fontsize=5.1, color=PALETTE["muted"])
        if idx < len(example_cards) - 1:
            ax.annotate(
                "",
                xy=(x0 + width + 0.020, 0.53),
                xytext=(x0 + width + 0.004, 0.53),
                arrowprops=dict(arrowstyle="-|>", lw=0.65, color="#9eb3c7"),
            )
        x0 += width + 0.026
    ax.text(
        0.5,
        0.08,
        "These five choices are visualized with true uRNA and cell coordinates in Figure 2.",
        ha="center",
        va="center",
        fontsize=5.8,
        color=PALETTE["ink"],
    )
    ax.set_title("Bridge to the spatial evidence figure", fontsize=7.6, pad=4)
    panel_label(ax, "g", x=-0.02, y=1.08)

    fig.subplots_adjust(left=0.073, right=0.986, top=0.963, bottom=0.052)
    fig.savefig(FIGURES / "figure_1_integrated_urna_evidence_chain.png", bbox_inches="tight")
    fig.savefig(FIGURES / "figure_1_integrated_urna_evidence_chain.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_supplementary_discordance(summary):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), dpi=220)
    df = summary.copy()

    ax = axes[0]
    discord = df.sort_values("sss_spearman_vs_full", ascending=True).head(15).sort_values("sss_spearman_vs_full")
    ax.barh(discord["gene"], discord["sss_spearman_vs_full"], color=PALETTE["red"])
    ax.set_xlabel("Spearman rho")
    ax.set_title("Lowest uRNA/full SSS concordance")
    style_axes(ax)

    ax = axes[1]
    discord["pair"] = discord["full_best_cluster"].astype(str) + "\n-> " + discord["urna_best_cluster"].astype(str)
    ax.scatter(discord["unassigned_fraction"], discord["unassigned_count"], color=PALETTE["red"], s=55, alpha=0.80)
    for _, row in discord.iterrows():
        ax.text(row["unassigned_fraction"], row["unassigned_count"] * 1.03, row["gene"], fontsize=7, color=PALETTE["ink"], ha="center")
    ax.set_yscale("log")
    ax.set_xlabel("uRNA fraction")
    ax.set_ylabel("uRNA count")
    ax.set_title("Discordance is not only low coverage")
    style_axes(ax)

    fig.tight_layout()
    fig.savefig(FIGURES / "supplementary_discordant_genes.png", bbox_inches="tight")
    fig.savefig(FIGURES / "supplementary_discordant_genes.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    setup()
    metadata, counts, summary, urna_sss = load_inputs()
    summary = add_marker_groups(summary)
    overall, cluster_counts, _by_marker, _marker_rows = save_tables(metadata, counts, summary)
    figure_graphical_abstract(metadata, overall)
    figure_combined_lead_in(metadata, counts, summary, overall, cluster_counts)
    figure_supplementary_discordance(summary)


if __name__ == "__main__":
    main()
