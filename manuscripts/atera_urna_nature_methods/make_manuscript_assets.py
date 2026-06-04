from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), dpi=220)

    ax = axes[0, 0]
    labels = ["all rows", "gene rows", "QV>=20 gene", "QV>=20 uRNA"]
    values = [metadata["total_rows"], metadata["gene_rows"], metadata["qv_gene_rows"], metadata["unassigned_rows"]]
    colors = [PALETTE["gray"], "#9eb3c7", PALETTE["blue"], PALETTE["red"]]
    ax.bar(labels, np.array(values) / 1e6, color=colors)
    ax.set_ylabel("transcripts (millions)")
    ax.set_title("Atera WTA transcript scale")
    ax.tick_params(axis="x", rotation=20)
    style_axes(ax)

    ax = axes[0, 1]
    bins = np.linspace(0, counts["unassigned_fraction"].quantile(0.995), 50)
    ax.hist(counts["unassigned_fraction"], bins=bins, color=PALETTE["teal"], alpha=0.85)
    ax.axvline(counts["unassigned_fraction"].median(), color=PALETTE["ink"], linewidth=1.3)
    ax.set_xlabel("uRNA fraction per gene")
    ax.set_ylabel("genes")
    ax.set_title("uRNA fraction across 18,028 genes")
    style_axes(ax)

    ax = axes[1, 0]
    top = counts.sort_values("unassigned_count", ascending=False).head(15).sort_values("unassigned_count")
    ax.barh(top["gene"], top["unassigned_count"] / 1e3, color=PALETTE["gold"])
    ax.set_xlabel("uRNA count (thousands)")
    ax.set_title("Highest-coverage uRNA genes")
    style_axes(ax)

    ax = axes[1, 1]
    cc = cluster_counts.head(10).sort_values("n_genes")
    ax.barh(cc["urna_best_cluster"], cc["n_genes"], color=PALETTE["green"])
    ax.set_xlabel("genes")
    ax.set_title("uRNA-only best COSTE cluster")
    style_axes(ax)

    fig.tight_layout()
    fig.savefig(FIGURES / "figure_1_dataset_and_attribution.png", bbox_inches="tight")
    fig.savefig(FIGURES / "figure_1_dataset_and_attribution.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_concordance(summary, overall):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), dpi=220)
    df = summary.copy()

    ax = axes[0, 0]
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
    ax.set_title("Coverage and concordance")
    style_axes(ax)

    ax = axes[0, 1]
    ax.hist(df["sss_spearman_vs_full"], bins=np.linspace(-0.1, 1.05, 35), color=PALETTE["blue"], alpha=0.85)
    ax.axvline(overall["median_spearman_vs_full"], color=PALETTE["ink"], linewidth=1.3, label="median")
    ax.set_xlabel("Spearman rho")
    ax.set_ylabel("genes")
    ax.set_title("Per-gene SSS concordance")
    style_axes(ax)

    ax = axes[1, 0]
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
    ax.set_title("Concordance thresholds")
    style_axes(ax)

    ax = axes[1, 1]
    grouped = (
        df.groupby("urna_best_module")
        .agg(n=("gene", "size"), median_rho=("sss_spearman_vs_full", "median"))
        .reset_index()
        .sort_values("n", ascending=False)
    )
    ax.scatter(grouped["n"], grouped["median_rho"], s=80, color=PALETTE["red"], alpha=0.85)
    for _, row in grouped.iterrows():
        ax.text(row["n"] + 2, row["median_rho"], row["urna_best_module"], fontsize=8, va="center", color=PALETTE["ink"])
    ax.set_xlabel("genes assigned to module")
    ax.set_ylabel("median rho")
    ax.set_title("Module-level concordance")
    ax.set_xlim(0, max(grouped["n"]) + 60)
    ax.set_ylim(0, 1.05)
    style_axes(ax)

    fig.tight_layout()
    fig.savefig(FIGURES / "figure_2_concordance.png", bbox_inches="tight")
    fig.savefig(FIGURES / "figure_2_concordance.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_marker_validation(summary, urna_sss):
    marker_order = []
    marker_labels = []
    for group, genes in MARKER_GENES.items():
        for gene in genes:
            if gene in urna_sss.index:
                marker_order.append(gene)
                marker_labels.append(f"{gene}  [{group}]")

    heat = urna_sss.loc[marker_order, REPRESENTATIVE_CLUSTERS].copy()
    heat = 1 - heat.clip(0, 1)

    fig, ax = plt.subplots(figsize=(10.5, 9.2), dpi=220)
    image = ax.imshow(heat.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(REPRESENTATIVE_CLUSTERS)))
    ax.set_xticklabels(REPRESENTATIVE_CLUSTERS, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(marker_labels)))
    ax.set_yticklabels(marker_labels, fontsize=7)
    ax.set_title("uRNA-only COSTE proximity for curated markers")
    ax.set_xlabel("cell group")
    ax.set_ylabel("marker gene")
    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("1 - SSS (higher means closer)")
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURES / "figure_3_marker_validation.png", bbox_inches="tight")
    fig.savefig(FIGURES / "figure_3_marker_validation.pdf", bbox_inches="tight")
    plt.close(fig)

    marker_summary = summary[summary["gene"].isin(marker_order)][
        ["gene", "paper_marker_group", "urna_best_cluster", "full_best_cluster", "best_cluster_match", "sss_spearman_vs_full"]
    ].copy()
    marker_summary.to_csv(TABLES / "figure_3_marker_callouts.csv", index=False)


def figure_discordance(summary):
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
    fig.savefig(FIGURES / "figure_4_discordant_genes.png", bbox_inches="tight")
    fig.savefig(FIGURES / "figure_4_discordant_genes.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    setup()
    metadata, counts, summary, urna_sss = load_inputs()
    summary = add_marker_groups(summary)
    overall, cluster_counts, _by_marker, _marker_rows = save_tables(metadata, counts, summary)
    figure_graphical_abstract(metadata, overall)
    figure_overview(metadata, counts, summary, cluster_counts)
    figure_concordance(summary, overall)
    figure_marker_validation(summary, urna_sss)
    figure_discordance(summary)


if __name__ == "__main__":
    main()
