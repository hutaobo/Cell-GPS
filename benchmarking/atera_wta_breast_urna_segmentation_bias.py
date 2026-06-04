#!/usr/bin/env python
"""uRNA-only COSTE/SSS segmentation-bias check for Atera WTA breast data."""

from __future__ import annotations

import argparse
import json
import logging
import math
import time
from collections import Counter
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq
import seaborn as sns
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors


LOGGER = logging.getLogger("atera_urna_bias")


DEFAULT_MARKER_GROUPS: dict[str, list[str]] = {
    "11q13": ["KCNJ3", "ELOVL2", "CCND1", "FGF19"],
    "Tumor shifted to basal/DCIS": ["ANO1", "FGF3", "CTTN", "FADD", "ERBB2", "ESR1", "PGR", "FOXA1"],
    "Basal-like": ["SOSTDC1", "KLK5", "ITGB6", "DSC3", "KRT23", "KLK7", "MMP7"],
    "Macrophage": ["C1QA", "CD163", "CSF1R", "SIGLEC1", "C3"],
    "Plasma": ["IGHA1", "IGHA2", "JCHAIN", "IGHM"],
    "Vascular": ["CDH5", "MMRN2", "EPAS1"],
    "Rare epithelial": ["TAT", "HSPB8", "CLIC6", "PIP"],
}


MODULE_GROUPS: dict[str, list[str]] = {
    "Tumor": [
        "11q13 Invasive Tumor Cells",
        "11q13 Invasive Tumor Cells (G1/S)",
        "11q13 Invasive Tumor Cells (Mitotic)",
        "CAFs, Invasive Associated",
    ],
    "Vascular/Stromal": [
        "Endothelial Cells",
        "Pericytes",
        "CXCL14+ Fibroblasts",
        "CAFs, DCIS Associated",
        "Myoepithelial Cells",
    ],
    "Immune": [
        "B Cells",
        "T Lymphocytes",
        "Dendritic Cells",
        "Macrophages",
        "Plasma Cells",
        "Myeloid Cells",
        "Mast Cells",
    ],
    "Rare epithelial": [
        "Apocrine Cells",
        "Luminal-like Amorphous DCIS Cells",
        "Basal-like Structured DCIS Cells",
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transcripts", type=Path, required=True, help="Atera transcripts.parquet path.")
    parser.add_argument("--cells", type=Path, required=True, help="Atera cells.csv.gz or cells.parquet path.")
    parser.add_argument("--cell-groups", type=Path, required=True, help="Atera cell group CSV path.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--full-tbc", type=Path, default=None, help="Existing all-transcript t_and_c CSV for comparison.")
    parser.add_argument("--qv-min", type=float, default=20.0)
    parser.add_argument("--min-urna-count", type=int, default=100)
    parser.add_argument("--top-genes", type=int, default=500)
    parser.add_argument("--extra-genes", nargs="*", default=[])
    parser.add_argument("--force", action="store_true", help="Ignore cached count/coordinate tables.")
    return parser.parse_args()


def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s")


def marker_genes() -> list[str]:
    genes: list[str] = []
    for group in DEFAULT_MARKER_GROUPS.values():
        genes.extend(group)
    return sorted(set(genes))


def marker_group_lookup() -> dict[str, str]:
    lookup: dict[str, str] = {}
    for group_name, genes in DEFAULT_MARKER_GROUPS.items():
        for gene in genes:
            lookup[gene] = group_name
    return lookup


def load_cells(cells_path: Path, cell_groups_path: Path) -> pd.DataFrame:
    if cells_path.suffix == ".parquet":
        cells = pd.read_parquet(cells_path)
    else:
        cells = pd.read_csv(cells_path)

    groups = pd.read_csv(cell_groups_path)
    rename_map = {}
    if "Barcode" in groups.columns:
        rename_map["Barcode"] = "cell_id"
    if "Clusters" in groups.columns:
        rename_map["Clusters"] = "celltype"
    if "group" in groups.columns:
        rename_map["group"] = "celltype"
    groups = groups.rename(columns=rename_map)
    required_groups = {"cell_id", "celltype"}
    if not required_groups.issubset(groups.columns):
        raise ValueError(f"Cell group file is missing {required_groups - set(groups.columns)}")

    coord_rename = {}
    if "x_centroid" in cells.columns:
        coord_rename["x_centroid"] = "x"
    if "y_centroid" in cells.columns:
        coord_rename["y_centroid"] = "y"
    cells = cells.rename(columns=coord_rename)
    required_cells = {"cell_id", "x", "y"}
    if not required_cells.issubset(cells.columns):
        raise ValueError(f"Cells file is missing {required_cells - set(cells.columns)}")

    merged = cells[["cell_id", "x", "y"]].merge(groups[["cell_id", "celltype"]], on="cell_id", how="left")
    merged = merged.dropna(subset=["celltype", "x", "y"]).copy()
    merged["celltype"] = merged["celltype"].astype(str)
    merged = merged.loc[merged["celltype"] != "Unassigned"].reset_index(drop=True)
    return merged


def value_counts_feature_name(table) -> Counter:
    counts: Counter = Counter()
    if table.num_rows == 0:
        return counts
    for item in pc.value_counts(table["feature_name"]).to_pylist():
        counts[item["values"]] += item["counts"]
    return counts


def count_transcripts_by_gene(transcripts_path: Path, qv_min: float, output_dir: Path, force: bool) -> pd.DataFrame:
    cache_path = output_dir / "transcript_gene_counts.csv"
    if cache_path.exists() and not force:
        LOGGER.info("Reusing cached gene counts: %s", cache_path)
        return pd.read_csv(cache_path)

    parquet = pq.ParquetFile(transcripts_path)
    total_rows = 0
    gene_rows = 0
    qv_gene_rows = 0
    unassigned_rows = 0
    assigned_counts: Counter = Counter()
    unassigned_counts: Counter = Counter()
    start = time.time()

    columns = ["cell_id", "feature_name", "qv", "is_gene"]
    for row_group in range(parquet.metadata.num_row_groups):
        table = parquet.read_row_group(row_group, columns=columns)
        total_rows += table.num_rows
        gene_mask = pc.equal(table["is_gene"], True)
        gene_table = table.filter(gene_mask)
        gene_rows += gene_table.num_rows
        qv_gene_table = gene_table.filter(pc.greater_equal(gene_table["qv"], qv_min))
        qv_gene_rows += qv_gene_table.num_rows
        unassigned_table = qv_gene_table.filter(pc.equal(qv_gene_table["cell_id"], "UNASSIGNED"))
        assigned_table = qv_gene_table.filter(pc.not_equal(qv_gene_table["cell_id"], "UNASSIGNED"))
        unassigned_rows += unassigned_table.num_rows
        unassigned_counts.update(value_counts_feature_name(unassigned_table))
        assigned_counts.update(value_counts_feature_name(assigned_table))

        if (row_group + 1) % 100 == 0:
            LOGGER.info(
                "Counted row groups %s/%s: %.1fM uRNA rows, %s uRNA genes, %.1fs",
                row_group + 1,
                parquet.metadata.num_row_groups,
                unassigned_rows / 1e6,
                len(unassigned_counts),
                time.time() - start,
            )

    genes = sorted(set(assigned_counts) | set(unassigned_counts))
    rows = []
    for gene in genes:
        assigned = int(assigned_counts.get(gene, 0))
        unassigned = int(unassigned_counts.get(gene, 0))
        total_qv_gene = assigned + unassigned
        rows.append(
            {
                "gene": gene,
                "assigned_count": assigned,
                "unassigned_count": unassigned,
                "total_qv_gene_count": total_qv_gene,
                "unassigned_fraction": unassigned / total_qv_gene if total_qv_gene else np.nan,
            }
        )
    counts = pd.DataFrame(rows).sort_values("unassigned_count", ascending=False)
    counts.attrs["total_rows"] = total_rows
    counts.attrs["gene_rows"] = gene_rows
    counts.attrs["qv_gene_rows"] = qv_gene_rows
    counts.attrs["unassigned_rows"] = unassigned_rows
    counts.to_csv(cache_path, index=False)
    (output_dir / "transcript_gene_count_summary.json").write_text(
        json.dumps(
            {
                "total_rows": total_rows,
                "gene_rows": gene_rows,
                "qv_gene_rows": qv_gene_rows,
                "unassigned_rows": unassigned_rows,
                "n_genes": len(genes),
                "n_genes_with_unassigned": int((counts["unassigned_count"] > 0).sum()),
                "elapsed_seconds": time.time() - start,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return counts


def choose_genes(counts: pd.DataFrame, top_genes: int, min_count: int, extra_genes: Iterable[str]) -> list[str]:
    markers = set(marker_genes())
    extras = {gene.strip() for gene in extra_genes if gene.strip()}
    eligible = counts.loc[counts["unassigned_count"] >= min_count, "gene"].tolist()
    selected = set(eligible[:top_genes]) | markers | extras
    available = set(counts.loc[counts["unassigned_count"] > 0, "gene"])
    return sorted(selected & available)


def load_selected_urna_coords(
    transcripts_path: Path,
    genes: list[str],
    qv_min: float,
    output_dir: Path,
    force: bool,
) -> pd.DataFrame:
    cache_path = output_dir / f"selected_urna_coords_{len(genes)}genes.parquet"
    if cache_path.exists() and not force:
        LOGGER.info("Reusing cached selected uRNA coordinates: %s", cache_path)
        return pd.read_parquet(cache_path)

    gene_set = set(genes)
    parquet = pq.ParquetFile(transcripts_path)
    frames: list[pd.DataFrame] = []
    start = time.time()
    columns = ["cell_id", "feature_name", "x_location", "y_location", "qv", "is_gene"]
    gene_array = pa_array(list(gene_set))

    for row_group in range(parquet.metadata.num_row_groups):
        table = parquet.read_row_group(row_group, columns=columns)
        mask = pc.and_(
            pc.equal(table["cell_id"], "UNASSIGNED"),
            pc.and_(
                pc.greater_equal(table["qv"], qv_min),
                pc.and_(pc.equal(table["is_gene"], True), pc.is_in(table["feature_name"], value_set=gene_array)),
            ),
        )
        sub = table.filter(mask)
        if sub.num_rows:
            frame = sub.select(["feature_name", "x_location", "y_location"]).to_pandas()
            frame = frame.rename(columns={"x_location": "x", "y_location": "y"})
            frames.append(frame)
        if (row_group + 1) % 100 == 0:
            kept = sum(len(frame) for frame in frames)
            LOGGER.info(
                "Loaded row groups %s/%s: %.1fM selected uRNA rows, %.1fs",
                row_group + 1,
                parquet.metadata.num_row_groups,
                kept / 1e6,
                time.time() - start,
            )

    coords = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["feature_name", "x", "y"])
    coords["feature_name"] = coords["feature_name"].astype(str)
    coords["x"] = coords["x"].astype("float32")
    coords["y"] = coords["y"].astype("float32")
    coords.to_parquet(cache_path, index=False)
    return coords


def pa_array(values: list[str]):
    import pyarrow as pa

    return pa.array(values)


class GeneSssComputer:
    def __init__(self, cells: pd.DataFrame):
        self.cells = cells.copy()
        self.celltypes = sorted(self.cells["celltype"].unique())
        self.codes = pd.Categorical(self.cells["celltype"], categories=self.celltypes).codes
        self.cell_coords = self.cells[["x", "y"]].to_numpy(dtype=np.float32)
        self.cell_trees: dict[str, NearestNeighbors] = {}
        self.base_distance = self._compute_cell_distance_matrix()

    def _compute_cell_distance_matrix(self) -> pd.DataFrame:
        LOGGER.info("Fitting nearest-neighbor trees for %s cell types", len(self.celltypes))
        n_types = len(self.celltypes)
        matrix = np.full((n_types, n_types), np.nan, dtype=float)
        for target_idx, target in enumerate(self.celltypes):
            target_coords = self.cell_coords[self.codes == target_idx]
            tree = NearestNeighbors(n_neighbors=1, algorithm="auto").fit(target_coords)
            self.cell_trees[target] = tree
            distances = tree.kneighbors(self.cell_coords, return_distance=True)[0][:, 0]
            sums = np.bincount(self.codes, weights=distances, minlength=n_types)
            counts = np.bincount(self.codes, minlength=n_types)
            matrix[:, target_idx] = np.divide(sums, counts, out=np.full(n_types, np.nan), where=counts > 0)
        return pd.DataFrame(matrix, index=self.celltypes, columns=self.celltypes)

    def compute_gene_row(self, gene: str, gene_coords: np.ndarray) -> pd.Series:
        labels = self.celltypes + [gene]
        n_types = len(self.celltypes)
        matrix = np.zeros((n_types + 1, n_types + 1), dtype=float)
        matrix[:n_types, :n_types] = self.base_distance.loc[self.celltypes, self.celltypes].to_numpy()

        for target_idx, target in enumerate(self.celltypes):
            distances = self.cell_trees[target].kneighbors(gene_coords, return_distance=True)[0][:, 0]
            matrix[n_types, target_idx] = float(np.mean(distances))

        gene_tree = NearestNeighbors(n_neighbors=1, algorithm="auto").fit(gene_coords)
        cell_to_gene = gene_tree.kneighbors(self.cell_coords, return_distance=True)[0][:, 0]
        sums = np.bincount(self.codes, weights=cell_to_gene, minlength=n_types)
        counts = np.bincount(self.codes, minlength=n_types)
        matrix[:n_types, n_types] = np.divide(sums, counts, out=np.full(n_types, np.nan), where=counts > 0)
        matrix[n_types, n_types] = 0.0

        distance_matrix = pd.DataFrame(matrix, index=labels, columns=labels)
        row_coph = compute_row_cophenetic(distance_matrix)
        return row_coph.loc[gene, self.celltypes]


def compute_row_cophenetic(distance_matrix: pd.DataFrame, method: str = "average") -> pd.DataFrame:
    values = distance_matrix.to_numpy(dtype=float)
    row_linkage = linkage(values, method=method)
    _, row_coph_condensed = cophenet(row_linkage, pdist(values))
    square = squareform(row_coph_condensed)
    row_coph = pd.DataFrame(square, index=distance_matrix.index, columns=distance_matrix.index)
    dmin = float(np.nanmin(square))
    dmax = float(np.nanmax(square))
    if math.isclose(dmin, dmax):
        return row_coph
    return (row_coph - dmin) / (dmax - dmin)


def compute_urna_sss(cells: pd.DataFrame, coords: pd.DataFrame, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    computer = GeneSssComputer(cells)
    computer.base_distance.to_csv(output_dir / "celltype_directed_distance_matrix.csv")
    compute_row_cophenetic(computer.base_distance).to_csv(output_dir / "celltype_structuremap_from_directed_distances.csv")

    rows: list[pd.Series] = []
    start = time.time()
    grouped = coords.groupby("feature_name", sort=True)
    for idx, (gene, group) in enumerate(grouped, start=1):
        gene_coords = group[["x", "y"]].to_numpy(dtype=np.float32)
        row = computer.compute_gene_row(gene, gene_coords)
        row.name = gene
        rows.append(row)
        if idx % 25 == 0:
            LOGGER.info("Computed uRNA SSS for %s/%s genes in %.1fs", idx, len(grouped), time.time() - start)

    sss = pd.DataFrame(rows)
    sss.to_csv(output_dir / "urna_only_t_and_c_result.csv")
    return sss, computer.base_distance


def summarize_gene_assignments(
    urna_sss: pd.DataFrame,
    counts: pd.DataFrame,
    full_tbc: pd.DataFrame | None,
) -> pd.DataFrame:
    count_lookup = counts.set_index("gene")
    marker_lookup = marker_group_lookup()
    rows = []
    for gene, row in urna_sss.iterrows():
        best_cluster = row.astype(float).idxmin()
        best_sss = float(row[best_cluster])
        module_scores = {}
        for module, members in MODULE_GROUPS.items():
            present = [member for member in members if member in row.index]
            module_scores[module] = float(row[present].mean()) if present else np.nan
        best_module = min(module_scores, key=lambda key: module_scores[key] if not np.isnan(module_scores[key]) else np.inf)

        full_best_cluster = None
        full_best_sss = np.nan
        spearman = np.nan
        if full_tbc is not None and gene in full_tbc.index:
            common = [col for col in urna_sss.columns if col in full_tbc.columns]
            full_row = full_tbc.loc[gene, common].astype(float)
            full_best_cluster = full_row.idxmin()
            full_best_sss = float(full_row[full_best_cluster])
            rho = spearmanr(row[common].astype(float), full_row, nan_policy="omit")
            spearman = float(rho.statistic)

        counts_row = count_lookup.loc[gene] if gene in count_lookup.index else None
        rows.append(
            {
                "gene": gene,
                "marker_group": marker_lookup.get(gene, ""),
                "unassigned_count": int(counts_row["unassigned_count"]) if counts_row is not None else np.nan,
                "assigned_count": int(counts_row["assigned_count"]) if counts_row is not None else np.nan,
                "unassigned_fraction": float(counts_row["unassigned_fraction"]) if counts_row is not None else np.nan,
                "urna_best_cluster": best_cluster,
                "urna_best_sss": best_sss,
                "urna_best_module": best_module,
                "full_best_cluster": full_best_cluster,
                "full_best_sss": full_best_sss,
                "best_cluster_match": bool(best_cluster == full_best_cluster) if full_best_cluster is not None else np.nan,
                "sss_spearman_vs_full": spearman,
                **{f"urna_module_mean_{module}": value for module, value in module_scores.items()},
            }
        )
    return pd.DataFrame(rows).sort_values(["marker_group", "unassigned_count"], ascending=[True, False])


def make_figures(
    urna_sss: pd.DataFrame,
    summary: pd.DataFrame,
    counts: pd.DataFrame,
    output_dir: Path,
    full_tbc: pd.DataFrame | None,
) -> None:
    sns.set_theme(style="white", context="paper")
    marker_order = [gene for gene in marker_genes() if gene in urna_sss.index]
    if marker_order:
        fig, ax = plt.subplots(figsize=(12, max(6, 0.25 * len(marker_order))))
        sns.heatmap(urna_sss.loc[marker_order].astype(float), cmap="vlag", vmin=0, vmax=1, ax=ax)
        ax.set_title("uRNA-only COSTE SSS: marker genes vs breast cell groups")
        ax.set_xlabel("Cell group")
        ax.set_ylabel("Marker gene")
        fig.tight_layout()
        fig.savefig(output_dir / "marker_gene_urna_sss_heatmap.png", dpi=300)
        fig.savefig(output_dir / "marker_gene_urna_sss_heatmap.pdf")
        plt.close(fig)

    top_counts = counts.head(40).copy()
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.barplot(data=top_counts, y="gene", x="unassigned_count", color="#4C78A8", ax=ax)
    ax.set_title("Top genes by high-quality unassigned RNA count")
    ax.set_xlabel("uRNA count")
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(output_dir / "top_urna_gene_counts.png", dpi=300)
    fig.savefig(output_dir / "top_urna_gene_counts.pdf")
    plt.close(fig)

    if full_tbc is not None:
        common_genes = [gene for gene in urna_sss.index if gene in full_tbc.index]
        common_cols = [col for col in urna_sss.columns if col in full_tbc.columns]
        if common_genes and common_cols:
            x = full_tbc.loc[common_genes, common_cols].astype(float).to_numpy().ravel()
            y = urna_sss.loc[common_genes, common_cols].astype(float).to_numpy().ravel()
            fig, ax = plt.subplots(figsize=(6.5, 6.5))
            ax.scatter(x, y, s=4, alpha=0.2, linewidths=0)
            ax.plot([0, 1], [0, 1], color="black", lw=1)
            rho = spearmanr(x, y, nan_policy="omit").statistic
            ax.set_title(f"All-transcript vs uRNA-only SSS (Spearman {rho:.3f})")
            ax.set_xlabel("All-transcript gene-to-cell SSS")
            ax.set_ylabel("uRNA-only gene-to-cell SSS")
            fig.tight_layout()
            fig.savefig(output_dir / "full_vs_urna_sss_scatter.png", dpi=300)
            fig.savefig(output_dir / "full_vs_urna_sss_scatter.pdf")
            plt.close(fig)

            fig, ax = plt.subplots(figsize=(7, 4))
            sns.histplot(summary["sss_spearman_vs_full"].dropna(), bins=40, color="#59A14F", ax=ax)
            ax.set_title("Per-gene SSS profile concordance")
            ax.set_xlabel("Spearman correlation: uRNA-only vs all-transcript SSS")
            fig.tight_layout()
            fig.savefig(output_dir / "per_gene_sss_concordance_histogram.png", dpi=300)
            fig.savefig(output_dir / "per_gene_sss_concordance_histogram.pdf")
            plt.close(fig)


def write_markdown_report(
    output_dir: Path,
    counts: pd.DataFrame,
    urna_sss: pd.DataFrame,
    summary: pd.DataFrame,
    full_tbc: pd.DataFrame | None,
    metadata: dict,
) -> None:
    marker_summary = summary.loc[summary["marker_group"] != ""].copy()
    lines = [
        "# Atera WTA breast uRNA segmentation-bias check",
        "",
        "## Scope",
        "",
        "This analysis follows the BGPT suggestion by treating high-quality unassigned RNAs as independent transcript nodes,",
        "then recomputing COSTE/SSS gene-to-cell profiles and comparing them with the existing all-transcript gene-to-cell SSS table.",
        "",
        "## Dataset facts",
        "",
        f"- Total parquet transcript rows: {metadata.get('total_rows', 'NA'):,}",
        f"- High-quality gene transcript rows: {metadata.get('qv_gene_rows', 'NA'):,}",
        f"- High-quality unassigned gene transcript rows: {metadata.get('unassigned_rows', 'NA'):,}",
        f"- Genes with at least one high-quality uRNA: {int((counts['unassigned_count'] > 0).sum()):,}",
        f"- Genes analyzed with uRNA-only COSTE/SSS: {len(urna_sss):,}",
        "",
        "## Main readout",
        "",
    ]
    if full_tbc is not None and summary["sss_spearman_vs_full"].notna().any():
        lines.extend(
            [
                f"- Median per-gene Spearman concordance vs all-transcript SSS: {summary['sss_spearman_vs_full'].median():.3f}",
                f"- Best-cluster exact match rate vs all-transcript SSS: {summary['best_cluster_match'].mean():.3f}",
                "",
            ]
        )
    lines.extend(
        [
            "Lower SSS means the uRNA point pattern is spatially closer to that cell group in the COSTE hierarchy.",
            "A concordant uRNA-only profile supports the idea that unassigned transcripts preserve compartment/proximity structure rather than behaving as pure segmentation noise.",
            "",
            "## Marker-gene examples",
            "",
        ]
    )
    if not marker_summary.empty:
        cols = [
            "gene",
            "marker_group",
            "unassigned_count",
            "unassigned_fraction",
            "urna_best_cluster",
            "urna_best_sss",
            "full_best_cluster",
            "sss_spearman_vs_full",
        ]
        lines.append(markdown_table(marker_summary[cols]))
    else:
        lines.append("No marker genes passed the uRNA coordinate filter.")
    lines.extend(
        [
            "",
            "## Output files",
            "",
            "- `transcript_gene_counts.csv`: assigned/uRNA counts and uRNA fraction by gene.",
            "- `urna_only_t_and_c_result.csv`: uRNA-only gene-to-cell SSS matrix.",
            "- `urna_segmentation_bias_gene_summary.csv`: best cluster/module and all-transcript concordance summary.",
            "- `marker_gene_urna_sss_heatmap.*`, `full_vs_urna_sss_scatter.*`, `per_gene_sss_concordance_histogram.*`: summary figures.",
            "",
        ]
    )
    (output_dir / "Atera_WTA_breast_uRNA_segmentation_bias_report.md").write_text("\n".join(lines), encoding="utf-8")


def markdown_table(df: pd.DataFrame) -> str:
    def fmt(value) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, (float, np.floating)):
            return f"{value:.3f}"
        return str(value)

    headers = list(df.columns)
    rows = [[fmt(value) for value in row] for row in df.to_numpy()]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(cell.replace("|", "/") for cell in row) + " |")
    return "\n".join(lines)


def main() -> None:
    configure_logging()
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Loading cells and annotations")
    cells = load_cells(args.cells, args.cell_groups)
    LOGGER.info("Loaded %s annotated cells across %s cell groups", len(cells), cells["celltype"].nunique())

    counts = count_transcripts_by_gene(args.transcripts, args.qv_min, args.output_dir, args.force)
    selected_genes = choose_genes(counts, args.top_genes, args.min_urna_count, args.extra_genes)
    pd.Series(selected_genes, name="gene").to_csv(args.output_dir / "selected_urna_genes.csv", index=False)
    LOGGER.info("Selected %s genes for uRNA-only SSS", len(selected_genes))

    coords = load_selected_urna_coords(args.transcripts, selected_genes, args.qv_min, args.output_dir, args.force)
    LOGGER.info("Loaded %.1fM selected uRNA coordinate rows", len(coords) / 1e6)
    urna_sss, _ = compute_urna_sss(cells, coords, args.output_dir)

    full_tbc = None
    if args.full_tbc is not None and args.full_tbc.exists():
        full_tbc = pd.read_csv(args.full_tbc, index_col=0)
        full_tbc = full_tbc.loc[:, [col for col in urna_sss.columns if col in full_tbc.columns]]
        LOGGER.info("Loaded all-transcript t_and_c reference: %s genes x %s cell groups", *full_tbc.shape)
    elif args.full_tbc is not None:
        LOGGER.warning("All-transcript t_and_c reference not found: %s", args.full_tbc)

    summary = summarize_gene_assignments(urna_sss, counts, full_tbc)
    summary.to_csv(args.output_dir / "urna_segmentation_bias_gene_summary.csv", index=False)
    make_figures(urna_sss, summary, counts, args.output_dir, full_tbc)

    meta_path = args.output_dir / "transcript_gene_count_summary.json"
    metadata = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.exists() else {}
    metadata.update(
        {
            "cells": str(args.cells),
            "cell_groups": str(args.cell_groups),
            "transcripts": str(args.transcripts),
            "full_tbc": str(args.full_tbc) if args.full_tbc else None,
            "qv_min": args.qv_min,
            "min_urna_count": args.min_urna_count,
            "top_genes": args.top_genes,
            "n_selected_genes": len(selected_genes),
            "n_urna_sss_genes": len(urna_sss),
        }
    )
    (args.output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    write_markdown_report(args.output_dir, counts, urna_sss, summary, full_tbc, metadata)
    LOGGER.info("Done: %s", args.output_dir)


if __name__ == "__main__":
    main()
