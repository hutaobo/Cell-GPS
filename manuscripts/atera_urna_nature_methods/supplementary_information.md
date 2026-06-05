# Supplementary Information

## Supplementary Note 1: Scope of the current analysis

This study was designed as an orthogonal segmentation-bias check for the public 10x Genomics Atera whole-transcriptome FFPE breast cancer dataset. The analysis does not attempt to reassign unassigned RNAs to cells. Instead, high-quality unassigned RNAs are retained as spatial molecular observations and are compared with annotated cell-group coordinates using COSTE-derived Spatial Separation Score profiles.

The current draft should be read as a single-dataset proof-of-concept. A final Nature Methods Article should add multi-dataset validation and comparisons to transcript-aware segmentation or reassignment methods.

## Supplementary Note 2: Dataset scale

The source transcript table contained 740,442,119 rows. Of these, 738,289,569 were gene rows and 624,095,990 were high-quality gene transcripts after applying QV >= 20. The high-quality unassigned gene transcript pool contained 122,716,153 molecules, equal to 19.7% of high-quality gene transcripts.

uRNAs were present for 18,023 of 18,028 genes. The median gene-level uRNA fraction across the full panel was 0.160. The analyzed set consisted of the top 500 genes by uRNA count plus curated marker genes, yielding 530 genes for uRNA-only COSTE analysis.

## Supplementary Note 3: Concordance with all-transcript COSTE profiles

Across the 530 analyzed genes, the median Spearman concordance between uRNA-only and all-transcript gene-to-cell SSS profiles was 0.649, with an interquartile range of 0.226 to 0.998. The mean concordance was 0.555. Exact best-cell-group agreement was observed for 50.2% of genes.

The broad distribution is informative. High-concordance genes show that uRNAs can recover the same major tissue compartments as all-transcript profiles. Lower-concordance genes indicate that uRNAs are not merely a diluted copy of cell-assigned signal and may emphasize cell boundaries, extracellular or extrasomatic compartments, highly abundant stromal regions, or segmentation-sensitive interfaces.

## Supplementary Note 4: Marker groups

Curated markers were used as positive controls because their expected tissue compartments are interpretable in the Atera breast cancer annotation.

The 11q13 marker group consisted of CCND1, ELOVL2, KCNJ3 and FGF19. All four genes were closest to 11q13 invasive tumor cells in uRNA-only COSTE profiles, with a median Spearman concordance of 0.999 against all-transcript profiles.

Basal-like DCIS markers KRT23, DSC3, SOSTDC1, KLK5, KLK7, ITGB6 and MMP7 all localized to basal-like structured DCIS cells. Macrophage markers C3, C1QA, CSF1R, CD163 and SIGLEC1 localized to macrophages. Plasma-cell markers IGHA1, IGHM, IGHA2 and JCHAIN localized to plasma cells. Vascular markers EPAS1, CDH5 and MMRN2 localized to pericyte or endothelial compartments.

Tumor and DCIS interface markers were more heterogeneous. FOXA1, CTTN, ANO1, ESR1, FADD and PGR localized to basal-like structured DCIS cells, whereas FGF3 localized to 11q13 invasive tumor cells and ERBB2 shifted toward CXCL14+ fibroblast proximity in the uRNA-only analysis. These shifted cases are useful candidates for spatial inspection because they may mark boundary-sensitive or stromal-adjacent signal rather than pure cell-intrinsic expression.

## Supplementary Note 5: Proposed figures

Graphical Abstract. Workflow showing extraction of high-quality uRNAs, uRNA-only COSTE and orthogonal comparison to all-transcript SSS profiles.

Figure 1. Dataset scale and uRNA-only cell-group attribution. Panel A shows transcript counts through filtering. Panel B shows the distribution of gene-level uRNA fractions. Panel C shows genes with the highest uRNA counts. Panel D shows the uRNA-only best COSTE cell group for analyzed genes.

Figure 2. Concordance between uRNA-only and all-transcript COSTE profiles. Panel A relates uRNA count to per-gene Spearman concordance. Panel B shows the concordance distribution. Panel C summarizes concordance thresholds. Panel D summarizes module-level concordance.

Figure 3. Marker validation heatmap. Curated marker genes are shown against representative cell groups using 1 - SSS, where higher values indicate closer COSTE proximity.

Figure 4. Lowest-concordance genes. The figure highlights genes for which uRNA-only and all-transcript profiles diverge, showing that discordance is not only explained by low uRNA coverage.

Figure 5. Real spatial examples. The figure combines a whole-section overview, true local uRNA/cell-coordinate zoom maps, SSS heatmaps, concordance context, marker-validation statistics and per-example calls. CCND1, C1QA, JCHAIN and CDH5 are concordant positive examples; ERBB2 is a discordant example whose uRNA-only best group differs from the all-transcript best group.

## Supplementary Tables

Supplementary Table 1 (`tables/summary_statistics.csv`). Overall transcript, gene coverage and concordance statistics.

Supplementary Table 2 (`tables/celltype_attribution_counts.csv`). uRNA-only best COSTE cell-group counts for analyzed genes.

Supplementary Table 3 (`tables/marker_group_summary.csv`). Marker group-level uRNA counts, fractions, concordance and best-cluster match rates.

Supplementary Table 4 (`tables/marker_validation_table.csv`). Per-marker uRNA-only and all-transcript best cell groups.

Supplementary Table 5 (`tables/lowest_concordance_genes.csv`). The 30 lowest-concordance genes.

Supplementary Table 6 (`tables/top_50_uRNA_genes.csv`). The 50 genes with the highest high-quality uRNA counts.

## Supplementary Note 6: Limitations

This analysis does not prove that every unassigned RNA is biologically meaningful. It shows that the uRNA pool contains structured spatial information that survives stringent quality filtering and can recover expected tissue compartments for marker genes.

The analysis uses a single public dataset, one cell-group annotation set and one all-transcript COSTE reference table. It does not test all possible segmentation algorithms, segmentation parameter settings, tissue types or imaging platforms. It also does not classify uRNAs into technical, extracellular, extrasomatic, diffusion-related or missed-cell categories.

The next experimental step is to repeat the same diagnostic on multiple tissues and platforms, then add synthetic segmentation perturbations and direct comparisons with transcript-aware cell segmentation or reassignment methods.
