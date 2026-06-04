# Atera WTA breast uRNA segmentation-bias check

## Scope

This analysis follows the BGPT suggestion by treating high-quality unassigned RNAs as independent transcript nodes,
then recomputing COSTE/SSS gene-to-cell profiles and comparing them with the existing all-transcript gene-to-cell SSS table.

## Dataset facts

- Total parquet transcript rows: 740,442,119
- High-quality gene transcript rows: 624,095,990
- High-quality unassigned gene transcript rows: 122,716,153
- Genes with at least one high-quality uRNA: 18,023
- Genes analyzed with uRNA-only COSTE/SSS: 530

## Main readout

- Median per-gene Spearman concordance vs all-transcript SSS: 0.649
- Best-cluster exact match rate vs all-transcript SSS: 0.502

Lower SSS means the uRNA point pattern is spatially closer to that cell group in the COSTE hierarchy.
A concordant uRNA-only profile supports the idea that unassigned transcripts preserve compartment/proximity structure rather than behaving as pure segmentation noise.

## Marker-gene examples

| gene | marker_group | unassigned_count | unassigned_fraction | urna_best_cluster | urna_best_sss | full_best_cluster | sss_spearman_vs_full |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CCND1 | 11q13 | 703247 | 0.205 | 11q13 Invasive Tumor Cells | 0.073 | 11q13 Invasive Tumor Cells | 0.999 |
| ELOVL2 | 11q13 | 26361 | 0.145 | 11q13 Invasive Tumor Cells | 0.044 | 11q13 Invasive Tumor Cells | 1.000 |
| KCNJ3 | 11q13 | 6118 | 0.156 | 11q13 Invasive Tumor Cells | 0.037 | 11q13 Invasive Tumor Cells | 1.000 |
| FGF19 | 11q13 | 1208 | 0.081 | 11q13 Invasive Tumor Cells | 0.187 | 11q13 Invasive Tumor Cells | 0.997 |
| KRT23 | Basal-like | 2847 | 0.170 | Basal-like Structured DCIS Cells | 0.104 | Basal-like Structured DCIS Cells | 1.000 |
| DSC3 | Basal-like | 2037 | 0.161 | Basal-like Structured DCIS Cells | 0.147 | Basal-like Structured DCIS Cells | 0.649 |
| SOSTDC1 | Basal-like | 1680 | 0.171 | Basal-like Structured DCIS Cells | 0.084 | Basal-like Structured DCIS Cells | 1.000 |
| KLK5 | Basal-like | 1514 | 0.174 | Basal-like Structured DCIS Cells | 0.081 | Basal-like Structured DCIS Cells | 1.000 |
| KLK7 | Basal-like | 893 | 0.177 | Basal-like Structured DCIS Cells | 0.099 | Basal-like Structured DCIS Cells | 1.000 |
| ITGB6 | Basal-like | 808 | 0.115 | Basal-like Structured DCIS Cells | 0.099 | Basal-like Structured DCIS Cells | 1.000 |
| MMP7 | Basal-like | 528 | 0.140 | Basal-like Structured DCIS Cells | 0.130 | Basal-like Structured DCIS Cells | 1.000 |
| C3 | Macrophage | 30692 | 0.424 | Macrophages | 0.052 | Macrophages | 0.997 |
| C1QA | Macrophage | 26608 | 0.489 | Macrophages | 0.042 | Macrophages | 1.000 |
| CSF1R | Macrophage | 11023 | 0.323 | Macrophages | 0.039 | Macrophages | 0.698 |
| CD163 | Macrophage | 9289 | 0.346 | Macrophages | 0.047 | Macrophages | 0.698 |
| SIGLEC1 | Macrophage | 1686 | 0.221 | Macrophages | 0.081 | Macrophages | 0.332 |
| IGHA1 | Plasma | 18438 | 0.409 | Plasma Cells | 0.091 | Plasma Cells | 1.000 |
| IGHM | Plasma | 16588 | 0.489 | Plasma Cells | 0.091 | Plasma Cells | 1.000 |
| IGHA2 | Plasma | 12960 | 0.436 | Plasma Cells | 0.082 | Plasma Cells | 1.000 |
| JCHAIN | Plasma | 5925 | 0.353 | Plasma Cells | 0.064 | Plasma Cells | 1.000 |
| PIP | Rare epithelial | 15324 | 0.296 | Luminal-like Amorphous DCIS Cells | 0.260 | Luminal-like Amorphous DCIS Cells | 1.000 |
| HSPB8 | Rare epithelial | 14247 | 0.300 | Luminal-like Amorphous DCIS Cells | 0.174 | Luminal-like Amorphous DCIS Cells | 1.000 |
| CLIC6 | Rare epithelial | 10209 | 0.258 | Luminal-like Amorphous DCIS Cells | 0.200 | Luminal-like Amorphous DCIS Cells | 1.000 |
| TAT | Rare epithelial | 8377 | 0.219 | Apocrine Cells | 0.158 | Apocrine Cells | 1.000 |
| FOXA1 | Tumor shifted to basal/DCIS | 148561 | 0.207 | Basal-like Structured DCIS Cells | 0.119 | Basal-like Structured DCIS Cells | 0.649 |
| CTTN | Tumor shifted to basal/DCIS | 48369 | 0.160 | Basal-like Structured DCIS Cells | 0.108 | Basal-like Structured DCIS Cells | 1.000 |
| ANO1 | Tumor shifted to basal/DCIS | 47247 | 0.146 | Basal-like Structured DCIS Cells | 0.089 | Basal-like Structured DCIS Cells | 1.000 |
| ESR1 | Tumor shifted to basal/DCIS | 24679 | 0.145 | Basal-like Structured DCIS Cells | 0.116 | Basal-like Structured DCIS Cells | 1.000 |
| FADD | Tumor shifted to basal/DCIS | 10152 | 0.214 | Basal-like Structured DCIS Cells | 0.119 | Basal-like Structured DCIS Cells | 1.000 |
| PGR | Tumor shifted to basal/DCIS | 6300 | 0.156 | Basal-like Structured DCIS Cells | 0.104 | Basal-like Structured DCIS Cells | 0.649 |
| ERBB2 | Tumor shifted to basal/DCIS | 1947 | 0.184 | CXCL14+ Fibroblasts | 0.166 | Basal-like Structured DCIS Cells | 0.227 |
| FGF3 | Tumor shifted to basal/DCIS | 326 | 0.140 | 11q13 Invasive Tumor Cells | 0.184 | Basal-like Structured DCIS Cells | 0.869 |
| EPAS1 | Vascular | 64127 | 0.349 | Pericytes | 0.034 | Pericytes | 1.000 |
| CDH5 | Vascular | 6827 | 0.242 | Endothelial Cells | 0.042 | Endothelial Cells | 1.000 |
| MMRN2 | Vascular | 5475 | 0.247 | Endothelial Cells | 0.040 | Endothelial Cells | 1.000 |

## Output files

- `transcript_gene_counts.csv`: assigned/uRNA counts and uRNA fraction by gene.
- `urna_only_t_and_c_result.csv`: uRNA-only gene-to-cell SSS matrix.
- `urna_segmentation_bias_gene_summary.csv`: best cluster/module and all-transcript concordance summary.
- `marker_gene_urna_sss_heatmap.*`, `full_vs_urna_sss_scatter.*`, `per_gene_sss_concordance_histogram.*`: summary figures.
