# Unassigned RNAs report segmentation-orthogonal tissue architecture in whole-transcriptome imaging data

Taobo Hu, Mengping Long, Mats Nilsson and colleagues

Affiliations: [to be completed]

Correspondence: [to be completed]

## Abstract

Image-based spatial transcriptomics relies on assigning detected RNA molecules to segmented cells, yet large fractions of high-quality transcripts remain outside cell masks and are commonly discarded. We tested whether unassigned RNAs (uRNAs) can provide an orthogonal readout of tissue architecture by applying Cophenetic Spatial Topology Embedding (COSTE) to uRNA coordinates from a 10x Genomics Atera whole-transcriptome FFPE breast cancer dataset. The dataset contained 624,095,990 high-quality gene transcripts, of which 122,716,153 were unassigned. uRNAs were detected for 18,023 of 18,028 genes. Across 530 high-coverage or marker genes, uRNA-only COSTE profiles showed moderate concordance with all-transcript profiles and recovered expected tumor, stromal, immune and vascular compartments. Canonical markers for 11q13-amplified tumor cells, basal-like DCIS, macrophages, plasma cells and endothelial/pericyte structures retained cell-type-specific spatial proximity. uRNA-only topology therefore functions as a segmentation-orthogonal diagnostic for structured signal and boundary-sensitive biology.

Spatial transcriptomics has shifted tissue analysis from dissociated cell profiles to molecular maps that preserve anatomical context [1]. Modern platforms combine high molecular plexity, subcellular coordinates and large tissue fields, creating data structures in which the primary observation is an RNA molecule with a physical location rather than a cell with a count vector. Community tools such as Giotto and Squidpy have made spatial analysis scalable and accessible [2,3], but many downstream analyses still require a cell-level matrix constructed by assigning molecules to segmented cells. This dependency creates a practical vulnerability: when segmentation is uncertain, the biological interpretation of cell states, neighborhoods and ligand-receptor relationships can change.

Cell segmentation in imaging-based spatial transcriptomics is difficult because tissue sections are two-dimensional cuts through three-dimensional objects, cell bodies can overlap, nuclei can be absent from the section plane and transcript distributions can extend beyond visible cell masks. Transcript-aware segmentation methods have addressed parts of this problem by incorporating molecular composition and spatial likelihood into boundary inference [4]. However, even improved segmentation does not eliminate a large pool of molecules that remain unassigned to any cell. These unassigned RNAs are often removed as a preprocessing step, implicitly treating them as optical noise, decoding errors, diffusion or segmentation waste.

Recent work argues that this assumption is incomplete. uRNAs can be technically confounded, but they may also contain extrasomatic RNA localization, missed cellular fragments and structured signals near cell-cell interfaces [5,6]. The analytical question is therefore not whether all uRNAs are biological. It is whether the uRNA pool contains enough reproducible tissue structure to be useful as an independent diagnostic. A useful diagnostic should not rely on the same cell assignment that it is meant to evaluate. It should instead ask whether the spatial arrangement of excluded molecules independently points back to known tissue compartments.

We tested this idea using COSTE, a topology-based framework that converts directed nearest-neighbor distance profiles into Spatial Separation Score (SSS) profiles [7,8]. COSTE can be applied to cells, cell types or single-transcript inputs because it operates on spatial coordinates and labels rather than requiring a cell-by-gene count matrix. This property makes it suitable for a segmentation-orthogonal uRNA analysis: high-quality uRNAs can be treated as independent molecular nodes, and their proximity to annotated cell groups can be summarized without reassigning them to cells.

Here we apply uRNA-only COSTE to a public 10x Genomics Atera whole-transcriptome FFPE breast cancer dataset [9]. The assay contains 18,028 genes and 624.1 million high-quality gene transcripts, providing enough uRNA density to evaluate gene-level topology. We compare uRNA-only SSS profiles with an all-transcript COSTE reference and use curated breast cancer, stromal, immune and vascular markers as positive controls. The results show that uRNAs are not a homogeneous waste stream. They preserve interpretable tissue architecture while also highlighting cases where excluded molecules diverge from cell-assigned signal.

## Results

### Whole-transcriptome uRNAs form a large spatial molecular layer

The Atera breast cancer dataset contains 170,057 detected cells and 18,028 genes in a human FFPE breast cancer section with annotated tumor, stromal, immune and rare epithelial compartments [9]. The transcript table contained 740,442,119 rows, including 738,289,569 gene rows. After applying a quality threshold of QV >= 20, 624,095,990 gene transcripts remained. Of these, 122,716,153 transcripts were assigned to the `UNASSIGNED` cell identifier, representing 19.7% of high-quality gene transcripts (Fig. 1a).

The uRNA pool was broad rather than restricted to a few genes. At least one high-quality uRNA was detected for 18,023 of 18,028 genes, corresponding to 99.97% of the panel. The median gene-level uRNA fraction across the full panel was 0.160. The distribution was not uniform: highly expressed and spatially abundant genes such as CCND1, NHERF1, XBP1 and OR4F17 contributed hundreds of thousands of high-quality uRNAs, whereas many genes had lower but still non-zero uRNA support (Fig. 1c). This scale makes the dataset suitable for testing whether uRNAs retain gene-specific tissue structure.

To avoid over-interpreting sparse genes, we selected the top 500 genes by high-quality uRNA count and added curated marker genes for interpretable breast cancer compartments. This produced 530 genes for uRNA-only COSTE analysis (Fig. 1b). For each gene, uRNA coordinates were treated as a separate spatial population, and COSTE computed SSS values between that gene population and each annotated cell group. Lower SSS indicates closer spatial proximity in the COSTE hierarchy. This design uses cell-group coordinates as reference landmarks but does not assign uRNAs to cells. Thus, all-panel summaries describe the 18,028 assayed genes, whereas COSTE and concordance results refer to the 530 analyzed genes.

The resulting uRNA-only best-cell-group calls were concentrated in expected high-density tissue compartments. Basal-like structured DCIS cells and CXCL14+ fibroblasts each were the closest cell group for 189 genes, followed by macrophages for 100 genes, endothelial cells for 16 genes and pericytes for 8 genes (Fig. 1d). This distribution indicates that uRNAs carry strong compartmental structure but also warns that high-density stromal and interface regions can dominate global uRNA topology.

### uRNA-only COSTE is moderately concordant with all-transcript topology

We next compared uRNA-only SSS profiles with all-transcript gene-to-cell SSS profiles previously generated from the same sample. Across the 530 analyzed genes, the median per-gene Spearman concordance was 0.649, with an interquartile range of 0.226 to 0.998 and a mean of 0.555 (Fig. 1e). Exact best-cell-group agreement between uRNA-only and all-transcript profiles occurred for 50.2% of genes. In total, 60.0% of genes had Spearman concordance >= 0.5, 37.7% had concordance >= 0.8 and 21.9% had concordance effectively equal to 1.0 (Fig. 1e).

These values are important because the uRNA-only analysis is intentionally deprived of cell assignment information. Perfect agreement would suggest that uRNA topology is only a redundant copy of the assigned-cell signal. No agreement would suggest that uRNAs behave as random or uninformative noise. The observed middle ground supports a more useful interpretation: uRNAs contain enough spatial structure to recover many major compartments while also emphasizing a different molecular layer.

The concordance pattern was not uniform across broad modules. Tumor and rare epithelial marker groups showed high median concordance when curated markers were considered, whereas global vascular/stromal attributions showed lower median concordance in the top-500 gene set. This pattern is consistent with the biological and technical heterogeneity of unassigned molecules. Some genes are excluded because transcripts lie near true compartment boundaries or cellular protrusions; others may reflect segmentation errors, missed fragments, local tissue geometry or abundant molecules near stromal regions. uRNA-only COSTE therefore functions as a diagnostic rather than a replacement for segmentation-aware expression analysis.

### Curated markers recover tumor, immune and vascular compartments from uRNAs alone

Curated marker genes provided positive controls with interpretable expected compartments. The 11q13 marker group, including CCND1, ELOVL2, KCNJ3 and FGF19, localized to 11q13 invasive tumor cells in the uRNA-only SSS profiles. The median Spearman concordance between uRNA-only and all-transcript profiles for this marker group was 0.999, and all four genes matched the all-transcript best cell group (Fig. 1f).

Basal-like DCIS markers also recovered their expected compartment. KRT23, DSC3, SOSTDC1, KLK5, KLK7, ITGB6 and MMP7 were closest to basal-like structured DCIS cells in uRNA-only COSTE. The group had a best-cluster match rate of 1.0 and a median Spearman concordance of 1.0. These genes are useful controls because basal-like structured DCIS is an anatomically and molecularly defined region in the dataset, and recovery from uRNAs alone argues against a purely random uRNA pool.

Immune and vascular markers showed the same principle. C3, C1QA, CSF1R, CD163 and SIGLEC1 were closest to macrophages. IGHA1, IGHM, IGHA2 and JCHAIN were closest to plasma cells. EPAS1 was closest to pericytes, whereas CDH5 and MMRN2 were closest to endothelial cells. The vascular and plasma marker groups had median concordance values of 1.0; the macrophage marker group had a median concordance of 0.698, with all five genes retaining macrophage as the best uRNA-only cell group.

Rare epithelial markers were also recovered in a biologically interpretable way. PIP, HSPB8 and CLIC6 were closest to luminal-like amorphous DCIS cells, whereas TAT localized to apocrine cells. These results are consistent with the Atera dataset annotation strategy, in which apocrine cells were identified by histology and PIP expression and tumor substructures were annotated using molecular and spatial evidence related to breast atlas and tumor-microenvironment resources [9-11]. Together, the marker analyses show that uRNAs retain gene-specific spatial proximity to expected cell groups across tumor, stromal, immune and vascular compartments.

To verify that these COSTE proximity calls corresponded to real spatial arrangements rather than only matrix-level statistics, we generated local position maps for selected examples. The selection path first moved from the 18,028-gene panel to the 530-gene uRNA-only COSTE analysis, then to curated marker controls with interpretable expected compartments, and finally to four concordant positive examples plus one discordant boundary-sensitive control (Fig. 1b,f,g). A whole-section overview first marks the five sampled windows, and the corresponding zoom panels then show the true uRNA and cell-centroid positions (Fig. 2a-f). For CCND1, C1QA, JCHAIN and CDH5, the plotted uRNAs overlapped the corresponding 11q13 tumor, macrophage, plasma-cell and endothelial landmarks in 850-um windows chosen from the true transcript and cell-coordinate data. The same examples were highlighted in the uRNA-only SSS heatmap and global concordance distribution, showing that their local spatial patterns matched their low-SSS target cell groups and high all-transcript concordance (Fig. 2g-j). ERBB2 was included as a discordant example: its uRNAs were closest to CXCL14+ fibroblasts in the uRNA-only analysis, whereas its all-transcript profile favored basal-like structured DCIS cells (Fig. 1g and Fig. 2f-j).

### Discordant uRNA profiles highlight boundary-sensitive biology

The strongest argument for using uRNAs as a diagnostic is not that every uRNA-only profile matches the all-transcript profile. It is that discordance is interpretable and actionable. Some genes with substantial uRNA support showed low or negative Spearman concordance with all-transcript SSS profiles (Supplementary Fig. 1). These discordant genes were not simply the lowest-coverage cases; several had tens to hundreds of thousands of high-quality uRNAs and uRNA fractions in the range observed for well-concordant genes.

Tumor/DCIS interface markers illustrate this point. FOXA1, CTTN, ANO1, ESR1, FADD and PGR were closest to basal-like structured DCIS cells in both analyses or retained high concordance. FGF3, however, shifted toward 11q13 invasive tumor cells in uRNA-only COSTE, whereas ERBB2 shifted toward CXCL14+ fibroblast proximity. These shifts should not be interpreted automatically as true expression reassignment. Instead, they nominate genes and regions for spatial inspection: uRNAs may accumulate near compartment boundaries, stromal interfaces, protrusions or regions where segmentation under-captures the relevant cellular morphology.

This diagnostic behavior is useful for segmentation-bias assessment. If a marker gene has a cell-assigned profile pointing to one compartment but its excluded molecules strongly point to an adjacent or different compartment, that gene becomes a candidate for reviewing masks, transcript density, local tissue morphology and cell-group annotation. Conversely, if excluded molecules reproduce the expected compartment, they provide orthogonal support that the spatial signal is robust to cell assignment. uRNA-only COSTE therefore converts a usually discarded molecular pool into a structured quality-control and discovery readout.

## Discussion

This analysis shows that high-quality unassigned RNAs in a whole-transcriptome imaging dataset preserve tissue architecture. In the Atera FFPE breast cancer section, one fifth of high-quality gene transcripts were outside cell assignments, and nearly every assayed gene had at least one uRNA. When analyzed without assigning those molecules to cells, uRNAs recovered 11q13 tumor, basal-like DCIS, macrophage, plasma-cell and vascular compartments. This finding supports a practical conclusion: uRNAs should not be treated by default as a homogeneous technical waste stream.

The result also clarifies what uRNA analysis can and cannot do. uRNA-only COSTE is not a replacement for cell segmentation, cell typing or differential expression on cell-level matrices. It does not decide whether a molecule should be assigned to a specific neighboring cell. Instead, it asks whether the molecules that were excluded from cell masks point to coherent tissue compartments. This makes it useful as an orthogonal diagnostic: it reuses annotated cell groups as landmarks but does not reuse transcript-to-cell assignments as evidence.

The moderate global concordance with all-transcript SSS profiles is central to the interpretation. If the uRNA signal were perfectly concordant, it would add little beyond existing cell-assigned analysis. If it were random, it would fail as a diagnostic. The observed mixture of high-concordance markers and discordant high-coverage genes suggests that uRNAs encode both robust compartmental signal and boundary-sensitive variation. For practical spatial transcriptomics workflows, this means uRNAs can help prioritize where segmentation masks, cell-type annotations and transcript-level spatial patterns should be inspected more carefully.

Several limitations remain. The present draft is based on one public Atera breast cancer dataset and one annotation set. It does not benchmark multiple tissues, platforms or segmentation algorithms. It also does not decompose uRNAs into technical artifacts, missed cells, diffusion-related molecules, extracellular RNAs or true extrasomatic localization. A complete Nature Methods Article should therefore add multi-dataset validation, synthetic segmentation perturbations, comparison against transcript-aware segmentation/reassignment methods and prospective biological validation of selected discordant genes.

Despite these limitations, the analysis establishes a simple and scalable route for exploiting a large molecular layer that is usually discarded. As whole-transcriptome in situ assays increase transcript density, uRNA pools will become large enough to support gene-level spatial topology. Treating these molecules as independent spatial evidence can improve segmentation-bias diagnostics and may reveal tissue-interface biology that is partially hidden by cell-body-centric analysis.

## Online Methods

### Data source

We analyzed the public 10x Genomics Atera In Situ Gene Expression preview dataset for FFPE human breast cancer [9]. The dataset was generated with a pre-commercial Atera whole-transcriptome assay containing 18,028 genes and includes cell segmentation generated by the 10x Genomics in situ multimodal segmentation solution. The 10x page reports 170,057 detected cells, 624,095,990 total high-quality decoded transcripts and a region area of 58,944,371.2 um2. Cell-group annotations were provided with the dataset and were based on graph-based clustering, marker genes, histology and published breast atlas resources [9-11].

### uRNA definition and filtering

Transcript rows were read from `transcripts.parquet`. We retained rows with `is_gene=True` and `qv >= 20`. uRNAs were defined as high-quality gene transcripts with `cell_id == UNASSIGNED`. Assigned transcripts were defined as high-quality gene transcripts with any other `cell_id`. Gene-level uRNA counts, assigned counts and uRNA fractions were computed across the full transcript table.

### Gene selection

Genes were eligible for uRNA-only COSTE if they had at least one high-quality uRNA. For the main analysis, we selected the top 500 genes by high-quality uRNA count and added curated marker genes for 11q13-amplified tumor cells, basal-like DCIS, macrophages, plasma cells, vascular structures, rare epithelial compartments and tumor/DCIS interface genes. This yielded 530 analyzed genes.

### Cell-group landmarks

Cell centroids were read from the Atera cell table and merged with the cell-group annotation file. Cells labeled as `Unassigned` in the annotation file were removed from the cell landmark set. The remaining cells were grouped by annotated cell type and used as COSTE reference populations.

### uRNA-only COSTE

For each analyzed gene, the coordinates of its high-quality uRNAs were treated as an independent spatial population. Mean nearest-neighbor distances were computed from each cell group to each other cell group, from the gene-uRNA population to each cell group, and from each cell group to the gene-uRNA population. The resulting directed distance matrix was transformed into a COSTE cophenetic representation using average-linkage hierarchical clustering. SSS values were min-max normalized to the range 0 to 1 within the cophenetic matrix, where lower values indicate closer spatial proximity. The gene row of this matrix yielded a gene-to-cell-group uRNA-only SSS profile.

### Comparison with all-transcript COSTE

uRNA-only SSS profiles were compared with a previously generated all-transcript gene-to-cell SSS table for the same dataset. For each gene, we computed Spearman correlation between the uRNA-only and all-transcript SSS vectors over common cell groups. We also compared the best cell group in each profile, defined as the cell group with the lowest SSS.

### Marker analysis

Curated marker groups were defined before analysis. The 11q13 group included CCND1, ELOVL2, KCNJ3 and FGF19. Basal-like DCIS markers included KRT23, DSC3, SOSTDC1, KLK5, KLK7, ITGB6 and MMP7. Macrophage markers included C3, C1QA, CSF1R, CD163 and SIGLEC1. Plasma-cell markers included IGHA1, IGHM, IGHA2 and JCHAIN. Rare epithelial markers included PIP, HSPB8, CLIC6 and TAT. Tumor/DCIS interface genes included FOXA1, CTTN, ANO1, ESR1, FADD, PGR, ERBB2 and FGF3. Vascular markers included EPAS1, CDH5 and MMRN2.

### Spatial example windows

Spatial example windows were extracted from the cached uRNA-coordinate parquet file generated on the A100 server. We selected CCND1, C1QA, JCHAIN and CDH5 as concordant examples and ERBB2 as a discordant example. For each gene, an 850-um square window was chosen by scoring grid bins for joint enrichment of the gene's uRNA coordinates and the target cell group. Target groups were 11q13 invasive tumor cells for CCND1, macrophages for C1QA, plasma cells for JCHAIN, endothelial cells for CDH5 and CXCL14+ fibroblasts for ERBB2. The Figure 2 overview uses a reproducible sample of 120,000 annotated cell centroids to show where the windows sit in the whole section. The zoom maps show true uRNA coordinates, true target cell centroids and sampled non-target cell centroids from each selected window. To maintain readability, uRNAs were downsampled to at most 7,000 plotted points per example and non-target cells to at most 3,500 plotted points; all target cells in the window were retained.

### Statistics and visualization

All summary statistics were computed in Python from the generated CSV files. Concordance was summarized using Spearman correlation, median, mean and interquartile range. Exact best-cell-group agreement was reported as a fraction of genes. Figures were generated with Matplotlib from the committed benchmark outputs using `make_manuscript_assets.py`.

### Software and reproducibility

The large transcript-level uRNA analysis was run with `benchmarking/atera_wta_breast_urna_segmentation_bias.py`. The manuscript-level tables and figures were generated with `manuscripts/atera_urna_nature_methods/make_manuscript_assets.py`. The latter script can be rerun from the repository root without reprocessing the full transcript parquet file.

### Use of large language models

OpenAI ChatGPT/Codex was used to help organize the manuscript draft, write summary prose and generate manuscript-asset code from author-provided analysis results. All numerical results are derived from the cited data files and scripts, and the authors remain responsible for verification, interpretation and final text.

### Data availability

The Atera FFPE human breast cancer dataset is publicly available from 10x Genomics. Derived summary files used in this draft are included in the repository under `benchmarking/atera_wta_breast_urna_segmentation_bias_results/` and `manuscripts/atera_urna_nature_methods/tables/`. Large raw transcript files are not redistributed in this repository.

### Code availability

The analysis scripts are available in this repository at `benchmarking/atera_wta_breast_urna_segmentation_bias.py` and `manuscripts/atera_urna_nature_methods/make_manuscript_assets.py`. A DOI-backed software release should be created before submission.

## Figure Legends

### Graphical Abstract

High-quality uRNAs are extracted from whole-transcriptome imaging data and analyzed as independent molecular nodes. uRNA-only COSTE produces gene-to-cell-group SSS profiles that can be compared with all-transcript profiles to diagnose segmentation-orthogonal spatial structure.

### Figure 1. Integrated uRNA evidence chain before spatial validation

a, Transcript counts through gene, quality and uRNA filtering. b, Analysis-set definition from the 18,028-gene panel to the 530 genes analyzed by uRNA-only COSTE and concordance. c, Distribution of gene-level uRNA fractions across all panel genes, with the 530 analyzed genes overlaid. d, uRNA-only best COSTE cell groups for the 530 analyzed genes. e, Compact summary of global concordance between uRNA-only and all-transcript SSS profiles across the 530 analyzed genes. f, Curated marker-control logic showing expected compartments recovered from uRNAs alone and the heterogeneous tumor/DCIS interface control. g, Spatial-example bridge showing the four concordant positive examples and the discordant ERBB2 control carried into Figure 2. Detailed five-gene SSS heatmaps, concordance placement and per-example statistics are shown only in Figure 2.

### Figure 2. Real spatial examples support uRNA-only COSTE calls

a, Whole-section overview of sampled cell centroids with colored boxes marking the five 850-um windows. b-f, True spatial windows for CCND1, C1QA, JCHAIN, CDH5 and ERBB2. Black points are high-quality uRNAs for the plotted gene, open colored circles are the target cell group and pale gray points are other cells in the same window. g, uRNA-only COSTE proximity for the five examples and representative cell groups, plotted as 1 - SSS so higher values indicate closer proximity. Black boxes mark the plotted target cell group. h, Global distribution of Spearman concordance between uRNA-only and all-transcript SSS profiles, with selected examples highlighted. i, Marker group validation summary. j, Per-example statistics, including uRNA-only SSS, concordance, best-cell-group match and the uRNA-only versus all-transcript best cell group.

## References

1. Tian, L., Chen, F. & Macosko, E. Z. The expanding vistas of spatial transcriptomics. Nat. Biotechnol. 41, 773-782 (2023). https://doi.org/10.1038/s41587-022-01448-2
2. Dries, R., Zhu, Q., Dong, R. et al. Giotto: a toolbox for integrative analysis and visualization of spatial expression data. Genome Biol. 22, 78 (2021). https://doi.org/10.1186/s13059-021-02286-2
3. Palla, G., Spitzer, H., Klein, M. et al. Squidpy: a scalable framework for spatial omics analysis. Nat. Methods 19, 171-178 (2022). https://doi.org/10.1038/s41592-021-01358-2
4. Petukhov, V., Xu, R. J., Soldatov, R. A. et al. Cell segmentation in imaging-based spatial transcriptomics. Nat. Biotechnol. 40, 345-354 (2022). https://doi.org/10.1038/s41587-021-01044-w
5. Marco Salas, S., Dammann, M., Rubens, R. K. et al. Exploration of RNA outside segmented cells in spatial transcriptomics reveals extrasomatic RNA organization. Preprint at bioRxiv (2025). https://doi.org/10.64898/2025.12.07.692889
6. Yuan, L., Zheng, Y. P. et al. Reconstructing biologically coherent cellular profiles from imaging-based spatial transcriptomics. Preprint at bioRxiv (2026). https://doi.org/10.64898/2026.03.08.710395
7. Long, M., Hu, T., Sountoulidis, A., Samakovlis, C. & Nilsson, M. Cophenetic Spatial Topology Embedding reveals multiscale tissue architecture in spatial omics. Preprint at bioRxiv (2026). https://doi.org/10.64898/2026.05.26.727847
8. Hu, T., Long, M. & Nilsson, M. Cell-GPS (COSTE Windows Application): Cophenetic Spatial Topology Embedding for multiscale tissue architecture analysis in spatial omics. Zenodo (2026). https://doi.org/10.5281/zenodo.19482685
9. 10x Genomics. Preview Data: Atera In Situ Gene Expression, FFPE Human Breast Cancer. https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-breast-cancer
10. Kumar, T., Nee, K., Wei, R. et al. A spatially resolved single-cell genomic atlas of the adult human breast. Nature 620, 181-191 (2023). https://doi.org/10.1038/s41586-023-06252-9
11. Janesick, A., Shelansky, R., Gottscho, A. D. et al. High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis. Nat. Commun. 14, 8353 (2023). https://doi.org/10.1038/s41467-023-43458-x
