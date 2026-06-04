# Cover letter draft

Dear Editors,

We submit the manuscript "Unassigned RNAs report segmentation-orthogonal tissue architecture in whole-transcriptome imaging data" for consideration as a Nature Methods Article.

Cell segmentation is a central analytical step in imaging-based spatial transcriptomics, yet a substantial fraction of high-quality molecules remains outside cell masks and is often removed before downstream analysis. We present a simple segmentation-orthogonal diagnostic in which unassigned RNAs are treated as independent spatial nodes and analyzed with Cophenetic Spatial Topology Embedding (COSTE). Applied to the public 10x Genomics Atera whole-transcriptome FFPE breast cancer dataset, this approach evaluates whether molecules excluded from cell-level count matrices preserve tissue architecture rather than behaving as random technical noise.

The study analyzes 624.1 million high-quality gene transcripts, including 122.7 million high-quality unassigned RNAs. uRNAs are detected for 18,023 of 18,028 assayed genes. Across 530 high-coverage or curated marker genes, uRNA-only COSTE profiles show moderate concordance with all-transcript COSTE profiles and recover expected tumor, stromal, immune and vascular compartments. Canonical markers for 11q13-amplified tumor cells, basal-like DCIS, macrophages, plasma cells and vascular structures retain cell-type-specific spatial proximity from uRNA coordinates alone.

The work is appropriate for Nature Methods because it addresses a practical problem in spatial omics analysis: how to evaluate segmentation bias without treating cell assignment as the only ground truth. It also provides a reusable framework for whole-transcriptome imaging datasets in which unassigned transcript pools are large enough to support gene-level spatial topology analysis.

Related work from the authors includes the COSTE/Cell-GPS method and software. This manuscript extends that framework to a segmentation-orthogonal uRNA diagnostic and should be considered alongside the disclosed related manuscript or preprint.

Sincerely,

[Corresponding author name]
