# di-sim / COSTE Atera benchmark: correspondence memo

## Purpose

This memo packages the Atera WTA breast benchmark into a concise form suitable for
emailing the corresponding authors of the PNAS di-sim paper:

Rohe, K., Qin, T. & Yu, B. (2016). *Co-clustering directed graphs to discover
asymmetries and directional communities*. PNAS 113(45), 12679-12684.
DOI: https://doi.org/10.1073/pnas.1525793113

The message should be respectful and clear:

- di-sim is the prior conceptual method for directed asymmetric co-clustering.
- COSTE was developed independently for spatial omics, but it uses a closely related
  directed Searcher-to-Findee mathematical object.
- The benchmark is not framed as "COSTE beats di-sim." It is framed as "di-sim is
  an elegant spectral way to expose asymmetry; COSTE is a hierarchical/cophenetic
  way to preserve the full directed profile geometry."
- The useful result is that di-sim can screen the whole 18,028-gene Atera panel,
  after which COSTE gives a detailed cophenetic hierarchy for the 1,024 strongest
  directional genes.

## Suggested recipients

The PNAS article lists Karl Rohe and Bin Yu as corresponding authors. Current public
pages also show active contact information.

- Karl Rohe: karlrohe@stat.wisc.edu
- Bin Yu: binyu@berkeley.edu or binyu@stat.berkeley.edu

## One-paragraph technical summary

We built a directed gene-to-gene matrix from the public 10x Genomics Atera WTA FFPE
breast cancer Xenium dataset. Each entry summarizes how the spatial distribution of
one gene relates directionally to another gene's spatial distribution. On this same
asymmetric matrix, di-sim computes row/source and column/target spectral embeddings
from a regularized directed graph Laplacian, whereas COSTE clusters the directed
row and column profiles and converts those hierarchies into normalized cophenetic
distance matrices. In this benchmark, full-panel di-sim was practical for all
18,028 genes, and a di-sim importance score selected 1,024 genes for COSTE. The
selected COSTE result recovered interpretable breast cancer spatial programs and
supported downstream HistoSeg-style spatial domains.

## Benchmark evidence to mention

Data:

- Public 10x Genomics Atera WTA FFPE human breast cancer dataset.
- 170,057 cells.
- 18,028 predesigned gene targets.
- 624,095,990 high-quality decoded gene transcripts.

Whole-matrix runs:

| run | matrix | matrix build sec | COSTE sec | di-sim sec | peak build GB | peak COSTE GB | peak di-sim GB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| celltype_full | 20 x 20 | 3.92 | 0.15 | 0.33 | 0.23 | 0.19 | 0.19 |
| gene_top1024_grid96 | 1024 x 1024 | 31.72 | 0.86 | 3.40 | 5.30 | 0.90 | 0.92 |
| gene_top2048_grid96 | 2048 x 2048 | 41.82 | 22.25 | 4.08 | 5.82 | 1.04 | 1.01 |
| gene_top4096_grid96 | 4096 x 4096 | 44.39 | 189.03 | 11.48 | 6.19 | 1.63 | 1.70 |

Full di-sim then selected COSTE:

| step | items | seconds | peak GB |
| --- | ---: | ---: | ---: |
| matrix build | 18,028 genes | 52.60 | 6.98 |
| full di-sim | 18,028 genes | 60.20 | 14.02 |
| selected COSTE cophenetic | 1,024 genes | 0.86 | 2.13 |

Cell-free transcript-point version:

| step | items | seconds | peak GB |
| --- | ---: | ---: | ---: |
| transcript metadata scan | 18,028 genes | 97.87 | 0.98 |
| transcript grid aggregation | 18,028 genes | 62.27 | 1.34 |
| cell-free directed matrix | 18,028 genes | 4.03 | 3.63 |
| full di-sim | 18,028 genes | 28.28 | 7.27 |
| selected COSTE cophenetic | 1,024 genes | 0.82 | 3.69 |
| transcript-point HistoSeg assignment | 60,906,532 points | 356.29 | 0.94 |

Selected-gene evidence:

- The 1,024 di-sim-selected genes include interpretable breast cancer and
  microenvironment markers such as SLC30A8, CLIC6, HSPB8, GRIA2, NIBAN1, PIP,
  SERPINA6, MYBPC1, TAT, KCNQ3, MSMB, BMPER, and SERPINA1.
- COSTE over the selected genes produced 16 gene modules.
- HistoSeg-style domains from those modules recovered tumor-rich, immune/stroma,
  luminal/DCIS, and apocrine/mixed microenvironment regions.
- The transcript-point version assigns 60.9M selected-gene transcript points
  directly, without a spatial grid and without relying on cells as the analysis unit.

## Interpretation in plain language

di-sim is useful as a fast first pass over the whole directed gene network. It finds
which genes behave differently as "sources" versus "targets" in the asymmetric spatial
matrix. COSTE is useful after that because it keeps the full directed profile shape and
turns it into a hierarchical cophenetic geometry. So the two methods are not redundant:
di-sim is a spectral asymmetry detector, and COSTE is a topology-preserving hierarchical
readout.

The most natural combined workflow is:

1. Build one directed spatial matrix.
2. Run full-panel di-sim on all 18,028 genes.
3. Select genes with strong directional/asymmetric behavior.
4. Run COSTE on the selected subset.
5. Use COSTE modules for spatial interpretation and HistoSeg-style domains.

This is exactly why the comparison is scientifically friendly rather than competitive:
the PNAS di-sim idea gives an elegant way to expose directionality at scale; COSTE
uses that same directed-matrix view but emphasizes hierarchical spatial topology.

## Caveats to be transparent about

- The spatial gene-to-gene matrix is not the original binary directed graph used in
  the PNAS examples. It is a weighted directed matrix derived from spatial expression.
- The current benchmark uses one public Atera breast dataset.
- COSTE's exact cophenetic step scales less favorably than truncated di-sim SVD, so
  full-panel COSTE is not the right comparison at 18,028 genes unless approximate or
  blockwise strategies are added.
- The HistoSeg interpretation is molecular/pathology-guided, not a formal diagnostic
  pathology label.
- The H&E overlays use a sampled visual preview for readability, while the statistics
  use all 60,906,532 selected-gene transcript points.

## Suggested attachments

Attach these in the email:

1. `../../benchmarking/atera_breast_disim_coste_benchmark_summary.md`
2. The three Nature-style story figures:
   - `figures/nature_pathology_story_tumor_rich_steps.png`
   - `figures/nature_pathology_story_immune_stroma_steps.png`
   - `figures/nature_pathology_story_mixed_luminal_apocrine_steps.png`
3. Optionally, the source code scripts:
   - `atera_breast_disim_coste_benchmark.py`
   - `atera_breast_transcript_point_histoseg.py`
   - `atera_breast_point_vs_cell_histoseg_comparison.py`
   - `atera_breast_nature_pathology_figures.py`

## Draft email

Subject: di-sim and COSTE on spatial omics directed matrices

Dear Professor Rohe and Professor Yu,

I recently read your PNAS paper on di-sim, "Co-clustering directed graphs to
discover asymmetries and directional communities." I found the method genuinely
beautiful. The idea of preserving the row/source and column/target roles in a
directed matrix, instead of symmetrizing away the interesting part, is exactly
the kind of idea that feels obvious only after someone very clever has made it
visible.

I am writing because I have been developing a spatial omics method called COSTE
(Cophenetic Spatial Topology Embedding). I developed it independently and only
recently read your di-sim paper in detail. After reading it, I realized that the
mathematical spirit is closely related: both methods start from an asymmetric
directed object with different source and target roles. di-sim uses a
regularized directed graph Laplacian and SVD to produce left/right embeddings
and co-clusters; COSTE instead clusters the directed row and column profiles and
uses the resulting cophenetic distances as a spatial topology readout.

To make the relationship concrete, I ran a small benchmark on a public 10x
Genomics Atera WTA FFPE breast cancer Xenium dataset. The dataset contains
170,057 cells, 18,028 genes, and about 624M high-quality decoded gene
transcripts. I built a directed gene-to-gene spatial matrix and applied a
di-sim-like implementation and COSTE to the same matrix.

The most useful workflow was:

1. run full-panel di-sim on all 18,028 genes;
2. select 1,024 genes with strong directional/asymmetric behavior;
3. run COSTE on those selected genes;
4. use COSTE modules for HistoSeg-style spatial interpretation.

In this run, full di-sim over all 18,028 genes took about 60 seconds on the
cell-derived matrix, and the selected 1,024-gene COSTE cophenetic step took
about 0.86 seconds. In a cell-free transcript-point version, full di-sim took
about 28 seconds and selected COSTE took about 0.82 seconds. The selected genes
included interpretable breast cancer and microenvironment markers such as
SLC30A8, CLIC6, HSPB8, PIP, SERPINA6, KCNQ3, MSMB, BMPER, and SERPINA1. The
downstream COSTE/HistoSeg maps recovered tumor-rich, immune/stroma-interface,
luminal/DCIS, and apocrine/mixed spatial regions.

My current interpretation is that di-sim and COSTE are not competitors so much
as two complementary views of the same directed-matrix idea. di-sim is a fast
and elegant spectral asymmetry detector; COSTE preserves more of the full
directed profile geometry as a hierarchical/cophenetic spatial readout.

I would be very grateful for your thoughts on whether this is a reasonable way
to connect the two methods, especially whether you see any theoretical issue
with adapting the di-sim view to weighted directed spatial matrices of genes.
If useful, I would be happy to share the benchmark code, summary tables, and
figures. I would also be glad to briefly discuss by email or Zoom.

With sincere appreciation for your work,

[Your name]

## Source links checked

- PNAS / PMC article page: https://pmc.ncbi.nlm.nih.gov/articles/PMC5111689/
- DOI: https://doi.org/10.1073/pnas.1525793113
- Karl Rohe current page: https://stat.wisc.edu/staff/rohe-karl/
- Bin Yu current page: https://binyu.stat.berkeley.edu/
- 10x Genomics Atera WTA breast dataset page: https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-breast-cancer
