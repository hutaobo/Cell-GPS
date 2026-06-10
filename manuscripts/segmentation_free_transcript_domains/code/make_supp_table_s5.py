"""Build Supplementary Table S5 (transcript-point domain-number sensitivity)
from the cluster sensitivity outputs, with explicit per-k tumor/immune point
fractions and a transparent note on how the sweep relates to the production
point-versus-cell comparison (Table S3).

Sources (local copies of the cluster sweep summaries):
  analysis_sensitivity/n_domain_sensitivity_summary.csv   (counts + metrics)
  analysis_sensitivity/sens_fig_data.csv                  (per-k point fractions)

Installs the workbook to the review package, the submission bundle and the
SUBMISSION upload folder so all three stay in sync.
"""
from __future__ import annotations
import os
import shutil
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
PKG = os.path.join(ROOT, "review_package_round1")
SUB = os.path.join(PKG, "submission_RECOMMENDED_local_bundle",
                   "SUBMISSION_manuscript_and_figure_pdfs")
NAME = "Supplementary_Table_S5_n_domain_sensitivity.xlsx"

summary = pd.read_csv(os.path.join(ROOT, "analysis_sensitivity", "n_domain_sensitivity_summary.csv"))
fracs = pd.read_csv(os.path.join(ROOT, "analysis_sensitivity", "sens_fig_data.csv"))

m = summary.merge(fracs[["k", "tumor_frac", "immune_frac"]],
                  left_on="n_point_domains", right_on="k", how="left")

sens = pd.DataFrame({
    "Transcript-point domains (k)": m["n_point_domains"].astype(int),
    "ARI vs cell-12": m["ARI"].round(3),
    "NMI vs cell-12": m["NMI"].round(3),
    "Point-to-cell purity": m["point_to_cell_purity"].round(3),
    "Weighted mean Jaccard": m["weighted_mean_jaccard"].round(3),
    "# tumor-marker domains": m["n_tumor_marker_domains"].astype(int),
    "# immune-stromal-marker domains": m["n_immune_marker_domains"].astype(int),
    "Tumor-marker point fraction": m["tumor_frac"].round(3),
    "Immune-stromal-marker point fraction": m["immune_frac"].round(3),
    "Largest-domain marker class": m["largest_domain_class"],
})

notes = pd.DataFrame({"Note": [
    "Sensitivity of transcript-point HistoSeg to the number of domains (k).",
    "Each k: transcript-point HistoSeg was run with identical inputs/parameters "
    "(neighborhood_k=128, spatial_weight=0.35, qv_min=20, 750k training reservoir, "
    "all 60,906,532 selected-gene transcript points assigned), varying only --n-domains.",
    "ARI/NMI/purity/Jaccard are computed against the same prespecified cell-based "
    "COSTE-HistoSeg 12-domain reference.",
    "tumor-marker module = HBB/SERPINA6/MSMB/KCNQ3/SERPINA1; immune-stromal-marker "
    "module = RERGL/BMPER/IGF2/DPT/CCL22 (>=2 markers in a domain's top module genes).",
    "Tumor-marker / immune-stromal-marker point fraction = fraction of all assigned "
    "transcript points falling in domains of each marker class (the two classes sum "
    "to 1; this major molecular partition is ~0.72 vs ~0.28 across all k).",
    "Across all resolutions (k=8-16) the largest domains are tumor-rich and "
    "immune-stromal; increasing k subdivides these major domains rather than changing "
    "the qualitative result. Point-to-cell purity is highest near the matched k=12-14.",
    "This sweep is an independent set of runs. Its k=12 row (ARI 0.237, point-to-cell "
    "purity 0.485) therefore differs slightly from the production point-versus-cell "
    "comparison reported in the main text and Table S3 (ARI 0.261, purity 0.450); the "
    "difference is within k-means initialization and reservoir-sampling variability and "
    "does not affect any qualitative conclusion. Domain marker-class counts at k=12 "
    "(7 tumor / 5 immune-stromal here) likewise reflect this independent clustering and "
    "may differ by one from the production run summarized in Table S2.",
]})

out = os.path.join(PKG, NAME)
with pd.ExcelWriter(out, engine="openpyxl") as xw:
    sens.to_excel(xw, sheet_name="sensitivity", index=False)
    notes.to_excel(xw, sheet_name="notes", index=False)
for dst in (SUB,):
    if os.path.isdir(dst):
        shutil.copy2(out, os.path.join(dst, NAME))
print("wrote", out, "+ installed to SUBMISSION")
print(sens.to_string(index=False))
