"""Build Supplementary Figure 6: transcript-point HistoSeg n-domain sensitivity.

a) Point-vs-cell-12 agreement (ARI, NMI, point-to-cell purity, weighted Jaccard)
   across k = 8, 10, 12, 14, 16 transcript-point domains.
b) Fraction of transcript points falling in tumor-marker vs immune-stromal-marker
   domains across k (the major molecular partition is invariant to k).

Renders pdf/png/svg/tiff into figures/supplementary, and installs the formal
Supplementary_Figure_6.* into the bundle + all_figure_pdfs + SUBMISSION folder.
"""
from __future__ import annotations
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
SUPP = os.path.join(ROOT, "figures", "supplementary")
BUNDLE = os.path.join(ROOT, "review_package_round1", "submission_RECOMMENDED_local_bundle", "figures")
ALLPDF = os.path.join(ROOT, "review_package_round1", "all_figure_pdfs_RECOMMENDED")
SUBFOLDER = os.path.join(ROOT, "review_package_round1", "submission_RECOMMENDED_local_bundle",
                         "SUBMISSION_manuscript_and_figure_pdfs")
df = pd.read_csv(os.path.join(ROOT, "analysis_sensitivity", "sens_fig_data.csv"))

plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
INK = "#1f2933"; TUMOR = "#C2185B"; IMMUNE = "#008C8C"
k = df["k"].to_numpy()

fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.0, 2.9))
fig.subplots_adjust(left=0.08, right=0.985, bottom=0.16, top=0.88, wspace=0.32)

# panel a: agreement metrics vs k
# Okabe-Ito colour-blind-safe quartet, distinct from panel b's tumor/immune
# colours and with per-metric marker shapes for extra deuteranopia robustness.
metrics = [("ARI", "#0072B2", "o"), ("NMI", "#E69F00", "s"), ("purity", "#009E73", "^"), ("jaccard", "#CC79A7", "D")]
for name, col, mk in metrics:
    axa.plot(k, df[name], "-", marker=mk, color=col, lw=1.3, ms=4, label=name)
axa.set_xlabel("number of transcript-point domains (k)")
axa.set_ylabel("agreement vs cell-based 12")
axa.set_xticks(k); axa.set_ylim(0, 0.6)
axa.legend(frameon=False, fontsize=6.6, ncol=2, loc="upper right", handlelength=1.2, columnspacing=1.0)
for s in ("top", "right"):
    axa.spines[s].set_visible(False)
axa.tick_params(colors=INK, labelsize=7)
axa.set_title("a", loc="left", fontsize=8, fontweight="bold", x=-0.13, y=1.0)

# panel b: tumor vs immune point fraction (stacked, sums to 1)
axb.bar(k, df["tumor_frac"], width=1.2, color=TUMOR, label="tumor-marker domains")
axb.bar(k, df["immune_frac"], width=1.2, bottom=df["tumor_frac"], color=IMMUNE, label="immune-stromal-marker domains")
for xi, t in zip(k, df["tumor_frac"]):
    axb.text(xi, t / 2, f"{t:.2f}", ha="center", va="center", color="white", fontsize=6.4, fontweight="bold")
for xi, t, im in zip(k, df["tumor_frac"], df["immune_frac"]):
    axb.text(xi, t + im / 2, f"{im:.2f}", ha="center", va="center", color="white", fontsize=6.4, fontweight="bold")
axb.set_xlabel("number of transcript-point domains (k)")
axb.set_ylabel("fraction of transcript points")
axb.set_xticks(k); axb.set_ylim(0, 1.0)
axb.legend(frameon=False, fontsize=6.6, loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=2, handlelength=1.0)
for s in ("top", "right"):
    axb.spines[s].set_visible(False)
axb.tick_params(colors=INK, labelsize=7)
axb.set_title("b", loc="left", fontsize=8, fontweight="bold", x=-0.13, y=1.0)

stem = os.path.join(SUPP, "n_domain_sensitivity")
fig.savefig(stem + ".pdf", bbox_inches="tight")
fig.savefig(stem + ".svg", bbox_inches="tight")
fig.savefig(stem + ".png", dpi=600, bbox_inches="tight")
plt.close(fig)
with Image.open(stem + ".png") as im:
    im.convert("RGB").save(stem + ".tiff", dpi=(600, 600), compression="tiff_lzw")

# install as formal Supplementary_Figure_6.*
import shutil
for ext in ("pdf", "png", "svg", "tiff"):
    shutil.copy2(stem + "." + ext, os.path.join(BUNDLE, "Supplementary_Figure_6." + ext))
shutil.copy2(stem + ".pdf", os.path.join(ALLPDF, "Supplementary_Figure_6.pdf"))
shutil.copy2(stem + ".pdf", os.path.join(SUBFOLDER, "Supplementary_Figure_6.pdf"))
print("wrote Supplementary_Figure_6 (pdf/png/svg/tiff) + installed to bundle/all_pdfs/SUBMISSION")
