"""Regenerate Figure 2b (point_vs_cell_overlap_heatmap) from the local overlap
matrix CSV, with the in-cell annotation numbers enlarged by 2 pt (7 -> 9).

This is a faithful replica of plot_heatmap() in
atera_breast_point_vs_cell_histoseg_comparison.py (matplotlib defaults, viridis,
figsize 10x8, constrained_layout); the ONLY change is the annotation fontsize.
The matrix is read from figures/source_data/point_vs_cell_overlap_point_fraction.csv
(= the `row_fraction` DataFrame the original used).
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
CSV = os.path.join(ROOT, "figures", "source_data", "point_vs_cell_overlap_point_fraction.csv")
OUT = os.path.join(ROOT, "figures", "main", "point_vs_cell_overlap_heatmap")

ANNOT_FONTSIZE = 9   # was 7 in the source; "+2"


def main():
    row_fraction = pd.read_csv(CSV, index_col=0)

    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    image = ax.imshow(row_fraction.to_numpy(dtype=float), vmin=0,
                      vmax=max(0.5, float(row_fraction.max().max())), cmap="viridis")
    ax.set_xticks(np.arange(row_fraction.shape[1]))
    ax.set_xticklabels([c.replace("COSTE-HistoSeg-", "C") for c in row_fraction.columns],
                       rotation=45, ha="right")
    ax.set_yticks(np.arange(row_fraction.shape[0]))
    ax.set_yticklabels([i.replace("TranscriptPoint-HistoSeg-", "P") for i in row_fraction.index])
    ax.set_xlabel("cell-based COSTE-HistoSeg domain")
    ax.set_ylabel("transcript-point HistoSeg domain")
    ax.set_title("Transcript-point to cell-domain overlap (row fraction)")
    for i in range(row_fraction.shape[0]):
        for j in range(row_fraction.shape[1]):
            value = float(row_fraction.iloc[i, j])
            if value >= 0.05:
                ax.text(j, i, f"{value:.2f}", ha="center", va="center",
                        color="white" if value > 0.25 else "black",
                        fontsize=ANNOT_FONTSIZE)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("fraction of point-domain transcripts")
    fig.savefig(OUT + ".png", dpi=300)
    fig.savefig(OUT + ".svg")
    plt.close(fig)
    print("wrote", OUT + ".{png,svg}", f"(annotation fontsize={ANNOT_FONTSIZE})")


if __name__ == "__main__":
    main()
