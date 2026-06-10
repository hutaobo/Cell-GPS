"""Apply the circle-legend rebuild to the two MAIN-figure source SVGs that feed
Figure 1 (transcript-point HistoSeg map) and Figure 2 (point-vs-cell spatial
comparison), then re-render their PNG/PDF at the same DPI that the original
matplotlib output used.

Reads from figures/main_backup_20260608_circle/ (pristine source backups) and
writes to figures/main/ so a future `make_round1_revision_package.py` run picks
up the new colorbars.
"""

from __future__ import annotations

import os

import fitz
from PIL import Image

import relayout_supplementary_figures as RS

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "figures", "main_backup_20260608_circle"))
DST = os.path.normpath(os.path.join(HERE, "..", "figures", "main"))

# stem -> circle-legend params. height_k = legend column height relative to the
# pristine matplotlib colorbar height.
#   Figure 1b (transcript_point_histoseg_domains): 0.70x
#   Figure 2a (point_vs_cell_spatial_comparison):  0.50x  (user request 2026-06-08)
TARGETS = {
    "transcript_point_histoseg_domains": {"height_k": 0.70},
    "point_vs_cell_spatial_comparison":  {"height_k": 0.50},
}


def _render(svg_text, stem_out, dpi=300, make_pdf=True):
    open(stem_out + ".svg", "w", encoding="utf-8", newline="\n").write(svg_text)
    doc = fitz.open("svg", svg_text.encode("utf-8"))
    if make_pdf:
        open(stem_out + ".pdf", "wb").write(doc.convert_to_pdf())
    pix = doc[0].get_pixmap(matrix=fitz.Matrix(dpi / 72, dpi / 72), alpha=False)
    tmp = stem_out + ".tmp.png"
    pix.save(tmp)
    with Image.open(tmp) as im:
        im.convert("RGB").save(stem_out + ".png", dpi=(dpi, dpi))
    os.remove(tmp)


def main():
    for stem, params in TARGETS.items():
        src_svg = os.path.join(SRC, stem + ".svg")
        svg = open(src_svg, encoding="utf-8").read().replace("\r\n", "\n").replace("\r", "\n")
        # only the colorbar rebuild — main-figure layout is otherwise correct.
        svg = RS._redo_colorbar_as_circle_legend(
            svg, height_k=params.get("height_k", 0.70),
            circle_r=params.get("circle_r", 3.8),
            label_font_pt=params.get("label_font_pt", 8.0),
            label_gap=params.get("label_gap", 3.5))
        stem_out = os.path.join(DST, stem)
        _render(svg, stem_out, dpi=300, make_pdf=True)
        print("wrote", stem_out + ".{svg,png,pdf}")


if __name__ == "__main__":
    main()
