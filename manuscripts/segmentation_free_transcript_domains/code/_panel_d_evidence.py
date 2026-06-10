"""Synthesise the new Supplementary Figure 3/4/5 panel d ("evidence card")
as a self-contained SVG <g> block in local coordinates (0..W, 0..H), to be
translated into the d slot of each plate by relayout_supplementary_figures.py.

Pieces, top to bottom:
  1) Header banner — story title + claim, with a story-colored accent stripe.
  2) Module genes — label + 5 colored chips (the point-side marker module).
  3) Top P->C matches — label + up to 6 rows: P->C chip (colored by pathology
     label relation) + matched dominant celltype short label + colored dot.
  4) KPI tiles — ARI / NMI / point-to-cell purity, three rounded rectangles.

Data sources (relative to the manuscript root):
  figures/source_data/nature_pathology_source_domain_statistics.csv
  figures/source_data/point_vs_cell_global_metrics.json

The output is pure SVG with no defs/clipPaths, so it composes cleanly with the
host plate SVG (no id-collision risk).
"""

from __future__ import annotations

import csv
import json
import os
from html import escape

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "figures", "source_data")


STORIES = {
    "tumor_rich": {
        "title": "Tumor-rich molecular signal",
        "claim": "Transcript domains in this module recover invasive "
                 "tumor regions.",
        "color": "#C2185B",
        "accent": "#7A003C",
        "tint":   "#FBE4EC",
    },
    "immune_stroma": {
        "title": "Immune-stromal interface",
        "claim": "These point domains identify immune-stromal "
                 "interface regions.",
        "color": "#008C8C",
        "accent": "#004F57",
        "tint":   "#E0F2F1",
    },
    "mixed_luminal_apocrine": {
        "title": "Mixed luminal/apocrine context",
        "claim": "Immune-stromal-like transcript domains appear inside "
                 "luminal/DCIS and apocrine cell-domain contexts.",
        "color": "#E68613",
        "accent": "#8F4A00",
        "tint":   "#FFF3E0",
    },
}


RELATION_STYLE = {
    "same":      ("#2E7D32", "#C8E6C9", "same"),       # green
    "related":   ("#F57F17", "#FFE0B2", "related"),    # amber
    "different": ("#616161", "#E0E0E0", "different"),  # gray
}


CELLTYPE_SHORT = {
    "11q13 Invasive Tumor Cells":         "11q13 tumor",
    "CAFs, DCIS Associated":               "CAFs (DCIS)",
    "Apocrine Cells":                      "apocrine",
    "Luminal-like Amorphous DCIS Cells":   "luminal DCIS",
    "T Lymphocytes":                       "T cells",
    "Macrophages":                         "macrophage",
    "Plasma Cells":                        "plasma",
    "Endothelial Cells":                   "endothelial",
    "Basal-like Structured DCIS Cells":    "basal DCIS",
    "CXCL14+ Fibroblasts":                 "CXCL14+ fib.",
}


def _load_matches(story_key):
    path = os.path.join(DATA, "nature_pathology_source_domain_statistics.csv")
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["story"] == story_key:
                rows.append(r)
    rows.sort(key=lambda r: -float(r["point_to_cell_fraction"]))
    return rows


def _load_metrics():
    with open(os.path.join(DATA, "point_vs_cell_global_metrics.json"), encoding="utf-8") as f:
        return json.load(f)


def _short(s):
    return CELLTYPE_SHORT.get(s, s if len(s) <= 14 else s[:13] + ".")


def _module_genes(matches):
    """The point-side top genes are constant within a story (same module)."""
    return matches[0]["point_top_genes"].split("/")


def _wrap(text, max_chars):
    out, line = [], ""
    for w in text.split():
        if line and len(line) + 1 + len(w) > max_chars:
            out.append(line); line = w
        else:
            line = w if not line else line + " " + w
    if line:
        out.append(line)
    return out


# ---------------------------------------------------------------------------
# Drawing primitives
# ---------------------------------------------------------------------------

def _rect(x, y, w, h, fill, rx=0, stroke=None, sw=0.0):
    s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    rr = f' rx="{rx}" ry="{rx}"' if rx else ""
    return f'<rect x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" fill="{fill}"{rr}{s}/>'


def _text(x, y, s, *, size=6.5, weight="normal", color="#222", anchor="start",
          family="Arial, Helvetica, sans-serif"):
    w = f' font-weight="{weight}"' if weight != "normal" else ""
    return (f'<text x="{x:.3f}" y="{y:.3f}" font-family="{family}" '
            f'font-size="{size:.2f}" fill="{color}" text-anchor="{anchor}"{w}'
            f'>{escape(s)}</text>')


def _line(x1, y1, x2, y2, color="#CFCFCF", w=0.5):
    return f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="{color}" stroke-width="{w}"/>'


def _circle(cx, cy, r, fill, stroke=None, sw=0.0):
    s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    return f'<circle cx="{cx:.3f}" cy="{cy:.3f}" r="{r:.3f}" fill="{fill}"{s}/>'


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

def make_evidence_card(story_key, box_w, box_h):
    """Return an SVG `<g>` block for panel d, in local coords (0,0)..(box_w,box_h).

    Sections auto-fit the available box_h: header + module-gene chips + match
    list (up to 5 rows) + KPI tiles. Any leftover vertical space is distributed
    proportionally so the card always fills the slot without gaping at the
    bottom or running off the top.
    """
    story = STORIES[story_key]
    matches = _load_matches(story_key)
    metrics = _load_metrics()
    genes = _module_genes(matches)

    INK = "#1f2933"
    MUTED = "#6b7886"
    pad = 5.0

    # ----- design natural heights (will be expanded to fill box) -----------
    claim_lines = _wrap(story["claim"], 36)[:3]
    natural_header = 14.4 + len(claim_lines) * 7.0 + 4.0
    chip_h_nat = 9.5
    gene_chips_h_nat = chip_h_nat + 11  # one row of chips + label
    n_show = min(len(matches), 5)
    match_row_h_nat = 13.0
    matches_h_nat = 11 + n_show * match_row_h_nat + 3
    kpi_label_h = 8.0
    kpi_h_nat = 30.0
    kpi_section_h_nat = kpi_label_h + kpi_h_nat
    sep_h = 6.0
    natural_total = (natural_header + sep_h + gene_chips_h_nat + sep_h
                     + matches_h_nat + sep_h + kpi_section_h_nat + 2 * pad)
    slack = max(0.0, box_h - natural_total)
    # Distribute slack: header gets 5%, gene-chips 5%, matches 70%, kpi 20%.
    h_header = natural_header + slack * 0.05
    h_gene = gene_chips_h_nat + slack * 0.05
    extra_match_row = (slack * 0.70) / max(n_show, 1)
    match_row_h = match_row_h_nat + extra_match_row
    h_matches = 11 + n_show * match_row_h + 3
    h_kpi = kpi_h_nat + slack * 0.20

    parts = [f'<g id="panel_d_{story_key}">']

    # outer subtle border
    parts.append(_rect(0, 0, box_w, box_h, fill="#ffffff",
                       stroke="#E5E7EB", sw=0.4))

    # ---- 1. Header banner -------------------------------------------------
    y = pad
    parts.append(_text(pad, y + 6.2, "d", size=8.5, weight="bold", color=INK))
    parts.append(_rect(pad + 6.5, y, 1.4, max(22, h_header - 2),
                       fill=story["accent"]))
    title_x = pad + 12
    parts.append(_text(title_x, y + 7.0, story["title"],
                       size=7.0, weight="bold", color=story["color"]))
    for i, line in enumerate(claim_lines):
        parts.append(_text(title_x, y + 14.4 + i * 7.0, line,
                           size=5.6, color=INK))
    y += h_header
    parts.append(_line(pad, y, box_w - pad, y, color="#E5E7EB", w=0.5))
    y += sep_h

    # ---- 2. Module gene chips ---------------------------------------------
    parts.append(_text(pad, y + 5.4, "Point module marker genes",
                       size=5.8, weight="bold", color=MUTED))
    y += 9
    chip_h = chip_h_nat
    cx = pad
    for g in genes:
        # ~3.1 pt per bold char at 5.8px + comfortable padding
        w = max(22.0, 3.10 * len(g) + 8.0)
        if cx + w > box_w - pad:
            cx = pad
            y += chip_h + 2.5
        parts.append(_rect(cx, y, w, chip_h, fill=story["tint"], rx=2.0,
                           stroke=story["color"], sw=0.4))
        parts.append(_text(cx + w / 2, y + chip_h - 2.6, g,
                           size=5.8, weight="bold", color=story["accent"],
                           anchor="middle"))
        cx += w + 3.0
    y += chip_h + 5
    parts.append(_line(pad, y, box_w - pad, y, color="#E5E7EB", w=0.5))
    y += sep_h

    # ---- 3. Top P->C matches ----------------------------------------------
    parts.append(_text(pad, y + 5.4,
                       f"Top {n_show} P→C matches (ranked by overlap)",
                       size=5.8, weight="bold", color=MUTED))
    y += 11

    chip_w = 34.0
    for r in matches[:n_show]:
        pid = int(r["point_domain_id"]); cid = int(r["best_cell_domain_id"])
        rel = r.get("pathology_label_relation", "different")
        rel_dark, rel_light, rel_name = RELATION_STYLE.get(rel, RELATION_STYLE["different"])
        # cell-side pathology label varies across rows -> more informative
        cell_label = r.get("cell_pathology_label") or r.get("matched_cell_dominant_celltype", "")
        if len(cell_label) > 28:
            cell_label = cell_label[:27] + "."
        celltype_label = _short(r["matched_cell_dominant_celltype"])
        frac = float(r["point_to_cell_fraction"])

        row_inner = match_row_h - 2.0
        # P->C chip
        parts.append(_rect(pad, y, chip_w, row_inner, fill=rel_light,
                           rx=2.0, stroke=rel_dark, sw=0.45))
        parts.append(_text(pad + chip_w / 2, y + row_inner * 0.65,
                           f"P{pid:02d}→C{cid:02d}",
                           size=6.4, weight="bold", color=rel_dark,
                           anchor="middle"))
        # cell pathology label (line 1) + celltype + overlap (line 2)
        text_x = pad + chip_w + 4.5
        parts.append(_text(text_x, y + 5.6, cell_label,
                           size=5.7, color=INK, weight="bold"))
        parts.append(_text(text_x, y + 11.0,
                           f"{celltype_label} · ovlp {frac:.2f}",
                           size=5.2, color=MUTED))
        # Relation is already encoded by the chip fill/stroke colour; the
        # redundant right-edge dot was removed and the colour key (same/related/
        # different) is defined in the figure legend instead.
        y += match_row_h
    y += 3
    parts.append(_line(pad, y, box_w - pad, y, color="#E5E7EB", w=0.5))
    y += sep_h

    # ---- 4. KPI tiles -----------------------------------------------------
    parts.append(_text(pad, y + 5.0, "Global comparison (whole section)",
                       size=5.4, weight="bold", color=MUTED))
    y += kpi_label_h
    kpis = [
        ("ARI",    f"{metrics['adjusted_rand_index']:.2f}",       "#42A5F5"),
        ("NMI",    f"{metrics['normalized_mutual_info']:.2f}",    "#7E57C2"),
        ("purity", f"{metrics['point_to_cell_purity']:.2f}",      "#26A69A"),
    ]
    avail = box_w - 2 * pad
    gap = 3.0
    tw = (avail - 2 * gap) / 3
    th = max(22.0, min(h_kpi, box_h - y - pad))
    for i, (k, v, c) in enumerate(kpis):
        x0 = pad + i * (tw + gap)
        parts.append(_rect(x0, y, tw, th, fill="#FAFAFA", rx=2.5,
                           stroke=c, sw=0.7))
        parts.append(_rect(x0, y, 2.0, th, fill=c, rx=0))
        parts.append(_text(x0 + tw / 2, y + th * 0.55, v,
                           size=11.0, weight="bold", color="#202830",
                           anchor="middle"))
        parts.append(_text(x0 + tw / 2, y + th - 3.0, k,
                           size=5.6, color=MUTED, anchor="middle"))

    parts.append('</g>')
    return "\n".join(parts) + "\n"


if __name__ == "__main__":
    # quick smoke-test render to a standalone SVG file
    import sys
    story = sys.argv[1] if len(sys.argv) > 1 else "tumor_rich"
    W, H = 150.0, 340.0
    block = make_evidence_card(story, W, H)
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}pt" '
           f'height="{H}pt" viewBox="0 0 {W} {H}">\n{block}</svg>\n')
    out = os.path.join(HERE, "_relayout_preview", f"panel_d_{story}.svg")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out, "w", encoding="utf-8", newline="\n").write(svg)
    print("wrote", out)
    # also render PNG for visual check
    import fitz
    from PIL import Image
    pix = fitz.open("svg", svg.encode("utf-8"))[0].get_pixmap(
        matrix=fitz.Matrix(4, 4), alpha=False)
    pix.save(out.replace(".svg", ".png"))
    print("wrote", out.replace(".svg", ".png"))
