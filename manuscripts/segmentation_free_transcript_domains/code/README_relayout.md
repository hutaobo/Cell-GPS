# Supplementary figure relayout (Nature-style packing, evidence-card panel d)

`relayout_supplementary_figures.py` is a **post-processing** step that fixes the
panel packing, whitespace and typography of the five supplementary figures, and
**replaces the old text-only panel d on the pathology plates with a synthesised
visual evidence card**. It edits the rendered matplotlib SVGs and re-renders
PDF/PNG/SVG/TIFF with `fitz` (PyMuPDF). It does **not** regenerate figures from
data — the H&E OME-TIFF and the large transcript/cell parquet inputs live on
the cluster (`/data/taobo.hu/...`) and are unavailable locally, so layout is
fixed on the vector output instead. Helpers: `_svg_relayout.py` (geometry,
group transforms, text fonts), `_panel_d_evidence.py` (evidence card synth).

## Final supplementary figure set

| Final name           | Source SVG (pristine)                                                   | Story / content                                       |
| -------------------- | ----------------------------------------------------------------------- | ----------------------------------------------------- |
| Supplementary_Figure_1 | coste_histoseg_domains.svg                                              | Cell-based COSTE-HistoSeg reference map               |
| Supplementary_Figure_2 | cellfree_histoseg_domains.svg                                           | Superseded transcript-grid prototype                  |
| Supplementary_Figure_3 | nature_pathology_tumor_rich_plate.svg                                   | Tumor-rich pathology plate (was 3A)                   |
| Supplementary_Figure_4 | nature_pathology_immune_stroma_plate.svg                                | Immune-stromal interface pathology plate (was 3B)     |
| Supplementary_Figure_5 | nature_pathology_mixed_luminal_apocrine_plate.svg                       | Mixed luminal/apocrine context pathology plate (was 3C) |

## What it changes (content is never altered — only position/scale/typography)

* **Supplementary_Figure_3/4/5** (pathology plates):
  * Panels a (global H&E), b (zoom), c (overlap bars) are repacked as rigid
    blocks from their pristine matplotlib SVG to remove dead bottom-left
    whitespace. Supp 5 uses a tighter variant (a over c on the left; the two
    stacked zooms fill the right column).
  * **Panel d is replaced** with a synthesised "evidence card": story claim
    banner, point-side marker gene chips, top P->C matches coloured by
    pathology-label relation (same / related / different), and three KPI tiles
    (ARI / NMI / point-to-cell purity). Built from
    `figures/source_data/nature_pathology_source_domain_statistics.csv` and
    `figures/source_data/point_vs_cell_global_metrics.json`.
  * Small fonts raised to a Nature-compliant floor; canvas trimmed to content.
* **Supplementary_Figure_1/2** (scatter maps):
  * Baked-glyph fonts enlarged ~1.2x; colorbar widened ~1.45x (horizontal-only,
    anchored at its tick edge, frame + cells together); top/bottom letterbox
    whitespace cropped.

## Pipeline position (idempotent)

    cluster source scripts                                  (regenerate rarely)
        -> figures/supplementary_backup_orig_20260608/      <-- pristine input, NEVER overwritten
        -> patch_round1_pathology_svg_text.py               (caption text fixups)
        -> relayout_supplementary_figures.py apply           <-- THIS STEP (idempotent)
            * panel layout + evidence card panel d
            * writes figures/supplementary/*.{svg,pdf,png,tiff}
            * copies into the submission bundle as Supplementary_Figure_1..5
            * removes any legacy 3A/3B/3C files in the bundle
        -> _update_manuscript_supp_refs.py                   (one-off: docx caption edits)

The script reads from the **backup directory** (pristine matplotlib SVGs) and
writes to `figures/supplementary/`. Re-running is safe; the inputs are never
mutated.

## Usage

    python relayout_supplementary_figures.py preview            # -> code/_relayout_preview/ (no overwrite)
    python relayout_supplementary_figures.py preview plates     # only Supp 3/4/5
    python relayout_supplementary_figures.py preview scatter    # only Supp 1/2
    python relayout_supplementary_figures.py apply              # overwrite figures/supplementary, propagate to bundle

    python _update_manuscript_supp_refs.py                      # one-off docx caption edits

Backups created 2026-06-08:
* `figures/supplementary_backup_orig_20260608/`
* `review_package_round1/submission_RECOMMENDED_local_bundle/figures_backup_orig_20260608/`
* `review_package_round1/all_figure_pdfs_RECOMMENDED_backup_orig_20260608/`
* `Manuscript_submission_short_communication_BACKUP_20260608.docx`
* `Manuscript_round1_formatted_optionC_v5_BACKUP_20260608.docx`
