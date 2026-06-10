"""Relayout the segmentation-free supplementary figures for tighter, Nature-style
panel packing and larger fonts, operating on the rendered SVGs (source data is
not available locally) and re-rendering PDF/PNG/TIFF with fitz.

Run modes:
    python relayout_supplementary_figures.py preview   # render to _relayout_preview/
    python relayout_supplementary_figures.py apply      # overwrite finals (backups must exist)
"""

from __future__ import annotations

import os
import sys

import _svg_relayout as R
import _panel_d_evidence as PD

HERE = os.path.dirname(os.path.abspath(__file__))
# Pristine matplotlib output (never overwritten) is the canonical input.
# Outputs are written to the live figures/supplementary/ dir.
SRC = os.path.normpath(os.path.join(HERE, "..", "figures", "supplementary_backup_orig_20260608"))
SUPP = os.path.normpath(os.path.join(HERE, "..", "figures", "supplementary"))
PREVIEW = os.path.join(HERE, "_relayout_preview")

PLATES = {
    "nature_pathology_tumor_rich_plate":             ("two_panel_b", "tumor_rich"),
    "nature_pathology_immune_stroma_plate":          ("two_panel_b", "immune_stroma"),
    "nature_pathology_mixed_luminal_apocrine_plate": ("stacked_b",   "mixed_luminal_apocrine"),
}


def relayout_plate(stem: str, kind: str, story_key: str, out_dir: str,
                   make_tiff: bool = False) -> str:
    src = os.path.join(SRC, stem + ".svg")
    doc = R.load_svg(src)
    bb = R.measure_all(doc)
    gtext = {gid: t for gid, t in doc.children}

    W = doc.width
    m = 6.0
    sup_gid = [g for g, _ in doc.children if g.startswith("text_")][0]
    LEFT_R = 360.0       # right edge of the left (panel a/b/c) region
    D_X = 366.0          # left edge of the evidence card
    TOP = 24.0

    # Original d-text axes group, dropped from the new layout (replaced by the
    # generated evidence card). Kept here only so any potential downstream
    # consumer knows which group to ignore.
    drop_gid = "axes_5" if kind == "stacked_b" else "axes_4"  # noqa: F841

    placed: list[tuple[str, str]] = []
    ext = [1e9, 1e9, -1e9, -1e9]

    def add(gid, box, ax="center", ay="center", min_font=6.0, max_scale=None):
        tx, ty, s = R.fit_transform(bb[gid], box, ax, ay, max_scale=max_scale)
        placed.append((gid, R.place_group(gtext[gid], tx, ty, s, min_font_px=min_font)))
        x0, y0, x1, y1 = bb[gid]
        pb = (s * x0 + tx, s * y0 + ty, s * x1 + tx, s * y1 + ty)
        ext[0] = min(ext[0], pb[0]); ext[1] = min(ext[1], pb[1])
        ext[2] = max(ext[2], pb[2]); ext[3] = max(ext[3], pb[3])
        return pb

    # ---- panel a (global H&E): top-left ----
    a = add("axes_1", (m, TOP, 196.0, 168.0), ax="left", ay="top", min_font=6.0)

    if kind == "two_panel_b":
        b = add("axes_2", (202.0, TOP, LEFT_R, 168.0), ax="center", ay="top", min_font=6.0)
        row_bottom = max(a[3], b[3])
        c_box = (m, row_bottom + 6.0, LEFT_R, row_bottom + 6.0 + 188.0)
        add("axes_3", c_box, ax="center", ay="top", min_font=6.0)
    else:
        c = add("axes_4", (m, a[3] + 8.0, 188.0, a[3] + 8.0 + 196.0), ax="left", ay="top", min_font=6.0)
        col_x0, col_x1 = 200.0, LEFT_R
        top_h = (c[3] - TOP)
        add("axes_2", (col_x0, TOP, col_x1, TOP + top_h * 0.52), ax="center", ay="top", min_font=6.0)
        add("axes_3", (col_x0, TOP + top_h * 0.52 + 8.0, col_x1, c[3]), ax="center", ay="top", min_font=6.0)

    # ---- suptitle: top-left (placed BEFORE evidence card so it counts toward ext) ----
    add(sup_gid, (m, 3.0, LEFT_R, 19.0), ax="left", ay="center", min_font=8.0, max_scale=1.25)

    # ---- panel d: synthesised evidence card on the right ----
    d_box_w = (W - m) - D_X
    d_box_h = ext[3] - TOP                          # match height of left content
    card_inner = PD.make_evidence_card(story_key, d_box_w, d_box_h)
    card = (f'  <g id="panel_d_card" transform="translate({D_X:.4f} {TOP:.4f})">\n'
            f'{card_inner}'
            f'  </g>\n')
    placed.append(("panel_d_card", card))
    # extend ext to include the card so the canvas trims correctly
    ext[2] = max(ext[2], D_X + d_box_w)
    ext[3] = max(ext[3], TOP + d_box_h)

    # ---- finalize canvas to content extent + margin (auto-trim) ----
    Wn = ext[2] + m
    Hn = ext[3] + m
    bg = R.background_rect(Wn, Hn)
    children = [("patch_bg", bg)] + [(g, t if t.endswith("\n") else t + "\n") for g, t in placed]

    new_svg = R.serialize(doc, children=children, width=Wn, height=Hn, vb=(0.0, 0.0, Wn, Hn))
    os.makedirs(out_dir, exist_ok=True)
    stem_out = os.path.join(out_dir, stem)
    R.render_outputs(new_svg, stem_out, dpi=600, make_tiff=make_tiff)
    return stem_out + ".png"


# per-figure overrides:
#   cbar_style     = "bar" (12-cell categorical bar) or "circles" (one circle
#                    per data value + numeric label).
#   cbar_k         = horizontal scale (bar style only; ignored for circles).
#   cbar_height_k  = vertical scale relative to pristine bar height. Sets the
#                    total column height for either style.
# Supp 1/2 use the circle-legend style per user request 2026-06-08.
SCATTER = {
    "coste_histoseg_domains":    {"cbar_style": "circles", "cbar_height_k": 0.70, "font_k": 1.20},
    "cellfree_histoseg_domains": {"cbar_style": "circles", "cbar_height_k": 0.70, "font_k": 1.20},
}


# Colour-blind-safe 12-class categorical palette (domain id 1..12), selected by
# greedy max-min dispersion of CIELAB dE under normal/deuteranopia/protanopia
# vision (min pairwise dE 16.15). These are the exact colours the re-plotted
# domain maps use (ListedColormap+BoundaryNorm), so the circle legend matches
# the scatter. Replaces the old non-colour-blind-safe tab20 12-subset.
_TAB20_FOR_12_VALUES = [
    "#72190E", "#F7CB45", "#808000", "#90C987", "#117733", "#44AA99",
    "#88CCEE", "#4477AA", "#332288", "#F032E6", "#882255", "#EE3377",
]


def _replace_group(svg, gid, new_block):
    """Replace the named <g id="gid"> ... </g> group (and its leading indent
    + trailing newline) with `new_block`. `new_block` should be the complete
    replacement including its own indentation and trailing newline."""
    import re
    m = re.search(rf'( *)<g id="{re.escape(gid)}"[^>]*>', svg)
    if not m:
        return svg
    # locate matching </g> using generic tag depth counting
    depth = 0
    tag_re = re.compile(r'</?g\b[^>]*?(/?)>')
    i = m.start()
    while True:
        tm = tag_re.search(svg, i)
        if tm is None:
            return svg
        tag = tm.group(0)
        if tag.startswith('</g'):
            depth -= 1
            if depth == 0:
                end = tm.end()
                if end < len(svg) and svg[end] == '\n':
                    end += 1
                return svg[:m.start()] + new_block + svg[end:]
        elif tag.endswith('/>'):
            pass
        else:
            depth += 1
        i = tm.end()


def _redo_colorbar_as_12_cells(svg, cbar_k=1.15, height_k=0.70,
                                cell_colors=None, font_pt=8.0):
    """Replace matplotlib's default tab20-on-continuous-data colorbar
    (which renders as 20 cells with auto-placed ticks at value positions on a
    continuous scale, off-centre from any cell) with a true 12-cell categorical
    colorbar plus one numeric tick (1..12) at every cell centre.

    Applied to Supp Fig 1 and 2 (both `cmap='tab20'` over data values 1..12).
    Combines the horizontal widen (cbar_k) and vertical shrink (height_k)
    that were previously applied as separate transforms — the new bar is built
    from scratch at the final dimensions, so no group transforms are needed."""
    import re

    if cell_colors is None:
        cell_colors = _TAB20_FOR_12_VALUES
    n = len(cell_colors)

    # The matplotlib SVG uses random clipPath ids per render (e.g.
    # `pe0450804e0` for one figure, `p8367caef96` for the next). Identify the
    # colorbar's specific clipPath by reading the clip-path reference inside
    # QuadMesh_1 itself.
    qm_ref = re.search(r'<g id="QuadMesh_1">.*?clip-path="url\(#([^)]+)\)"',
                       svg, re.S)
    if not qm_ref:
        return svg
    clip_id = qm_ref.group(1)
    cm = re.search(
        rf'<clipPath id="{re.escape(clip_id)}">\s*<rect x="([0-9.]+)" y="([0-9.]+)" '
        rf'width="([0-9.]+)" height="([0-9.]+)"', svg)
    if not cm:
        return svg
    orig_xl = float(cm.group(1)); orig_y_top = float(cm.group(2))
    orig_w = float(cm.group(3)); orig_h = float(cm.group(4))
    orig_xr = orig_xl + orig_w
    orig_cy = orig_y_top + orig_h / 2.0

    # Final bar geometry: widen horizontally at right edge, shrink vertically
    # at vertical centre.
    new_w = orig_w * cbar_k
    new_xl = orig_xr - new_w
    new_xr = orig_xr
    new_h = orig_h * height_k
    new_y_top = orig_cy - new_h / 2.0
    cell_h = new_h / n

    # 1. clipPath
    svg = re.sub(
        rf'(<clipPath id="{clip_id}">\s*<rect x=")[0-9.]+(" y=")[0-9.]+'
        rf'(" width=")[0-9.]+(" height=")[0-9.]+',
        rf'\g<1>{new_xl:.5f}\g<2>{new_y_top:.5f}\g<3>{new_w:.5f}\g<4>{new_h:.5f}',
        svg, count=1)

    # 2. QuadMesh_1: 12 colour rects, top cell = v=12, bottom cell = v=1
    cell_rects = []
    for i in range(n):
        v = n - i
        y = new_y_top + i * cell_h
        color = cell_colors[v - 1]
        cell_rects.append(
            f'    <rect x="{new_xl:.5f}" y="{y:.5f}" width="{new_w:.5f}" '
            f'height="{cell_h:.5f}" fill="{color}" stroke="none"/>')
    new_qm = ('  <g id="QuadMesh_1">\n   <g clip-path="url(#' + clip_id + ')">\n'
              + '\n'.join(cell_rects)
              + '\n   </g>\n  </g>\n')
    svg = _replace_group(svg, 'QuadMesh_1', new_qm)

    # 3. matplotlib.axis_4: 12 ticks, each a short outward line + numeric label
    #    at the centre of its cell. Use real <text> for the label so it's not
    #    locked to matplotlib's baked-glyph font scaling.
    #
    #    IMPORTANT: matplotlib emits the rotated "histoseg domain" axis label
    #    (text_23) as the LAST child of matplotlib.axis_4, not as a direct
    #    child of axes_2. We extract it verbatim before replacing axis_4 so it
    #    keeps its baked-glyph rendering and its position at the bar centre.
    t23_open = re.search(r'( *)<g id="text_23"', svg)
    t23_block = ''
    if t23_open:
        t23_start = t23_open.start()
        depth = 0; i = t23_start
        for tm in re.finditer(r'</?g\b[^>]*?(/?)>', svg[i:]):
            tag = tm.group(0)
            if tag.startswith('</g'):
                depth -= 1
                if depth == 0:
                    end = i + tm.end()
                    if end < len(svg) and svg[end] == '\n':
                        end += 1
                    t23_block = svg[t23_start:end]
                    break
            elif tag.endswith('/>'):
                pass
            else:
                depth += 1

    tick_lines = []
    for i in range(n):
        v = n - i
        y_center = new_y_top + (i + 0.5) * cell_h
        tick_lines.append(
            f'   <line x1="{new_xr:.5f}" y1="{y_center:.5f}" '
            f'x2="{new_xr + 2.6:.5f}" y2="{y_center:.5f}" '
            f'stroke="#000000" stroke-width="0.6"/>')
        # nudge text y by ~0.36*font so it appears vertically centred
        tick_lines.append(
            f'   <text x="{new_xr + 4.0:.5f}" y="{y_center + font_pt * 0.36:.5f}" '
            f'font-family="Arial, Helvetica, sans-serif" '
            f'font-size="{font_pt}" fill="#000000">{v}</text>')
    new_axis4 = ('  <g id="matplotlib.axis_4">\n'
                 + '\n'.join(tick_lines)
                 + '\n'
                 + t23_block          # keep the rotated axis label
                 + '  </g>\n')
    svg = _replace_group(svg, 'matplotlib.axis_4', new_axis4)

    # 4. patch_8: replace with a clean rectangular frame at the new bar bounds.
    new_frame = (f'  <g id="patch_8">\n'
                 f'   <rect x="{new_xl:.5f}" y="{new_y_top:.5f}" '
                 f'width="{new_w:.5f}" height="{new_h:.5f}" '
                 f'fill="none" stroke="#000000" stroke-width="0.6"/>\n'
                 f'  </g>\n')
    svg = _replace_group(svg, 'patch_8', new_frame)

    # 5. patch_7 (axes background) is left alone but its lingering cbar_k
    #    transform from the previous step (if any) would now stretch a
    #    non-existent rectangle. Strip any transform attribute on patch_7.
    svg = re.sub(r'(<g id="patch_7")\s+transform="[^"]*"', r'\1', svg, count=1)

    # 6. text_23 (rotated "histoseg domain" label) stays at its original
    #    position which equals the unchanged bar centre.
    return svg


def _find_group_span(svg, start_idx):
    """Given an index pointing at a `<g ...>` opening tag, return (start, end)
    span (inclusive of opening + matching </g> + trailing newline)."""
    import re
    depth = 0; i = start_idx
    for m in re.finditer(r'</?g\b[^>]*?(/?)>', svg[i:]):
        tag = m.group(0)
        if tag.startswith('</g'):
            depth -= 1
            if depth == 0:
                end = i + m.end()
                if end < len(svg) and svg[end] == '\n':
                    end += 1
                return (start_idx, end)
        elif tag.endswith('/>'):
            pass
        else:
            depth += 1
    return (start_idx, len(svg))


def _redo_colorbar_as_circle_legend(svg, height_k=0.70, cell_colors=None,
                                     circle_r=3.8, label_font_pt=8.0,
                                     label_gap=3.5):
    """Replace the matplotlib colorbar with a clean vertical CIRCLE LEGEND:
    one filled coloured circle per data value (1..12) with the numeric label
    just to its right. Removes the bar's outline frame. Total column height
    follows `height_k` (relative to the pristine bar).

    Generalised over which axes_N contains the colorbar: detects it from
    QuadMesh_1's enclosing axes group, and preserves whatever axis-level label
    text (e.g. "histoseg domain", "domain id") is present rather than
    hard-coding a specific text_N."""
    import re

    if cell_colors is None:
        cell_colors = _TAB20_FOR_12_VALUES
    n = len(cell_colors)

    # 1. clip rect (random id per render — resolve via QuadMesh_1's reference).
    qm_match = re.search(r'<g id="QuadMesh_1">.*?clip-path="url\(#([^)]+)\)"', svg, re.S)
    if not qm_match:
        return svg
    clip_id = qm_match.group(1)
    cm = re.search(
        rf'<clipPath id="{re.escape(clip_id)}">\s*<rect x="([0-9.]+)" y="([0-9.]+)" '
        rf'width="([0-9.]+)" height="([0-9.]+)"', svg)
    if not cm:
        return svg
    orig_xl = float(cm.group(1)); orig_y_top = float(cm.group(2))
    orig_w = float(cm.group(3)); orig_h = float(cm.group(4))
    orig_cy = orig_y_top + orig_h / 2.0
    orig_cx = orig_xl + orig_w / 2.0
    new_h = orig_h * height_k
    new_y_top = orig_cy - new_h / 2.0
    row_h = new_h / n
    circle_cx = orig_cx
    label_x = circle_cx + circle_r + label_gap

    # 2. Identify the colorbar axes_N (enclosing QuadMesh_1) and the matching
    #    matplotlib.axis_K with the y-ticks. The axes_N may be axes_2 (single-
    #    panel scatter) or axes_4 (multi-panel scatter), etc.
    qm_pos = svg.index('<g id="QuadMesh_1">')
    axes_before = svg[:qm_pos]
    last_axes = max(((m.group(1), m.start()) for m in
                     re.finditer(r'<g id="(axes_\d+)">', axes_before)),
                    key=lambda x: x[1])
    cbar_axes_id = last_axes[0]
    ax_start = svg.index(f'<g id="{cbar_axes_id}">')
    _, ax_end = _find_group_span(svg, ax_start)
    ax_chunk = svg[ax_start:ax_end]

    # Pick the matplotlib.axis_K that has the ytick groups (i.e. the y-axis
    # for a vertical bar). Prefer the higher-numbered one (matplotlib's
    # convention: axis_{2N-1}=x, axis_{2N}=y).
    mp_ids = re.findall(r'<g id="(matplotlib\.axis_\d+)">', ax_chunk)
    cbar_axis_id = None
    for mid in sorted(mp_ids, key=lambda s: -int(s.split('_')[-1])):
        m_at = ax_chunk.index(f'<g id="{mid}">')
        _, m_end = _find_group_span(ax_chunk, m_at)
        mp_chunk = ax_chunk[m_at:m_end]
        if '<g id="ytick_' in mp_chunk:
            cbar_axis_id = mid
            break
    if cbar_axis_id is None:
        cbar_axis_id = mp_ids[-1] if mp_ids else None

    # 3. Extract axis-level label(s): any text_N that is a child of
    #    matplotlib.axis_K but NOT inside any ytick_M.
    label_blocks = []
    if cbar_axis_id:
        cb_pos = svg.index(f'<g id="{cbar_axis_id}">')
        _, cb_end = _find_group_span(svg, cb_pos)
        cb_chunk = svg[cb_pos:cb_end]
        # Build set of byte ranges occupied by ytick children
        ytick_ranges = []
        for tm in re.finditer(r'<g id="ytick_\d+"[^>]*>', cb_chunk):
            yr = _find_group_span(cb_chunk, tm.start())
            ytick_ranges.append(yr)
        for tm in re.finditer(r'<g id="text_\d+"', cb_chunk):
            tp = tm.start()
            if any(lo <= tp < hi for lo, hi in ytick_ranges):
                continue
            _, te = _find_group_span(cb_chunk, tp)
            label_blocks.append(cb_chunk[tp:te])

    # 4. Build new QuadMesh_1 (12 circles).
    circles = []
    for i in range(n):
        v = n - i
        cy = new_y_top + (i + 0.5) * row_h
        color = cell_colors[v - 1]
        circles.append(
            f'    <circle cx="{circle_cx:.5f}" cy="{cy:.5f}" r="{circle_r:.3f}" '
            f'fill="{color}" stroke="#202830" stroke-width="0.35"/>')
    new_qm = ('  <g id="QuadMesh_1">\n'
              + '\n'.join(circles)
              + '\n  </g>\n')
    svg = _replace_group(svg, 'QuadMesh_1', new_qm)

    # 5. Build new matplotlib.axis_K (12 labels + preserved axis label(s)).
    if cbar_axis_id:
        label_lines = []
        for i in range(n):
            v = n - i
            cy = new_y_top + (i + 0.5) * row_h
            label_lines.append(
                f'   <text x="{label_x:.5f}" y="{cy + label_font_pt * 0.36:.5f}" '
                f'font-family="Arial, Helvetica, sans-serif" '
                f'font-size="{label_font_pt}" fill="#000000">{v}</text>')
        new_axis = (f'  <g id="{cbar_axis_id}">\n'
                    + '\n'.join(label_lines)
                    + '\n' + ''.join(label_blocks) + '  </g>\n')
        svg = _replace_group(svg, cbar_axis_id, new_axis)

    # 6. clip rect: collapse to the circle column footprint.
    svg = re.sub(
        rf'(<clipPath id="{re.escape(clip_id)}">\s*<rect x=")[0-9.]+(" y=")[0-9.]+'
        rf'(" width=")[0-9.]+(" height=")[0-9.]+',
        rf'\g<1>{circle_cx - circle_r:.5f}\g<2>{new_y_top:.5f}'
        rf'\g<3>{2*circle_r:.5f}\g<4>{new_h:.5f}',
        svg, count=1)

    # 7. Remove the bar's outline frame and background patch transforms inside
    #    the colorbar axes — detect them dynamically rather than hard-coding ids.
    cbar_axes_start = svg.index(f'<g id="{cbar_axes_id}">')
    _, cbar_axes_end = _find_group_span(svg, cbar_axes_start)
    cbar_axes_chunk = svg[cbar_axes_start:cbar_axes_end]
    patch_ids = re.findall(r'<g id="(patch_\d+)"', cbar_axes_chunk)
    for pid in patch_ids:
        # If patch carries a stale transform from an earlier widen step, strip it.
        svg = re.sub(rf'(<g id="{pid}")\s+transform="[^"]*"', r'\1', svg, count=1)
    # The LAST patch in a colorbar axes is the outline frame — blank it.
    if patch_ids:
        frame_pid = patch_ids[-1]
        svg = _replace_group(svg, frame_pid, f'  <g id="{frame_pid}"/>\n')

    return svg


def _shrink_colorbar_height(svg, height_k):
    """Shrink the colorbar bar vertically by `height_k`, centred, and ALSO snap
    each tick to the centre of the colormap cell it labels.

    Why snap? matplotlib's `cmap="tab20"` on continuous data renders the bar
    with 20 cells while the data range is 1..12, so the auto-placed tick at
    value v is at the *value-fraction* position on a continuous scale and does
    NOT line up with any cell centre (off by ~5-9pt in the original output).
    The user-visible cells are 20 stacked rectangles; we want labels 2/4/6/...
    centred on the cells they correspond to.

    We parse the cell rectangles from QuadMesh_1, then for each y-tick in
    axes_2 we find the nearest cell centre (in the pristine coords) and
    translate the tick to that cell centre's post-shrink y position."""
    import re
    # Find the bar clip rect (defines bar y_top, height); also gives clip id.
    cm = re.search(r'<clipPath id="(pe[0-9a-f]+)">\s*<rect x="([0-9.]+)" y="([0-9.]+)" '
                   r'width="([0-9.]+)" height="([0-9.]+)"', svg)
    if not cm:
        return svg
    clip_id, _x, y_str, _w, h_str = cm.groups()
    y_top = float(y_str); h_orig = float(h_str)
    cy = y_top + h_orig / 2.0
    new_h = h_orig * height_k
    new_y_top = cy - new_h / 2.0
    ty = (1.0 - height_k) * cy

    # 1. Update the clipPath rect.
    svg = re.sub(
        rf'(<clipPath id="{clip_id}">\s*<rect x="[0-9.]+" y=")[0-9.]+(" width="[0-9.]+" height=")[0-9.]+',
        rf'\g<1>{new_y_top:.5f}\g<2>{new_h:.5f}',
        svg, count=1)

    # 2. Bar + frame: compose horizontal cbar_k transform with vertical shrink.
    def _augment(m):
        tx = float(m.group(2)); sx = float(m.group(3))
        return (f'<g id="{m.group(1)}" transform='
                f'"translate({tx} {ty:.5f}) scale({sx} {height_k:.5f})"')
    pat = re.compile(
        r'<g id="(QuadMesh_1|patch_7|patch_8)" transform="translate\(([-0-9.]+) 0\) '
        r'scale\(([0-9.]+) 1\)"')
    svg = pat.sub(_augment, svg)

    # 3. Parse cell centres from QuadMesh_1 paths (each cell is a closed rect).
    qm_start = svg.index('<g id="QuadMesh_1"')
    qm_end = svg.index('</g>', qm_start)
    cell_centres = []
    for m in re.finditer(
            r'<path d="M [0-9.]+ ([0-9.]+)\s*\nL [0-9.]+ \1\s*\nL [0-9.]+ ([0-9.]+)',
            svg[qm_start:qm_end]):
        y_bot, y_top_c = float(m.group(1)), float(m.group(2))
        cell_centres.append((y_bot + y_top_c) / 2.0)
    cell_centres.sort()

    # 4. For each colorbar y-tick, snap to the nearest cell centre.
    for tm in re.finditer(r'<g id="(ytick_\d+)">', svg):
        ytick_id = tm.group(1)
        before = svg[:tm.start()]
        if before.rfind('<g id="axes_2">') < before.rfind('<g id="axes_1">'):
            continue
        chunk = svg[tm.end():tm.end() + 900]
        um = re.search(r'<use[^>]*y="([0-9.]+)"', chunk)
        if not um:
            continue
        y_tick = float(um.group(1))
        # nearest cell centre in original (pre-shrink) coords
        target_orig = min(cell_centres, key=lambda c: abs(c - y_tick))
        # transform that cell centre into post-shrink coords
        target_new = height_k * target_orig + ty
        delta = target_new - y_tick
        svg = svg.replace(
            f'<g id="{ytick_id}">',
            f'<g id="{ytick_id}" transform="translate(0 {delta:.5f})">',
            1)
    return svg


def relayout_scatter(stem: str, out_dir: str, font_k: float = 1.20,
                     cbar_k: float = 1.45, cbar_height_k: float = 1.0,
                     cbar_tick_font_pt: float = 8.0,
                     cbar_style: str = "bar",
                     make_tiff: bool = False) -> str:
    """Tidy a single-panel scatter map: enlarge baked-glyph fonts (font_k),
    and REBUILD the colorbar as one of two categorical legend styles:
      * `cbar_style="bar"`     - 12 stacked colour rects with numeric ticks
      * `cbar_style="circles"` - vertical column of 12 colour circles with
                                  numeric labels (cbar_k is ignored).
    The default matplotlib output rendered tab20-on-continuous-data as 20 cells
    with ticks that didn't align to any cell. Both rebuild styles fix that."""
    import re
    src = os.path.join(SRC, stem + ".svg")
    svg = open(src, encoding="utf-8").read().replace("\r\n", "\n").replace("\r", "\n")

    # --- enlarge baked-glyph fonts in everything OUTSIDE the colorbar ---
    # (The colorbar contents get fully rebuilt below with real <text>, so any
    #  glyph-scale bump there would be overwritten anyway.)
    svg = R.bump_glyph_scale(svg, font_k)

    # --- rebuild the colorbar ---
    if cbar_style == "circles":
        svg = _redo_colorbar_as_circle_legend(
            svg, height_k=cbar_height_k, label_font_pt=cbar_tick_font_pt)
    else:
        svg = _redo_colorbar_as_12_cells(
            svg, cbar_k=cbar_k, height_k=cbar_height_k, font_pt=cbar_tick_font_pt)

    # --- crop letterbox whitespace ---
    x0, y0, x1, y1 = R.content_bbox_full(svg, scale=2.0)
    pad = 4.0
    svg = R.set_canvas(svg, x0 - pad, y0 - pad, (x1 - x0) + 2 * pad, (y1 - y0) + 2 * pad)

    os.makedirs(out_dir, exist_ok=True)
    stem_out = os.path.join(out_dir, stem)
    R.render_outputs(svg, stem_out, dpi=600, make_tiff=make_tiff)
    return stem_out + ".png"


# source stem -> formal submission name (3A/B/C are now separate figures 3/4/5)
FORMAL = {
    "coste_histoseg_domains": "Supplementary_Figure_1",
    "cellfree_histoseg_domains": "Supplementary_Figure_2",
    "nature_pathology_tumor_rich_plate": "Supplementary_Figure_3",
    "nature_pathology_immune_stroma_plate": "Supplementary_Figure_4",
    "nature_pathology_mixed_luminal_apocrine_plate": "Supplementary_Figure_5",
}
# old formal names that may still exist in the package and need to be removed
LEGACY_FORMAL = {"Supplementary_Figure_3A", "Supplementary_Figure_3B", "Supplementary_Figure_3C"}
RP = os.path.normpath(os.path.join(HERE, "..", "review_package_round1"))
BUNDLE_FIGS = os.path.join(RP, "submission_RECOMMENDED_local_bundle", "figures")
ALL_PDFS = os.path.join(RP, "all_figure_pdfs_RECOMMENDED")


def propagate():
    """Copy figures/supplementary -> bundle + all_figure_pdfs. Skips files that
    can't be copied (e.g. locked by an open Word/PDF viewer) instead of
    crashing the run; callers can re-invoke once the lock clears."""
    import shutil
    skipped = []

    def _safe_copy(src, dst):
        try:
            shutil.copy2(src, dst)
            return True
        except PermissionError as e:
            skipped.append((dst, str(e)))
            return False

    # 1. Remove legacy Supplementary_Figure_3A/3B/3C files from the package.
    for stale in LEGACY_FORMAL:
        for d in (BUNDLE_FIGS, ALL_PDFS):
            for ext in ("svg", "pdf", "png", "tiff"):
                p = os.path.join(d, f"{stale}.{ext}")
                if os.path.exists(p):
                    try:
                        os.remove(p); print("removed legacy", os.path.basename(p))
                    except PermissionError:
                        skipped.append((p, "locked"))
    # 2. Copy new formal files.
    for src_stem, formal in FORMAL.items():
        for ext in ("svg", "pdf", "png", "tiff"):
            src = os.path.join(SUPP, f"{src_stem}.{ext}")
            if os.path.exists(src):
                _safe_copy(src, os.path.join(BUNDLE_FIGS, f"{formal}.{ext}"))
        pdf = os.path.join(BUNDLE_FIGS, f"{formal}.pdf")
        if os.path.exists(pdf):
            _safe_copy(pdf, os.path.join(ALL_PDFS, f"{formal}.pdf"))
        print("propagated", formal)

    if skipped:
        print()
        print(f"WARNING: {len(skipped)} file(s) skipped (likely open in another program):")
        for path, reason in skipped:
            print(f"  - {os.path.basename(path)}: {reason.splitlines()[0]}")
        print("Close the program and re-run `python relayout_supplementary_figures.py apply ...` to retry.")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "preview"
    target = sys.argv[2] if len(sys.argv) > 2 else "all"
    out_dir = PREVIEW if mode == "preview" else SUPP
    if target in ("all", "plates"):
        for stem, (kind, story_key) in PLATES.items():
            print("wrote", relayout_plate(stem, kind, story_key, out_dir,
                                          make_tiff=(mode == "apply")))
    if target in ("all", "scatter"):
        for stem, opts in SCATTER.items():
            print("wrote", relayout_scatter(stem, out_dir,
                                            make_tiff=(mode == "apply"), **opts))
    if mode == "apply":
        propagate()


if __name__ == "__main__":
    main()
