"""SVG relayout utilities for the segmentation-free supplementary figures.

The source figures are generated on a cluster (H&E OME-TIFF + large parquet data
that are not available locally), so layout/typography fixes are done by operating
on the already-rendered matplotlib SVGs and re-rendering PDF/PNG/TIFF with fitz.

Design facts this relies on (verified against the matplotlib 3.10 SVG output):
  * <svg> has two <defs> children (style + clipPaths) and one <g id="figure_1">.
  * figure_1's direct children are indented with exactly 2 spaces and are:
    patch_1 (white background), one <g id="axes_N"> per panel, and optional
    top-level <g id="text_N"> (suptitle). Each axes group is fully self-contained
    (image, ticks, title, scale bar, panel letter, table text all inside it),
    so a panel can be moved/scaled as a rigid block by wrapping it in a
    <g transform="..."> without touching its internal composition.
  * Text uses svg.fonttype:none, i.e. real <text> elements with an explicit
    `font-size: Npx` in their style attribute -> editable by regex.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import fitz
import numpy as np


# ---------------------------------------------------------------------------
# Parsing: split the SVG into header / figure children / trailer using the
# reliable 2-space indentation of figure_1's direct children.
# ---------------------------------------------------------------------------


@dataclass
class SvgDoc:
    path: str
    width: float
    height: float
    vb: tuple[float, float, float, float]
    header: str           # everything up to and including '<g id="figure_1">\n'
    children: list[tuple[str, str]]  # (id, raw_text_including_trailing_newline)
    figure_close: str     # ' </g>\n'
    trailer: str          # second <defs>...</defs> + '</svg>'
    raw: str = field(default="", repr=False)


def load_svg(path: str) -> SvgDoc:
    raw = open(path, "r", encoding="utf-8").read().replace("\r\n", "\n").replace("\r", "\n")

    m = re.search(r'<svg[^>]*\bwidth="([0-9.]+)pt"[^>]*\bheight="([0-9.]+)pt"[^>]*viewBox="([^"]+)"', raw)
    if not m:
        raise ValueError(f"could not parse <svg> header in {path}")
    width = float(m.group(1))
    height = float(m.group(2))
    vb = tuple(float(v) for v in m.group(3).replace(",", " ").split())

    fig_open_m = re.search(r'\n( *)<g id="figure_1">\n', raw)
    if not fig_open_m:
        raise ValueError("no figure_1 group")
    fig_body_start = fig_open_m.end()
    header = raw[:fig_body_start]

    # figure_1's direct children open with exactly two leading spaces.
    child_open_re = re.compile(r'^  <g id="([^"]+)">$', re.M)

    children: list[tuple[str, str]] = []
    pos = fig_body_start
    while True:
        om = child_open_re.search(raw, pos)
        if om is None:
            break
        gid = om.group(1)
        end = _match_group_end(raw, om.start())   # robust generic <g>/</g> matching
        children.append((gid, raw[om.start():end]))
        pos = end

    if not children:
        raise ValueError(f"no figure_1 children parsed in {path}")

    # figure_1 close is the 1-space-indented ' </g>' right after the last child
    # (exactly one leading space distinguishes it from a 2-space child close).
    fig_close_m = re.compile(r'^ </g>\n', re.M).search(raw, pos)
    if not fig_close_m:
        raise ValueError("no figure_1 close")
    figure_close = raw[fig_close_m.start():fig_close_m.end()]
    trailer = raw[fig_close_m.end():]

    return SvgDoc(path, width, height, vb, header, children, figure_close, trailer, raw)


def _match_group_end(raw: str, start: int) -> int:
    """Given index of a '<g' opening, return index just past its matching '</g>'."""
    i = start
    depth = 0
    tag_re = re.compile(r'</?g\b[^>]*?(/?)>')
    while True:
        m = tag_re.search(raw, i)
        if m is None:
            raise ValueError("unbalanced <g>")
        tag = m.group(0)
        if tag.startswith('</g'):
            depth -= 1
            if depth == 0:
                end = m.end()
                # include trailing newline
                if end < len(raw) and raw[end] == '\n':
                    end += 1
                return end
        elif tag.endswith('/>'):
            pass  # self-closed, no depth change
        else:
            depth += 1
        i = m.end()


def serialize(doc: SvgDoc, children: list[tuple[str, str]] | None = None,
              width: float | None = None, height: float | None = None,
              vb: tuple[float, float, float, float] | None = None) -> str:
    if children is None:
        children = doc.children
    header = doc.header
    trailer = doc.trailer
    if width is not None and height is not None and vb is not None:
        header = re.sub(r'width="[0-9.]+pt"', f'width="{width:.4f}pt"', header, count=1)
        header = re.sub(r'height="[0-9.]+pt"', f'height="{height:.4f}pt"', header, count=1)
        header = re.sub(r'viewBox="[^"]+"', f'viewBox="{vb[0]:.4f} {vb[1]:.4f} {vb[2]:.4f} {vb[3]:.4f}"', header, count=1)
    return header + "".join(t for _, t in children) + doc.figure_close + trailer


# ---------------------------------------------------------------------------
# Measurement: render one panel group in isolation and trim to its bbox.
# ---------------------------------------------------------------------------


def measure_group_bbox(doc: SvgDoc, gid: str, scale: float = 4.0,
                       alpha_thresh: int = 8) -> tuple[float, float, float, float]:
    """Return (x0, y0, x1, y1) visible bbox of group `gid` in pt (svg user units)."""
    target = [(cid, t) for cid, t in doc.children if cid == gid]
    if not target:
        raise KeyError(gid)
    svg = serialize(doc, children=target)
    pix = fitz.open("svg", svg.encode("utf-8"))[0].get_pixmap(
        matrix=fitz.Matrix(scale, scale), alpha=True)
    arr = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)
    alpha = arr[..., 3]
    mask = alpha > alpha_thresh
    if not mask.any():
        return (0.0, 0.0, 0.0, 0.0)
    ys, xs = np.where(mask)
    # pixmap origin corresponds to viewBox origin scaled; convert back to pt
    vx, vy = doc.vb[0], doc.vb[1]
    x0 = vx + xs.min() / scale
    x1 = vx + (xs.max() + 1) / scale
    y0 = vy + ys.min() / scale
    y1 = vy + (ys.max() + 1) / scale
    return (float(x0), float(y0), float(x1), float(y1))


def measure_all(doc: SvgDoc, scale: float = 4.0) -> dict[str, tuple]:
    out = {}
    for gid, _ in doc.children:
        if gid == "patch_1":
            continue
        out[gid] = measure_group_bbox(doc, gid, scale=scale)
    return out


# ---------------------------------------------------------------------------
# Editing: wrap a panel group in a translate+scale transform (rigid block move)
# and bump small fonts so the post-scale effective size meets a floor.
# ---------------------------------------------------------------------------


_FONT_RE = re.compile(r'font-size:\s*([0-9.]+)px')


def bump_fonts(group_text: str, min_effective_px: float, scale: float) -> str:
    """Raise raw font-size values so that raw*scale >= min_effective_px."""
    raw_floor = min_effective_px / scale

    def repl(m):
        val = float(m.group(1))
        if val < raw_floor:
            val = raw_floor
        return f"font-size: {val:.3f}px"

    return _FONT_RE.sub(repl, group_text)


def fit_transform(bbox, box, align_x="center", align_y="center",
                  max_scale: float | None = None):
    """translate+scale mapping the group's measured `bbox` into target `box`.

    bbox, box are (x0, y0, x1, y1). Uniform scale (contain). Returns (tx, ty, s).
    """
    bx0, by0, bx1, by1 = bbox
    bw, bh = bx1 - bx0, by1 - by0
    Bx0, By0, Bx1, By1 = box
    Bw, Bh = Bx1 - Bx0, By1 - By0
    s = min(Bw / bw, Bh / bh)
    if max_scale is not None:
        s = min(s, max_scale)
    # placed size
    pw, ph = s * bw, s * bh
    if align_x == "left":
        ox = Bx0
    elif align_x == "right":
        ox = Bx1 - pw
    else:
        ox = Bx0 + (Bw - pw) / 2
    if align_y == "top":
        oy = By0
    elif align_y == "bottom":
        oy = By1 - ph
    else:
        oy = By0 + (Bh - ph) / 2
    tx = ox - s * bx0
    ty = oy - s * by0
    return tx, ty, s


def place_group(group_text: str, tx: float, ty: float, s: float,
                min_font_px: float | None = None) -> str:
    """Insert transform into the group's opening tag; optionally bump fonts."""
    if min_font_px is not None and s > 0:
        group_text = bump_fonts(group_text, min_font_px, s)
    # opening tag is '  <g id="...">' -> add transform attribute
    new = re.sub(
        r'(<g id="[^"]+")>',
        rf'\1 transform="translate({tx:.4f} {ty:.4f}) scale({s:.5f})">',
        group_text, count=1)
    return new


def background_rect(width: float, height: float) -> str:
    return (f'  <g id="patch_bg">\n'
            f'   <path d="M 0 0 L {width:.3f} 0 L {width:.3f} {height:.3f} '
            f'L 0 {height:.3f} Z" style="fill: #ffffff"/>\n'
            f'  </g>\n')


def parse_viewbox(svg_text: str) -> tuple[float, float, float, float]:
    m = re.search(r'viewBox="([^"]+)"', svg_text)
    return tuple(float(v) for v in m.group(1).replace(",", " ").split())


def set_canvas(svg_text: str, x0: float, y0: float, w: float, h: float) -> str:
    svg_text = re.sub(r'width="[0-9.]+pt"', f'width="{w:.4f}pt"', svg_text, count=1)
    svg_text = re.sub(r'height="[0-9.]+pt"', f'height="{h:.4f}pt"', svg_text, count=1)
    svg_text = re.sub(r'viewBox="[^"]+"',
                      f'viewBox="{x0:.4f} {y0:.4f} {w:.4f} {h:.4f}"', svg_text, count=1)
    return svg_text


def content_bbox_full(svg_text: str, scale: float = 2.0, white_thresh: int = 248):
    """Content bbox (pt) of a whole SVG by rendering on white and trimming
    near-white pixels (robust to an opaque white background rect)."""
    vb = parse_viewbox(svg_text)
    pix = fitz.open("svg", svg_text.encode("utf-8"))[0].get_pixmap(
        matrix=fitz.Matrix(scale, scale), alpha=False)
    arr = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)[..., :3]
    mask = np.any(arr < white_thresh, axis=2)
    if not mask.any():
        return vb[0], vb[1], vb[0] + vb[2], vb[1] + vb[3]
    ys, xs = np.where(mask)
    return (vb[0] + xs.min() / scale, vb[1] + ys.min() / scale,
            vb[0] + (xs.max() + 1) / scale, vb[1] + (ys.max() + 1) / scale)


def bump_glyph_scale(svg_text: str, k: float) -> str:
    """Enlarge baked-glyph text by multiplying its `scale(0.X -0.X)` transforms."""
    def repl(m):
        a, b = float(m.group(1)), float(m.group(2))   # b is negative (glyph y-flip)
        if abs(a) < 0.5:
            return f"scale({a * k:.6f} {b * k:.6f})"
        return m.group(0)
    return re.sub(r'scale\((-?0\.\d+) (-0\.\d+)\)', repl, svg_text)


def render_outputs(svg_text: str, out_stem: str, dpi: int = 600,
                   make_tiff: bool = True) -> list[str]:
    """Write {stem}.svg/.pdf/.png(/.tiff) from svg_text using fitz."""
    from PIL import Image

    svg_path = out_stem + ".svg"
    open(svg_path, "w", encoding="utf-8", newline="\n").write(svg_text)
    doc = fitz.open("svg", svg_text.encode("utf-8"))
    page = doc[0]
    open(out_stem + ".pdf", "wb").write(doc.convert_to_pdf())
    pix = page.get_pixmap(matrix=fitz.Matrix(dpi / 72, dpi / 72), alpha=False)
    tmp = out_stem + ".tmp.png"
    pix.save(tmp)
    with Image.open(tmp) as im:
        rgb = im.convert("RGB")
        rgb.save(out_stem + ".png", dpi=(dpi, dpi))
        if make_tiff:
            rgb.save(out_stem + ".tiff", dpi=(dpi, dpi), compression="tiff_lzw")
    import os
    os.remove(tmp)
    outs = [out_stem + e for e in (".svg", ".pdf", ".png") + ((".tiff",) if make_tiff else ())]
    return outs


if __name__ == "__main__":
    import sys
    p = sys.argv[1]
    d = load_svg(p)
    print(f"canvas: {d.width:.2f} x {d.height:.2f} pt  viewBox={d.vb}")
    print(f"children: {[c for c, _ in d.children]}")
    bb = measure_all(d)
    for gid, (x0, y0, x1, y1) in bb.items():
        print(f"  {gid:10s} x[{x0:7.2f},{x1:7.2f}] w={x1-x0:7.2f}  y[{y0:7.2f},{y1:7.2f}] h={y1-y0:7.2f}")
