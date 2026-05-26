from __future__ import annotations

import importlib

import pandas as pd


def load_xenium_transcripts(folder: str) -> pd.DataFrame:
    try:
        px_io = importlib.import_module("pyXenium.io")
    except ImportError as exc:
        raise ImportError(
            "pyXenium>=0.4.3 is required for transcript-by-cell analysis. "
            "Install Cell-GPS[xenium] or install pyXenium in the active environment."
        ) from exc

    sdata = px_io.read_xenium(
        folder,
        as_="sdata",
        include_transcripts=True,
        include_boundaries=False,
        include_images=False,
    )
    if "transcripts" in getattr(sdata, "points", {}):
        coords = sdata.points["transcripts"]
    elif "transcripts" in getattr(sdata, "point_sources", {}):
        coords = sdata.point_sources["transcripts"].materialize()
    else:
        raise ValueError(f"No Xenium transcripts were found under {folder!r}.")
    return normalize_transcript_points(coords)


def normalize_transcript_points(coords) -> pd.DataFrame:
    try:
        coords = coords.compute()
    except AttributeError:
        pass

    frame = coords.copy()
    feature_col = next(
        (column for column in ("feature_name", "gene_name", "gene_identity") if column in frame.columns),
        None,
    )
    if feature_col is None:
        raise ValueError("Transcript table must contain feature_name, gene_name, or gene_identity.")
    missing = [column for column in ("x", "y") if column not in frame.columns]
    if missing:
        raise ValueError(f"Transcript table is missing coordinate column(s): {missing}")

    frame["feature_name"] = frame[feature_col].astype(str)
    if "cell_id" not in frame.columns:
        frame["cell_id"] = pd.NA

    frame = frame.loc[
        ~frame["feature_name"].str.contains("NegControl|Unassigned", na=False),
        ["x", "y", "feature_name", "cell_id"],
    ]
    return frame.reset_index(drop=True)


def iter_cellgps_gene_names(adata) -> list[str]:
    if "name" in adata.var.columns:
        return adata.var["name"].astype(str).tolist()
    return adata.var_names.astype(str).tolist()


iter_sfplot_gene_names = iter_cellgps_gene_names
