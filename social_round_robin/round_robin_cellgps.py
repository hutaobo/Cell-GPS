"""Convert social round-robin ratings into Cell-GPS searcher-findee inputs.

This module is intentionally outside ``src/cellgps`` and ``src/sfplot``. It
keeps the sociology-facing adapter separate from the spatial-omics package
while reusing the Cell-GPS cophenetic matrix routine when available.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Hashable, Iterable, Optional, Union

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class RoundRobinCellGPSResult:
    """Container for one round-robin group or dataset."""

    affinity_matrix: pd.DataFrame
    distance_matrix: pd.DataFrame
    row_cophenetic: pd.DataFrame
    col_cophenetic: pd.DataFrame


def build_searcher_findee_matrix(
    data: pd.DataFrame,
    actor_col: str = "actor",
    partner_col: str = "partner",
    score_col: str = "score",
    group_col: Optional[str] = None,
    aggfunc: Union[str, Any] = "mean",
    make_square: bool = True,
    sort_labels: bool = True,
) -> Union[pd.DataFrame, Dict[Hashable, pd.DataFrame]]:
    """Build an actor/searcher by partner/findee score matrix.

    Parameters
    ----------
    data
        Long-form round-robin table. Each row is one directed judgment from
        ``actor_col`` to ``partner_col``.
    actor_col, partner_col, score_col
        Column names for rater, rated person, and rating.
    group_col
        Optional group identifier. If provided, returns one matrix per group.
    aggfunc
        Aggregation used when multiple rows exist for one actor-partner pair.
    make_square
        Reindex rows and columns to the same participant set. This is the
        natural layout for standard round-robin designs.
    sort_labels
        Sort participant labels for stable output.
    """

    required = [actor_col, partner_col, score_col]
    if group_col is not None:
        required.append(group_col)
    _require_columns(data, required)

    if group_col is not None:
        return {
            group: build_searcher_findee_matrix(
                group_data,
                actor_col=actor_col,
                partner_col=partner_col,
                score_col=score_col,
                group_col=None,
                aggfunc=aggfunc,
                make_square=make_square,
                sort_labels=sort_labels,
            )
            for group, group_data in data.groupby(group_col, sort=sort_labels)
        }

    work = data[[actor_col, partner_col, score_col]].copy()
    work = work.dropna(subset=[actor_col, partner_col])
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    work = work.dropna(subset=[score_col])
    if work.empty:
        raise ValueError("No valid round-robin ratings remain after filtering.")

    matrix = work.pivot_table(
        index=actor_col,
        columns=partner_col,
        values=score_col,
        aggfunc=aggfunc,
    )

    if make_square:
        labels = _ordered_union(matrix.index, matrix.columns, sort_labels=sort_labels)
        matrix = matrix.reindex(index=labels, columns=labels)
    elif sort_labels:
        matrix = matrix.sort_index(axis=0).sort_index(axis=1)

    matrix.index.name = "searcher"
    matrix.columns.name = "findee"
    return matrix.astype(float)


def score_matrix_to_distance_matrix(
    affinity_matrix: pd.DataFrame,
    high_score_is_close: bool = True,
    min_score: Optional[float] = None,
    max_score: Optional[float] = None,
    normalize: bool = True,
    fill_missing: Union[str, float, int, None] = "neutral",
    diagonal: Union[str, float, int, None] = "neutral",
) -> pd.DataFrame:
    """Convert a score/affinity matrix into a Cell-GPS distance matrix.

    Round-robin ratings are often affinities: larger values mean more liking,
    trust, similarity, or positive evaluation. Cell-GPS distance matrices use
    the opposite convention: smaller values mean closer searcher-findee
    topology. Set ``high_score_is_close=True`` for that common case.

    Missing self-ratings on the diagonal are filled with a neutral value by
    default. This avoids introducing a fake zero self-distance signal.
    """

    score = affinity_matrix.apply(pd.to_numeric, errors="coerce").astype(float)
    values = score.to_numpy(dtype=float)
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        raise ValueError("Affinity matrix has no finite scores.")

    observed_min = float(np.nanmin(finite_values)) if min_score is None else float(min_score)
    observed_max = float(np.nanmax(finite_values)) if max_score is None else float(max_score)
    if observed_max < observed_min:
        raise ValueError("max_score must be greater than or equal to min_score.")

    if observed_max == observed_min:
        converted = np.where(np.isfinite(values), 0.0, np.nan)
    elif high_score_is_close:
        converted = observed_max - values
        if normalize:
            converted = converted / (observed_max - observed_min)
    else:
        converted = values - observed_min
        if normalize:
            converted = converted / (observed_max - observed_min)

    distance = pd.DataFrame(
        converted,
        index=affinity_matrix.index.copy(),
        columns=affinity_matrix.columns.copy(),
        dtype=float,
    )
    return fill_distance_matrix(distance, fill_missing=fill_missing, diagonal=diagonal)


def fill_distance_matrix(
    distance_matrix: pd.DataFrame,
    fill_missing: Union[str, float, int, None] = "neutral",
    diagonal: Union[str, float, int, None] = "neutral",
) -> pd.DataFrame:
    """Fill missing distances in a prepared searcher-findee matrix."""

    filled = distance_matrix.copy().astype(float)

    if diagonal != "keep":
        for label in _shared_labels(filled):
            filled.loc[label, label] = np.nan

    filled = _fill_missing(filled, fill_missing)

    if diagonal != "keep":
        diag_value = _resolve_scalar_fill(filled, diagonal)
        for label in _shared_labels(filled):
            filled.loc[label, label] = diag_value

    return filled


def compute_round_robin_cellgps(
    data: pd.DataFrame,
    actor_col: str = "actor",
    partner_col: str = "partner",
    score_col: str = "score",
    high_score_is_close: bool = True,
    min_score: Optional[float] = None,
    max_score: Optional[float] = None,
    normalize: bool = True,
    fill_missing: Union[str, float, int, None] = "neutral",
    diagonal: Union[str, float, int, None] = "neutral",
    aggfunc: Union[str, Any] = "mean",
    method: str = "average",
    show_corr: bool = False,
) -> RoundRobinCellGPSResult:
    """Run the round-robin to Cell-GPS adapter for one group or dataset."""

    affinity = build_searcher_findee_matrix(
        data,
        actor_col=actor_col,
        partner_col=partner_col,
        score_col=score_col,
        group_col=None,
        aggfunc=aggfunc,
        make_square=True,
    )
    distance = score_matrix_to_distance_matrix(
        affinity,
        high_score_is_close=high_score_is_close,
        min_score=min_score,
        max_score=max_score,
        normalize=normalize,
        fill_missing=fill_missing,
        diagonal=diagonal,
    )
    if distance.isna().any().any():
        raise ValueError(
            "Distance matrix still contains NaN values. Use a finite fill_missing "
            "and diagonal strategy before Cell-GPS cophenetic analysis."
        )
    if distance.shape[0] < 2 or distance.shape[1] < 2:
        raise ValueError("At least two searchers and two findees are required.")

    row_coph, col_coph = _compute_cellgps_cophenetic(
        distance,
        method=method,
        show_corr=show_corr,
    )
    return RoundRobinCellGPSResult(
        affinity_matrix=affinity,
        distance_matrix=distance,
        row_cophenetic=row_coph,
        col_cophenetic=col_coph,
    )


def compute_grouped_round_robin_cellgps(
    data: pd.DataFrame,
    group_col: str,
    actor_col: str = "actor",
    partner_col: str = "partner",
    score_col: str = "score",
    **kwargs: Any,
) -> Dict[Hashable, RoundRobinCellGPSResult]:
    """Run the adapter separately for each round-robin group."""

    _require_columns(data, [group_col, actor_col, partner_col, score_col])
    results: Dict[Hashable, RoundRobinCellGPSResult] = {}
    for group, group_data in data.groupby(group_col, sort=True):
        results[group] = compute_round_robin_cellgps(
            group_data,
            actor_col=actor_col,
            partner_col=partner_col,
            score_col=score_col,
            **kwargs,
        )
    return results


def write_result_tables(
    result: RoundRobinCellGPSResult,
    output_dir: Union[str, Path],
    prefix: str = "round_robin",
) -> Dict[str, Path]:
    """Write affinity, distance, and cophenetic tables to CSV files."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        "affinity": out / f"{prefix}_searcher_findee_affinity.csv",
        "distance": out / f"{prefix}_searcher_findee_distance.csv",
        "row_cophenetic": out / f"{prefix}_searcher_cophenetic.csv",
        "col_cophenetic": out / f"{prefix}_findee_cophenetic.csv",
    }
    result.affinity_matrix.to_csv(paths["affinity"])
    result.distance_matrix.to_csv(paths["distance"])
    result.row_cophenetic.to_csv(paths["row_cophenetic"])
    result.col_cophenetic.to_csv(paths["col_cophenetic"])
    return paths


def _compute_cellgps_cophenetic(
    distance_matrix: pd.DataFrame,
    method: str,
    show_corr: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    try:
        from cellgps import compute_cophenetic_from_distance_matrix
    except ImportError as exc:
        raise ImportError(
            "Cell-GPS is required for cophenetic analysis. Install the package "
            "or run from this repository with PYTHONPATH=src."
        ) from exc

    return compute_cophenetic_from_distance_matrix(
        distance_matrix,
        method=method,
        show_corr=show_corr,
    )


def _require_columns(data: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = [column for column in columns if column not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _ordered_union(left: Iterable[Any], right: Iterable[Any], sort_labels: bool) -> list[Any]:
    labels = list(left) + list(right)
    if sort_labels:
        return sorted(set(labels), key=lambda value: str(value))

    seen = set()
    ordered = []
    for label in labels:
        if label not in seen:
            ordered.append(label)
            seen.add(label)
    return ordered


def _shared_labels(matrix: pd.DataFrame) -> list[Any]:
    return [label for label in matrix.index if label in matrix.columns]


def _finite_values(matrix: pd.DataFrame, exclude_diagonal: bool = False) -> np.ndarray:
    values = matrix.to_numpy(dtype=float).copy()
    if exclude_diagonal:
        row_lookup = {label: pos for pos, label in enumerate(matrix.index)}
        for col_pos, label in enumerate(matrix.columns):
            row_pos = row_lookup.get(label)
            if row_pos is not None:
                values[row_pos, col_pos] = np.nan
    return values[np.isfinite(values)]


def _resolve_scalar_fill(matrix: pd.DataFrame, strategy: Union[str, float, int, None]) -> float:
    if strategy is None or strategy == "nan":
        return np.nan
    if isinstance(strategy, (int, float)):
        return float(strategy)

    finite = _finite_values(matrix, exclude_diagonal=True)
    if finite.size == 0:
        finite = _finite_values(matrix, exclude_diagonal=False)
    if finite.size == 0:
        return 0.0

    if strategy == "neutral":
        return float(np.median(finite))
    if strategy == "max_distance":
        return float(np.max(finite))
    if strategy == "min_distance":
        return float(np.min(finite))
    if strategy == "zero":
        return 0.0
    raise ValueError(
        "Unknown fill strategy. Use neutral, row_mean, column_mean, "
        "max_distance, min_distance, zero, keep, nan, None, or a number."
    )


def _fill_missing(
    matrix: pd.DataFrame,
    strategy: Union[str, float, int, None],
) -> pd.DataFrame:
    if strategy is None:
        return matrix
    if isinstance(strategy, (int, float)):
        return matrix.fillna(float(strategy))
    if strategy in {"neutral", "max_distance", "min_distance", "zero", "nan"}:
        return matrix.fillna(_resolve_scalar_fill(matrix, strategy))
    if strategy == "row_mean":
        return _fill_by_axis_mean(matrix, axis=1)
    if strategy == "column_mean":
        return _fill_by_axis_mean(matrix, axis=0)
    raise ValueError(
        "Unknown fill_missing strategy. Use neutral, row_mean, column_mean, "
        "max_distance, min_distance, zero, nan, None, or a number."
    )


def _fill_by_axis_mean(matrix: pd.DataFrame, axis: int) -> pd.DataFrame:
    filled = matrix.copy()
    fallback = _resolve_scalar_fill(filled, "neutral")
    if axis == 1:
        means = filled.mean(axis=1, skipna=True)
        for label in filled.index:
            value = means.loc[label]
            filled.loc[label, :] = filled.loc[label, :].fillna(
                fallback if pd.isna(value) else value
            )
        return filled

    means = filled.mean(axis=0, skipna=True)
    for label in filled.columns:
        value = means.loc[label]
        filled.loc[:, label] = filled.loc[:, label].fillna(
            fallback if pd.isna(value) else value
        )
    return filled
