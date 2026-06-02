"""Round-robin social-rating adapters for Cell-GPS-style analysis."""

from .round_robin_cellgps import (
    RoundRobinCellGPSResult,
    build_searcher_findee_matrix,
    compute_grouped_round_robin_cellgps,
    compute_round_robin_cellgps,
    score_matrix_to_distance_matrix,
    write_result_tables,
)

__all__ = [
    "RoundRobinCellGPSResult",
    "build_searcher_findee_matrix",
    "compute_grouped_round_robin_cellgps",
    "compute_round_robin_cellgps",
    "score_matrix_to_distance_matrix",
    "write_result_tables",
]
