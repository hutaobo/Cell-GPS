"""Toy round-robin analysis using the sociology adapter."""

import pandas as pd

from social_round_robin import compute_round_robin_cellgps, write_result_tables


def main() -> None:
    ratings = pd.DataFrame(
        {
            "actor": ["A", "A", "B", "B", "C", "C", "D", "D"],
            "partner": ["B", "C", "A", "D", "A", "D", "B", "C"],
            "liking": [5, 1, 4, 2, 1, 5, 4, 2],
        }
    )

    result = compute_round_robin_cellgps(
        ratings,
        actor_col="actor",
        partner_col="partner",
        score_col="liking",
        min_score=1,
        max_score=5,
        high_score_is_close=True,
    )

    print("Affinity matrix")
    print(result.affinity_matrix)
    print("\nDistance matrix")
    print(result.distance_matrix)
    print("\nSearcher cophenetic")
    print(result.row_cophenetic)
    print("\nFindee cophenetic")
    print(result.col_cophenetic)

    write_result_tables(result, "social_round_robin/output", prefix="toy")


if __name__ == "__main__":
    main()
