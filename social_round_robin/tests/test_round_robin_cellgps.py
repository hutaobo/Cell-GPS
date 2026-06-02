import unittest

import numpy as np
import pandas as pd

from social_round_robin import (
    build_searcher_findee_matrix,
    compute_grouped_round_robin_cellgps,
    compute_round_robin_cellgps,
    score_matrix_to_distance_matrix,
)


class RoundRobinAdapterTests(unittest.TestCase):
    def test_build_searcher_findee_matrix_uses_actor_rows_partner_columns(self):
        ratings = pd.DataFrame(
            {
                "actor": ["A", "A", "B", "B", "C", "C"],
                "partner": ["B", "C", "A", "C", "A", "B"],
                "score": [5, 1, 4, 2, 1, 5],
            }
        )

        matrix = build_searcher_findee_matrix(ratings)

        self.assertEqual(list(matrix.index), ["A", "B", "C"])
        self.assertEqual(list(matrix.columns), ["A", "B", "C"])
        self.assertEqual(matrix.loc["A", "B"], 5)
        self.assertEqual(matrix.loc["B", "A"], 4)
        self.assertTrue(np.isnan(matrix.loc["A", "A"]))

    def test_high_score_is_close_converts_to_low_distance(self):
        affinity = pd.DataFrame(
            [[np.nan, 5, 1], [4, np.nan, 2], [1, 5, np.nan]],
            index=["A", "B", "C"],
            columns=["A", "B", "C"],
        )

        distance = score_matrix_to_distance_matrix(
            affinity,
            min_score=1,
            max_score=5,
            high_score_is_close=True,
            fill_missing="neutral",
            diagonal="neutral",
        )

        self.assertEqual(distance.loc["A", "B"], 0.0)
        self.assertEqual(distance.loc["A", "C"], 1.0)
        self.assertTrue(np.isfinite(distance.loc["A", "A"]))

    def test_end_to_end_returns_searcher_and_findee_cophenetics(self):
        ratings = pd.DataFrame(
            {
                "actor": ["A", "A", "B", "B", "C", "C", "D", "D"],
                "partner": ["B", "C", "A", "D", "A", "D", "B", "C"],
                "score": [5, 1, 4, 2, 1, 5, 4, 2],
            }
        )

        result = compute_round_robin_cellgps(
            ratings,
            min_score=1,
            max_score=5,
            high_score_is_close=True,
        )

        self.assertEqual(result.affinity_matrix.shape, (4, 4))
        self.assertEqual(result.distance_matrix.shape, (4, 4))
        self.assertEqual(result.row_cophenetic.shape, (4, 4))
        self.assertEqual(result.col_cophenetic.shape, (4, 4))

    def test_grouped_round_robin_runs_per_group(self):
        ratings = pd.DataFrame(
            {
                "group": ["g1"] * 6 + ["g2"] * 6,
                "actor": ["A", "A", "B", "B", "C", "C"] * 2,
                "partner": ["B", "C", "A", "C", "A", "B"] * 2,
                "score": [5, 1, 4, 2, 1, 5, 2, 5, 1, 4, 5, 1],
            }
        )

        results = compute_grouped_round_robin_cellgps(
            ratings,
            group_col="group",
            min_score=1,
            max_score=5,
        )

        self.assertEqual(set(results), {"g1", "g2"})
        self.assertEqual(results["g1"].distance_matrix.shape, (3, 3))


if __name__ == "__main__":
    unittest.main()
