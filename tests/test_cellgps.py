#!/usr/bin/env python

"""Basic smoke tests for the Cell-GPS public package surfaces."""

import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import cellgps
import cellgps.analysis
import cellgps.plotting
import cellgps.preprocessing
import sfplot
import sfplot.analysis
import sfplot.plotting
import sfplot.preprocessing


class TestCellGps(unittest.TestCase):
    """Minimal tests that validate the public package surface."""

    def test_package_imports(self):
        """The top-level package should import successfully."""
        self.assertTrue(hasattr(cellgps, "__all__"))
        self.assertGreater(len(cellgps.__all__), 0)

    def test_legacy_import_still_works(self):
        """The historical sfplot namespace should remain compatible."""
        self.assertEqual(set(cellgps.__all__), set(sfplot.__all__) | {"analysis", "preprocessing", "plotting", "gui"})

    def test_expected_public_exports_exist(self):
        """Key reviewer-facing APIs should remain exposed at the top level."""
        expected = {
            "load_xenium_data",
            "load_xenium_table_bundle",
            "compute_cophenetic_distances_from_df",
            "compute_cophenetic_distances_from_adata",
            "compute_weighted_searcher_findee_distance_matrix_from_df",
            "compute_weighted_cophenetic_distances_from_df",
            "compute_entity_to_cell_topology",
            "ligand_receptor_topology_analysis",
            "pathway_topology_analysis",
            "plot_cophenetic_heatmap",
            "transcript_by_cell_analysis",
            "pick_batch_size",
            "compute_groupwise_average_distance_between_two_dfs",
        }
        self.assertTrue(expected.issubset(set(cellgps.__all__)))

    def test_subpackage_preprocessing_exports(self):
        """cellgps.preprocessing should expose its public APIs."""
        required = {"load_xenium_data", "load_xenium_table_bundle", "read_visium_bin"}
        self.assertTrue(required.issubset(set(cellgps.preprocessing.__all__)))
        if hasattr(cellgps.preprocessing, "merge_xenium_clusters_into_adata"):
            self.assertIn("merge_xenium_clusters_into_adata", cellgps.preprocessing.__all__)

    def test_subpackage_analysis_exports(self):
        """cellgps.analysis should expose its public APIs."""
        expected = {
            "compute_cophenetic_distances_from_df",
            "compute_cophenetic_distances_from_adata",
            "pick_batch_size",
            "split_B_by_distance_to_A",
            "transcript_by_cell_analysis",
            "compute_weighted_searcher_findee_distance_matrix_from_df",
            "compute_entity_to_cell_topology",
            "ligand_receptor_topology_analysis",
            "pathway_topology_analysis",
        }
        self.assertTrue(expected.issubset(set(cellgps.analysis.__all__)))

    def test_subpackage_plotting_exports(self):
        """cellgps.plotting should expose its public APIs."""
        expected = {
            "circle_heatmap",
            "plot_cophenetic_heatmap",
        }
        self.assertTrue(expected.issubset(set(cellgps.plotting.__all__)))

    def test_backwards_compat_circular_dendrogram(self):
        """sfplot.circular_dendrogram backwards-compat alias should still work when pycirclize is installed."""
        try:
            import pycirclize  # noqa: F401
        except ImportError:
            self.skipTest("pycirclize not installed")
        self.assertTrue(hasattr(sfplot, "circular_dendrogram"))
        self.assertTrue(hasattr(sfplot.circular_dendrogram, "plot_circular_dendrogram_pycirclize"))


if __name__ == "__main__":
    unittest.main()
