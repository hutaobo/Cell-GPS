from __future__ import annotations


def test_scverse_namespace_aliases_are_available():
    import cellgps
    import cellgps.pl
    import cellgps.pp
    import cellgps.tl

    assert cellgps.pp is cellgps.pp
    assert cellgps.tl is cellgps.tl
    assert cellgps.pl is cellgps.pl

    assert "load_xenium_data" in cellgps.pp.__all__
    assert "compute_cophenetic_distances_from_df" in cellgps.tl.__all__
    assert "plot_cophenetic_heatmap" in cellgps.pl.__all__
