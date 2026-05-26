from __future__ import annotations

import pathlib
import sys
import types

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

ad = pytest.importorskip("anndata")

from cellgps.analysis._xenium_io import normalize_transcript_points
from cellgps.preprocessing.data_processing import load_xenium_data, load_xenium_table_bundle


def _toy_adata(cluster_column: str = "cluster"):
    return ad.AnnData(
        X=np.asarray([[1.0, 0.0], [0.0, 2.0]]),
        obs=pd.DataFrame(
            {
                cluster_column: ["A", "B"],
                "x_centroid": [10.0, 20.0],
                "y_centroid": [30.0, 40.0],
            },
            index=["cell-a", "cell-b"],
        ),
        var=pd.DataFrame({"name": ["GeneA", "GeneB"]}, index=["gene-a", "gene-b"]),
    )


def _install_fake_pyxenium(monkeypatch, adata):
    calls = []
    io_module = types.ModuleType("pyXenium.io")

    def read_xenium(*args, **kwargs):
        calls.append((args, kwargs))
        return adata.copy()

    io_module.read_xenium = read_xenium
    package = types.ModuleType("pyXenium")
    package.io = io_module
    monkeypatch.setitem(sys.modules, "pyXenium", package)
    monkeypatch.setitem(sys.modules, "pyXenium.io", io_module)
    return calls


def test_load_xenium_data_delegates_to_pyxenium_and_preserves_cellgps_shape(monkeypatch):
    calls = _install_fake_pyxenium(monkeypatch, _toy_adata())

    adata = load_xenium_data("sample-xenium", normalize=False)

    assert calls[0][0] == ("sample-xenium",)
    assert calls[0][1]["as_"] == "anndata"
    assert calls[0][1]["include_transcripts"] is False
    assert calls[0][1]["include_boundaries"] is False
    assert adata.obs["cell_id"].tolist() == ["cell-a", "cell-b"]
    assert adata.obs["Cluster"].tolist() == ["A", "B"]
    np.testing.assert_allclose(adata.obsm["spatial"], [[10.0, 30.0], [20.0, 40.0]])
    assert adata.raw is not None


def test_load_xenium_data_normalizes_when_requested(monkeypatch):
    _install_fake_pyxenium(monkeypatch, _toy_adata())
    calls = []
    scanpy = types.SimpleNamespace(
        pp=types.SimpleNamespace(
            normalize_total=lambda adata, target_sum: calls.append(("normalize_total", target_sum)),
            log1p=lambda adata: calls.append(("log1p", None)),
            scale=lambda adata, max_value: calls.append(("scale", max_value)),
        )
    )
    monkeypatch.setitem(sys.modules, "scanpy", scanpy)

    load_xenium_data("sample-xenium", normalize=True)

    assert calls == [("normalize_total", 10000.0), ("log1p", None), ("scale", 10)]


def test_load_xenium_table_bundle_passes_explicit_artifacts_to_pyxenium(tmp_path, monkeypatch):
    calls = _install_fake_pyxenium(monkeypatch, _toy_adata("Clusters"))
    (tmp_path / "cells.parquet").write_text("", encoding="utf-8")
    (tmp_path / "sample_cell_groups.csv").write_text("", encoding="utf-8")
    (tmp_path / "cell_feature_matrix.h5").write_text("", encoding="utf-8")

    adata = load_xenium_table_bundle(tmp_path, normalize=False)

    kwargs = calls[0][1]
    assert kwargs["cells_path"] == tmp_path / "cells.parquet"
    assert kwargs["cell_groups_path"] == tmp_path / "sample_cell_groups.csv"
    assert kwargs["feature_matrix_path"] == tmp_path / "cell_feature_matrix.h5"
    assert kwargs["cluster_column_name"] == "Clusters"
    assert kwargs["clusters_relpath"] is None
    assert adata.obs["Clusters"].tolist() == ["A", "B"]
    assert adata.obs["Cluster"].tolist() == ["A", "B"]


def test_normalize_transcript_points_accepts_pyxenium_columns_without_cell_id():
    coords = pd.DataFrame(
        {
            "x": [1.0, 2.0, 3.0],
            "y": [4.0, 5.0, 6.0],
            "gene_name": ["GeneA", "NegControlProbe", "GeneB"],
        }
    )

    normalized = normalize_transcript_points(coords)

    assert normalized["feature_name"].tolist() == ["GeneA", "GeneB"]
    assert "cell_id" in normalized.columns
    assert normalized["cell_id"].isna().all()
