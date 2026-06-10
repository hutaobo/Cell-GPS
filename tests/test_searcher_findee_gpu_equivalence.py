#!/usr/bin/env python

"""GPU/CPU equivalence tests for the searcher-findee / cophenetic pipeline.

The GPU backend must be numerically equivalent to the scikit-learn CPU backend.
These tests run the GPU code path on ``device="cpu"`` (so they execute on any
machine, including CI without a GPU) and assert that, in ``float64``, the GPU
result reproduces the CPU result exactly. The same assertions hold on real CUDA
hardware by passing ``device="cuda"``.
"""

from __future__ import annotations

import pathlib
import sys
import tempfile
import types
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

try:
    import torch  # noqa: F401

    HAS_TORCH = True
except Exception:  # pragma: no cover - exercised only without torch installed
    HAS_TORCH = False

from sfplot.analysis.searcher_findee_score import (
    compute_cophenetic_distances_from_adata,
    compute_cophenetic_distances_from_df,
    compute_searcher_findee_distance_matrix_from_df,
)


def _make_df(n=1500, k=8, dims=2, seed=0, scale=100.0):
    rng = np.random.default_rng(seed)
    coords = rng.uniform(0, scale, size=(n, dims))
    cols = {"x": coords[:, 0], "y": coords[:, 1]}
    if dims == 3:
        cols["z"] = coords[:, 2]
    cols["celltype"] = rng.integers(0, k, size=n).astype(str)
    return pd.DataFrame(cols)


def _fake_adata(n=1500, k=8, seed=1, scale=100.0):
    rng = np.random.default_rng(seed)
    obs = pd.DataFrame(
        {
            "Cluster": rng.integers(0, k, size=n).astype(str),
            "cell_id": [f"c{i}" for i in range(n)],
        }
    )
    adata = types.SimpleNamespace()
    adata.obs = obs
    adata.obsm = {"spatial": rng.uniform(0, scale, size=(n, 2))}
    adata.uns = {}
    return adata


def _assert_frame_equal_values(cpu, gpu, atol):
    assert list(cpu.index) == list(gpu.index), "row labels differ"
    assert list(cpu.columns) == list(gpu.columns), "column labels differ"
    np.testing.assert_allclose(
        cpu.values.astype(float),
        gpu.values.astype(float),
        rtol=0.0,
        atol=atol,
        equal_nan=True,
    )


@unittest.skipUnless(HAS_TORCH, "PyTorch is not installed")
class TestSearcherFindeeGpuEquivalence(unittest.TestCase):
    """GPU backend must reproduce the CPU backend (float64 == exact)."""

    EXACT = 1e-9  # float64 GPU reproduces scikit-learn to floating-point precision

    def test_searcher_findee_matrix_2d(self):
        df = _make_df(dims=2, seed=2)
        cpu = compute_searcher_findee_distance_matrix_from_df(df, backend="cpu")
        gpu = compute_searcher_findee_distance_matrix_from_df(
            df, backend="gpu", device="cpu", gpu_dtype="float64"
        )
        _assert_frame_equal_values(cpu, gpu, self.EXACT)

    def test_searcher_findee_matrix_3d(self):
        df = _make_df(dims=3, seed=3)
        cpu = compute_searcher_findee_distance_matrix_from_df(df, z_col="z", backend="cpu")
        gpu = compute_searcher_findee_distance_matrix_from_df(
            df, z_col="z", backend="gpu", device="cpu", gpu_dtype="float64"
        )
        _assert_frame_equal_values(cpu, gpu, self.EXACT)

    def test_cophenetic_from_df(self):
        df = _make_df(dims=2, seed=4)
        r_cpu, c_cpu = compute_cophenetic_distances_from_df(df, backend="cpu")
        r_gpu, c_gpu = compute_cophenetic_distances_from_df(
            df, backend="gpu", device="cpu", gpu_dtype="float64"
        )
        _assert_frame_equal_values(r_cpu, r_gpu, self.EXACT)
        _assert_frame_equal_values(c_cpu, c_gpu, self.EXACT)

    def test_cophenetic_from_adata(self):
        with tempfile.TemporaryDirectory() as tmp:
            adata = _fake_adata(seed=5)
            r_cpu, c_cpu = compute_cophenetic_distances_from_adata(
                adata, cluster_col="Cluster", output_dir=tmp, backend="cpu"
            )
            r_gpu, c_gpu = compute_cophenetic_distances_from_adata(
                adata, cluster_col="Cluster", output_dir=tmp,
                backend="gpu", device="cpu", gpu_dtype="float64",
            )
            _assert_frame_equal_values(r_cpu, r_gpu, self.EXACT)
            _assert_frame_equal_values(c_cpu, c_gpu, self.EXACT)

    def test_auto_backend_matches_cpu(self):
        # On a CPU-only machine "auto" resolves to CPU and must match exactly.
        df = _make_df(dims=2, seed=6)
        cpu = compute_searcher_findee_distance_matrix_from_df(df, backend="cpu")
        auto = compute_searcher_findee_distance_matrix_from_df(
            df, backend="auto", device="cpu", gpu_dtype="float64"
        )
        _assert_frame_equal_values(cpu, auto, self.EXACT)

    def test_float32_is_close(self):
        # The fast path trades exactness for speed; it should still be close.
        df = _make_df(dims=2, seed=7)
        cpu = compute_searcher_findee_distance_matrix_from_df(df, backend="cpu")
        gpu32 = compute_searcher_findee_distance_matrix_from_df(
            df, backend="gpu", device="cpu", gpu_dtype="float32"
        )
        _assert_frame_equal_values(cpu, gpu32, atol=1e-2)

    def test_unknown_backend_raises(self):
        df = _make_df(dims=2, seed=8)
        with self.assertRaises(ValueError):
            compute_searcher_findee_distance_matrix_from_df(df, backend="quantum")

    def test_public_gpu_export(self):
        import cellgps

        self.assertTrue(
            hasattr(cellgps, "compute_searcher_findee_distance_matrix_from_df_gpu")
        )


if __name__ == "__main__":
    unittest.main()
