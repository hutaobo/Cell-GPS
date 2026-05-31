#!/usr/bin/env python
"""Verify + benchmark the GPU searcher-findee kernel against the CPU kernel.

This exercises the real Cell-GPS pipeline kernel
(``compute_searcher_findee_distance_matrix_from_df``) through its ``backend``
switch, so the numbers reflect the shipped code path, not a re-implementation.

It does the two things that need a real GPU (which CI / a CPU-only box cannot):

  [1] EQUIVALENCE  - confirms the GPU result equals the CPU result. In float64
                     this is exact (max|delta| ~ 0); float32 is the fast mode
                     (~1e-5 relative deviation).
  [2] BENCHMARK    - measures CPU-vs-GPU wall time across input sizes with proper
                     GPU timing (warm-up + torch.cuda.synchronize). The crossover,
                     if any, is where brute-force GPU cdist (O(N^2)) loses to the
                     CPU KD-tree (O(N log N)) -> that is when a cuML/cuVS index is
                     worth adding.

Examples
--------
    python examples/gpu_searcher_findee_benchmark.py
    python examples/gpu_searcher_findee_benchmark.py --max-n 1000000 --k 30 --dims 2
    python examples/gpu_searcher_findee_benchmark.py --dims 3 --dtype float32
    python examples/gpu_searcher_findee_benchmark.py --json results.json

Set --k (number of cell types) and --dims (coordinate dimensionality) to match
your real data so the timing reflects your workload.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd

try:
    from cellgps import compute_searcher_findee_distance_matrix_from_df as sf
except Exception:  # running from a source checkout without an install
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
    from cellgps import compute_searcher_findee_distance_matrix_from_df as sf

import torch


def make_df(n: int, k: int, dims: int, seed: int = 0, scale: float = 1000.0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    coords = rng.uniform(0, scale, size=(n, dims))
    cols = {"x": coords[:, 0], "y": coords[:, 1]}
    if dims == 3:
        cols["z"] = coords[:, 2]
    cols["celltype"] = rng.integers(0, k, size=n).astype(str)
    return pd.DataFrame(cols)


def max_abs_diff(a: pd.DataFrame, b: pd.DataFrame) -> float:
    a = a.reindex(index=b.index, columns=b.columns)
    av, bv = a.values.astype(float), b.values.astype(float)
    both = ~(np.isnan(av) | np.isnan(bv))
    return float(np.abs(av[both] - bv[both]).max()) if both.any() else 0.0


def time_call(fn, reps: int, device: str) -> float:
    fn()  # warm-up (kernel compile, allocator, caches)
    if device.startswith("cuda"):
        torch.cuda.synchronize()
    t0 = time.perf_counter()
    for _ in range(reps):
        fn()
    if device.startswith("cuda"):
        torch.cuda.synchronize()
    return (time.perf_counter() - t0) / reps


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--max-n", type=int, default=100_000, help="largest cell count to benchmark")
    p.add_argument("--k", type=int, default=10, help="number of cell types")
    p.add_argument("--dims", type=int, default=2, choices=(2, 3), help="coordinate dimensionality")
    p.add_argument("--dtype", default="float64", choices=("float64", "float32"),
                   help="GPU dtype for the benchmark (float64 = exact, float32 = fast)")
    p.add_argument("--device", default="cuda", help="torch device, e.g. cuda, cuda:0, cpu")
    p.add_argument("--max-memory-gb", type=float, default=4.0, help="GPU distance-block memory budget")
    p.add_argument("--json", default=None, help="optional path to dump results as JSON")
    args = p.parse_args()

    z_col = "z" if args.dims == 3 else None
    cuda_ok = torch.cuda.is_available()
    print(f"torch {torch.__version__}  device={args.device}  cuda_available={cuda_ok}")
    if args.device.startswith("cuda") and not cuda_ok:
        print("WARNING: CUDA not available; GPU timings will fall back to CPU and be meaningless.")
    print(f"params: k={args.k} dims={args.dims} gpu_dtype={args.dtype} max_n={args.max_n}\n")

    gpu_kw = dict(z_col=z_col, backend="gpu", device=args.device,
                  gpu_dtype=args.dtype, gpu_max_memory_gb=args.max_memory_gb)
    cpu_kw = dict(z_col=z_col, backend="cpu")

    # ---- [1] EQUIVALENCE -------------------------------------------------
    df = make_df(min(20_000, args.max_n), args.k, args.dims, seed=12345)
    cpu = sf(df, **cpu_kw)
    d64 = max_abs_diff(cpu, sf(df, z_col=z_col, backend="gpu", device=args.device,
                               gpu_dtype="float64", gpu_max_memory_gb=args.max_memory_gb))
    d32 = max_abs_diff(cpu, sf(df, z_col=z_col, backend="gpu", device=args.device,
                               gpu_dtype="float32", gpu_max_memory_gb=args.max_memory_gb))
    ok = d64 < 1e-9
    print("[1] EQUIVALENCE (n=%d)" % len(df))
    print(f"    float64 max|delta| = {d64:.3e}   -> {'OK (GPU == CPU)' if ok else 'MISMATCH'}")
    print(f"    float32 max|delta| = {d32:.3e}   (fast mode)\n")

    # ---- [2] BENCHMARK ---------------------------------------------------
    sizes = [n for n in (1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000,
                         200_000, 500_000, 1_000_000) if n <= args.max_n]
    print("[2] BENCHMARK (seconds/call; speedup = CPU/GPU)")
    print(f"    {'N':>9} {'CPU(s)':>10} {'GPU(s)':>10} {'speedup':>9}")
    rows = []
    for n in sizes:
        df = make_df(n, args.k, args.dims, seed=n)
        reps = max(1, int(200_000 // n))
        t_cpu = time_call(lambda: sf(df, **cpu_kw), reps, "cpu")
        t_gpu = time_call(lambda: sf(df, **gpu_kw), reps, args.device)
        speedup = t_cpu / t_gpu if t_gpu > 0 else float("nan")
        print(f"    {n:>9} {t_cpu:>10.4f} {t_gpu:>10.4f} {speedup:>8.2f}x")
        rows.append(dict(n=n, cpu_s=t_cpu, gpu_s=t_gpu, speedup=speedup, reps=reps))

    if args.json:
        payload = dict(
            torch=torch.__version__, device=args.device, cuda_available=cuda_ok,
            k=args.k, dims=args.dims, gpu_dtype=args.dtype,
            equivalence=dict(float64_max_abs_diff=d64, float32_max_abs_diff=d32, exact=ok),
            benchmark=rows,
        )
        pathlib.Path(args.json).write_text(json.dumps(payload, indent=2))
        print(f"\nWrote {args.json}")

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
