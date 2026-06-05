#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
import zarr
from sklearn.cluster import MiniBatchKMeans
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


DEFAULT_TRANSCRIPTS_ZARR = Path(
    "/data/taobo.hu/pyxenium_lr_benchmark_2026-04/data/source_cache/breast/"
    "WTA_Preview_FFPE_Breast_Cancer_outs/spatialdata.zarr/points/transcripts"
)
DEFAULT_CELLFREE_COSTE_DIR = Path(
    "/data/taobo.hu/atera_breast_cellfree_transcript_pipeline/full_transcript_grid96_qv20_select1024_geneonly"
)
DEFAULT_OUTPUT_DIR = Path("/data/taobo.hu/atera_breast_transcript_point_histoseg/coste1024_point_knn_domains12")


@dataclass
class TimedResult:
    step: str
    seconds: float
    peak_rss_gb: float
    payload: Any


class PeakMemoryTracker:
    def __init__(self, interval_seconds: float = 0.5) -> None:
        self.interval_seconds = interval_seconds
        self.process = psutil.Process()
        self.peak_rss = 0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def __enter__(self) -> "PeakMemoryTracker":
        self._thread = threading.Thread(target=self._sample, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2)
        self._sample_once()

    def _sample(self) -> None:
        while not self._stop.is_set():
            self._sample_once()
            self._stop.wait(self.interval_seconds)

    def _sample_once(self) -> None:
        rss = 0
        try:
            rss += self.process.memory_info().rss
            for child in self.process.children(recursive=True):
                try:
                    rss += child.memory_info().rss
                except psutil.Error:
                    continue
        except psutil.Error:
            return
        self.peak_rss = max(self.peak_rss, int(rss))


class PointReservoir:
    def __init__(self, max_points: int, seed: int, include_domain: bool = False) -> None:
        self.max_points = int(max_points)
        self.include_domain = include_domain
        self.rng = np.random.default_rng(seed)
        self.keys = np.empty(0, dtype=np.float32)
        self.x = np.empty(0, dtype=np.float32)
        self.y = np.empty(0, dtype=np.float32)
        self.gene_identity = np.empty(0, dtype=np.int32)
        self.module_id = np.empty(0, dtype=np.int16)
        self.domain_id = np.empty(0, dtype=np.int16)

    def add(
        self,
        x: np.ndarray,
        y: np.ndarray,
        gene_identity: np.ndarray,
        module_id: np.ndarray,
        domain_id: np.ndarray | None = None,
    ) -> None:
        if self.max_points <= 0 or len(x) == 0:
            return
        keys = self.rng.random(len(x), dtype=np.float32)
        self.keys = np.concatenate([self.keys, keys])
        self.x = np.concatenate([self.x, x.astype(np.float32, copy=False)])
        self.y = np.concatenate([self.y, y.astype(np.float32, copy=False)])
        self.gene_identity = np.concatenate([self.gene_identity, gene_identity.astype(np.int32, copy=False)])
        self.module_id = np.concatenate([self.module_id, module_id.astype(np.int16, copy=False)])
        if self.include_domain:
            if domain_id is None:
                raise ValueError("domain_id is required for a domain reservoir")
            self.domain_id = np.concatenate([self.domain_id, domain_id.astype(np.int16, copy=False)])
        if len(self.keys) > self.max_points:
            keep = np.argpartition(self.keys, self.max_points - 1)[: self.max_points]
            self._take(keep)

    def _take(self, keep: np.ndarray) -> None:
        order = np.argsort(self.keys[keep])
        keep = keep[order]
        self.keys = self.keys[keep]
        self.x = self.x[keep]
        self.y = self.y[keep]
        self.gene_identity = self.gene_identity[keep]
        self.module_id = self.module_id[keep]
        if self.include_domain:
            self.domain_id = self.domain_id[keep]

    def to_frame(self) -> pd.DataFrame:
        frame = pd.DataFrame(
            {
                "x": self.x,
                "y": self.y,
                "gene_identity": self.gene_identity,
                "coste_module_id": self.module_id.astype(int),
            }
        )
        if self.include_domain:
            frame["histoseg_structure_id"] = self.domain_id.astype(int)
            frame["histoseg_structure_name"] = [f"TranscriptPoint-HistoSeg-{x:02d}" for x in self.domain_id]
        return frame


def timed(step: str, func, *args, **kwargs) -> TimedResult:
    t0 = time.perf_counter()
    with PeakMemoryTracker() as tracker:
        payload = func(*args, **kwargs)
    return TimedResult(step=step, seconds=time.perf_counter() - t0, peak_rss_gb=tracker.peak_rss / (1024**3), payload=payload)


def _log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def transcript_range(n_rows: int, chunk_rows: int, max_transcripts: int) -> list[tuple[int, int]]:
    limit = n_rows if max_transcripts <= 0 else min(n_rows, max_transcripts)
    return [(start, min(start + chunk_rows, limit)) for start in range(0, limit, chunk_rows)]


def load_module_map(gene_summary_path: Path, modules_path: Path) -> dict[str, Any]:
    gene_summary = pd.read_csv(gene_summary_path)
    modules = pd.read_csv(modules_path)
    if "coste_module_id" not in modules.columns:
        raise ValueError(f"{modules_path} does not contain coste_module_id")
    table = modules[["gene", "coste_module_id", "coste_module_name", "importance_score"]].merge(
        gene_summary[["gene_identity", "gene", "transcript_count"]],
        on="gene",
        how="left",
    )
    missing = table.loc[table["gene_identity"].isna(), "gene"].astype(str).tolist()
    if missing:
        raise ValueError(f"{len(missing)} module genes have no gene_identity, e.g. {missing[:5]}")
    table["gene_identity"] = table["gene_identity"].astype(np.int32)
    table["coste_module_id"] = table["coste_module_id"].astype(np.int16)
    max_identity = int(table["gene_identity"].max())
    identity_to_module = np.full(max_identity + 1, -1, dtype=np.int16)
    identity_to_module[table["gene_identity"].to_numpy(dtype=np.int32)] = table["coste_module_id"].to_numpy(dtype=np.int16)
    return {
        "modules": modules,
        "module_gene_table": table,
        "identity_to_module": identity_to_module,
        "selected_identity_count": int(table["gene_identity"].nunique()),
        "selected_transcript_count_from_summary": int(table["transcript_count"].sum()),
    }


def collect_training_points(
    transcripts_zarr: Path,
    *,
    identity_to_module: np.ndarray,
    qv_min: float,
    chunk_rows: int,
    max_transcripts: int,
    max_train_points: int,
    random_state: int,
) -> dict[str, Any]:
    group = zarr.open_group(str(transcripts_zarr), mode="r")
    n_rows = int(group.attrs.get("n_rows", group["x"].shape[0]))
    chunks = transcript_range(n_rows, chunk_rows, max_transcripts)
    reservoir = PointReservoir(max_train_points, seed=random_state)
    selected_seen = 0
    kept_seen = 0
    module_counts: dict[int, int] = {}
    for chunk_idx, (start, stop) in enumerate(chunks, start=1):
        valid = np.asarray(group["valid"][start:stop], dtype=bool)
        qv = np.asarray(group["quality_score"][start:stop], dtype=np.float32)
        mask = valid & (qv >= qv_min)
        kept_seen += int(mask.sum())
        if not mask.any():
            continue
        raw_ids = np.asarray(group["gene_identity"][start:stop], dtype=np.int64)[mask]
        in_map = raw_ids < len(identity_to_module)
        raw_ids = raw_ids[in_map]
        module_ids = identity_to_module[raw_ids]
        selected_mask = module_ids > 0
        if not selected_mask.any():
            continue
        raw_ids = raw_ids[selected_mask].astype(np.int32, copy=False)
        module_ids = module_ids[selected_mask].astype(np.int16, copy=False)
        selected_seen += int(len(raw_ids))
        unique_modules, counts = np.unique(module_ids, return_counts=True)
        for module_id, count in zip(unique_modules, counts):
            module_counts[int(module_id)] = module_counts.get(int(module_id), 0) + int(count)
        x = np.asarray(group["x"][start:stop], dtype=np.float32)[mask][in_map][selected_mask]
        y = np.asarray(group["y"][start:stop], dtype=np.float32)[mask][in_map][selected_mask]
        reservoir.add(x, y, raw_ids, module_ids)
        if chunk_idx == 1 or chunk_idx % 25 == 0 or chunk_idx == len(chunks):
            _log(
                "point sample chunk "
                f"{chunk_idx:,}/{len(chunks):,}: selected {selected_seen:,}, train reservoir {len(reservoir.x):,}"
            )
    return {
        "n_rows": n_rows,
        "n_rows_scanned": chunks[-1][1] if chunks else 0,
        "n_valid_qv_transcripts_seen": int(kept_seen),
        "n_selected_module_transcripts_seen": int(selected_seen),
        "module_counts": module_counts,
        "training_points": reservoir.to_frame(),
    }


def build_point_features(
    coords: np.ndarray,
    module_ids: np.ndarray,
    *,
    module_values: list[int],
    neighborhood_k: int,
    feature_batch_size: int,
    spatial_weight: float,
) -> dict[str, Any]:
    module_to_col = {int(module_id): idx for idx, module_id in enumerate(module_values)}
    module_codes = np.array([module_to_col[int(module_id)] for module_id in module_ids], dtype=np.int16)
    k = int(min(max(1, neighborhood_k), len(coords)))
    _log(f"fitting point-neighborhood index for {len(coords):,} transcript points, k={k}")
    nbrs = NearestNeighbors(n_neighbors=k, algorithm="kd_tree", leaf_size=40, n_jobs=-1)
    nbrs.fit(coords)
    composition = np.zeros((len(coords), len(module_values)), dtype=np.float32)
    density = np.zeros((len(coords), 1), dtype=np.float32)
    for start in range(0, len(coords), feature_batch_size):
        stop = min(start + feature_batch_size, len(coords))
        distances, indices = nbrs.kneighbors(coords[start:stop], return_distance=True)
        neighbor_modules = module_codes[indices]
        hist = np.zeros((stop - start, len(module_values)), dtype=np.float32)
        np.add.at(hist, (np.repeat(np.arange(stop - start), k), neighbor_modules.ravel()), 1.0)
        hist /= float(k)
        radius = np.maximum(distances[:, -1], 1e-3)
        density[start:stop, 0] = np.log1p(k / (math.pi * radius * radius))
        composition[start:stop] = hist
        if start == 0 or stop == len(coords) or (start // feature_batch_size) % 10 == 0:
            _log(f"built local COSTE-neighborhood features {stop:,}/{len(coords):,}")
    expression_features = np.hstack([composition, density]).astype(np.float32)
    expression_scaler = StandardScaler()
    expression_scaled = expression_scaler.fit_transform(expression_features).astype(np.float32)
    spatial_scaler = StandardScaler()
    xy_scaled = spatial_scaler.fit_transform(coords).astype(np.float32) * float(spatial_weight)
    features = np.hstack([expression_scaled, xy_scaled]).astype(np.float32)
    return {
        "features": features,
        "neighbor_index": nbrs,
        "expression_scaler": expression_scaler,
        "spatial_scaler": spatial_scaler,
        "module_values": module_values,
        "module_codes": module_codes,
        "neighborhood_k": k,
    }


def fit_point_domains(
    training_points: pd.DataFrame,
    *,
    n_domains: int,
    neighborhood_k: int,
    feature_batch_size: int,
    spatial_weight: float,
    random_state: int,
) -> dict[str, Any]:
    coords = training_points[["x", "y"]].to_numpy(dtype=np.float32)
    module_ids = training_points["coste_module_id"].to_numpy(dtype=np.int16)
    module_values = sorted(int(x) for x in np.unique(module_ids))
    feature_model = build_point_features(
        coords,
        module_ids,
        module_values=module_values,
        neighborhood_k=neighborhood_k,
        feature_batch_size=feature_batch_size,
        spatial_weight=spatial_weight,
    )
    _log(f"fitting MiniBatchKMeans HistoSeg domains, n_domains={n_domains}")
    kmeans = MiniBatchKMeans(
        n_clusters=n_domains,
        random_state=random_state,
        batch_size=min(16384, max(1024, len(training_points) // 10)),
        n_init=30,
        reassignment_ratio=0.01,
    )
    labels = kmeans.fit_predict(feature_model["features"]).astype(np.int16) + 1
    out = training_points.copy()
    out["histoseg_structure_id"] = labels.astype(int)
    out["histoseg_structure_name"] = [f"TranscriptPoint-HistoSeg-{x:02d}" for x in labels]
    return {
        "training_domains": out,
        "feature_model": feature_model,
        "kmeans": kmeans,
    }


def assign_selected_points_by_nearest_training_domain(
    transcripts_zarr: Path,
    *,
    identity_to_module: np.ndarray,
    training_domains: pd.DataFrame,
    modules: pd.DataFrame,
    qv_min: float,
    chunk_rows: int,
    max_transcripts: int,
    max_assign_points: int,
    assignment_neighbors: int,
    preview_points: int,
    random_state: int,
    write_domain_chunks: bool,
    output_dir: Path,
) -> dict[str, Any]:
    coords_train = training_domains[["x", "y"]].to_numpy(dtype=np.float32)
    domain_train = training_domains["histoseg_structure_id"].to_numpy(dtype=np.int16)
    module_train = training_domains["coste_module_id"].to_numpy(dtype=np.int16)
    n_domains = int(training_domains["histoseg_structure_id"].max())
    n_modules = int(max(modules["coste_module_id"].max(), module_train.max()))
    neighbors = int(min(max(1, assignment_neighbors), len(coords_train)))
    _log(f"fitting nearest training-domain index for point assignment, k={neighbors}")
    nbrs = NearestNeighbors(n_neighbors=neighbors, algorithm="kd_tree", leaf_size=40, n_jobs=-1)
    nbrs.fit(coords_train)

    group = zarr.open_group(str(transcripts_zarr), mode="r")
    n_rows = int(group.attrs.get("n_rows", group["x"].shape[0]))
    chunks = transcript_range(n_rows, chunk_rows, max_transcripts)
    preview = PointReservoir(preview_points, seed=random_state + 1009, include_domain=True)
    domain_counts = np.zeros(n_domains + 1, dtype=np.int64)
    domain_module_counts = np.zeros((n_domains + 1, n_modules + 1), dtype=np.int64)
    assigned_total = 0
    chunk_output_dir = output_dir / "transcript_point_domain_chunks"
    if write_domain_chunks:
        chunk_output_dir.mkdir(parents=True, exist_ok=True)
    written_parts = 0

    for chunk_idx, (start, stop) in enumerate(chunks, start=1):
        if max_assign_points > 0 and assigned_total >= max_assign_points:
            break
        valid = np.asarray(group["valid"][start:stop], dtype=bool)
        qv = np.asarray(group["quality_score"][start:stop], dtype=np.float32)
        mask = valid & (qv >= qv_min)
        if not mask.any():
            continue
        raw_ids = np.asarray(group["gene_identity"][start:stop], dtype=np.int64)[mask]
        in_map = raw_ids < len(identity_to_module)
        raw_ids = raw_ids[in_map]
        module_ids = identity_to_module[raw_ids]
        selected_mask = module_ids > 0
        if not selected_mask.any():
            continue
        raw_ids = raw_ids[selected_mask].astype(np.int32, copy=False)
        module_ids = module_ids[selected_mask].astype(np.int16, copy=False)
        x = np.asarray(group["x"][start:stop], dtype=np.float32)[mask][in_map][selected_mask]
        y = np.asarray(group["y"][start:stop], dtype=np.float32)[mask][in_map][selected_mask]
        if max_assign_points > 0 and assigned_total + len(x) > max_assign_points:
            keep_n = max_assign_points - assigned_total
            raw_ids = raw_ids[:keep_n]
            module_ids = module_ids[:keep_n]
            x = x[:keep_n]
            y = y[:keep_n]
        coords = np.column_stack([x, y]).astype(np.float32)
        domains = majority_domain_predict(nbrs, coords, domain_train, n_domains)
        domain_counts += np.bincount(domains, minlength=n_domains + 1)
        flat = domains.astype(np.int64) * (n_modules + 1) + module_ids.astype(np.int64)
        binc = np.bincount(flat, minlength=(n_domains + 1) * (n_modules + 1))
        domain_module_counts += binc.reshape(n_domains + 1, n_modules + 1)
        preview.add(x, y, raw_ids, module_ids, domains)
        if write_domain_chunks:
            part = pd.DataFrame(
                {
                    "x": x,
                    "y": y,
                    "gene_identity": raw_ids,
                    "coste_module_id": module_ids.astype(int),
                    "histoseg_structure_id": domains.astype(int),
                    "histoseg_structure_name": [f"TranscriptPoint-HistoSeg-{d:02d}" for d in domains],
                }
            )
            part.to_parquet(chunk_output_dir / f"part_{written_parts:05d}.parquet", index=False)
            written_parts += 1
        assigned_total += int(len(domains))
        if chunk_idx == 1 or chunk_idx % 25 == 0 or chunk_idx == len(chunks):
            _log(f"assigned transcript-point domains chunk {chunk_idx:,}/{len(chunks):,}: {assigned_total:,} points")

    summary = summarize_point_domains(domain_counts, domain_module_counts, modules)
    preview_frame = preview.to_frame()
    return {
        "assigned_points": int(assigned_total),
        "domain_counts": domain_counts,
        "domain_module_counts": domain_module_counts,
        "summary": summary,
        "preview_points": preview_frame,
        "written_parts": int(written_parts),
    }


def majority_domain_predict(
    nbrs: NearestNeighbors,
    coords: np.ndarray,
    domain_train: np.ndarray,
    n_domains: int,
    *,
    batch_size: int = 100_000,
) -> np.ndarray:
    out = np.empty(len(coords), dtype=np.int16)
    for start in range(0, len(coords), batch_size):
        stop = min(start + batch_size, len(coords))
        indices = nbrs.kneighbors(coords[start:stop], return_distance=False)
        neighbor_domains = domain_train[indices]
        hist = np.zeros((stop - start, n_domains + 1), dtype=np.int16)
        np.add.at(hist, (np.repeat(np.arange(stop - start), indices.shape[1]), neighbor_domains.ravel()), 1)
        out[start:stop] = np.argmax(hist[:, 1:], axis=1).astype(np.int16) + 1
    return out


def summarize_point_domains(domain_counts: np.ndarray, domain_module_counts: np.ndarray, modules: pd.DataFrame) -> pd.DataFrame:
    top_gene_lookup = (
        modules.sort_values("importance_score", ascending=False)
        .groupby("coste_module_name")["gene"]
        .apply(lambda x: "/".join(x.head(5).astype(str)))
        .to_dict()
    )
    module_name_lookup = (
        modules.drop_duplicates("coste_module_id")
        .set_index("coste_module_id")["coste_module_name"]
        .astype(str)
        .to_dict()
    )
    total = int(domain_counts.sum())
    rows: list[dict[str, Any]] = []
    for domain_id in range(1, len(domain_counts)):
        count = int(domain_counts[domain_id])
        if count == 0:
            continue
        module_counts = domain_module_counts[domain_id]
        top_module_id = int(np.argmax(module_counts[1:]) + 1)
        top_module_name = module_name_lookup.get(top_module_id, f"CF_COSTE_M{top_module_id:02d}")
        rows.append(
            {
                "histoseg_structure_id": int(domain_id),
                "histoseg_structure_name": f"TranscriptPoint-HistoSeg-{domain_id:02d}",
                "n_transcript_points": count,
                "fraction_transcript_points": float(count / max(total, 1)),
                "top_coste_module": top_module_name,
                "top_coste_module_transcripts": int(module_counts[top_module_id]),
                "top_coste_module_fraction": float(module_counts[top_module_id] / max(count, 1)),
                "top_coste_module_genes": top_gene_lookup.get(top_module_name, ""),
            }
        )
    return pd.DataFrame(rows).sort_values("histoseg_structure_id")


def plot_point_domains(points: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
    scatter = ax.scatter(
        points["x"],
        points["y"],
        c=points["histoseg_structure_id"],
        s=0.35,
        cmap="tab20",
        linewidths=0,
        alpha=0.75,
        rasterized=True,
    )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("Transcript-point HistoSeg domains (no spatial grid)")
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("histoseg domain")
    fig.savefig(output_dir / "transcript_point_histoseg_domains.png", dpi=300)
    fig.savefig(output_dir / "transcript_point_histoseg_domains.svg")
    plt.close(fig)


def timed_row(result: TimedResult, method: str, n_items: int) -> dict[str, Any]:
    return {
        "step": result.step,
        "method": method,
        "n_items": int(n_items),
        "seconds": float(result.seconds),
        "peak_rss_gb": float(result.peak_rss_gb),
    }


def run(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    timed_rows: list[dict[str, Any]] = []

    module_map = load_module_map(args.gene_summary, args.modules)
    module_map["module_gene_table"].to_csv(args.output_dir / "transcript_point_coste_module_gene_map.csv", index=False)

    sample_result = timed(
        "collect_transcript_point_training_sample",
        collect_training_points,
        args.transcripts_zarr,
        identity_to_module=module_map["identity_to_module"],
        qv_min=args.qv_min,
        chunk_rows=args.chunk_rows,
        max_transcripts=args.max_transcripts,
        max_train_points=args.max_train_points,
        random_state=args.random_state,
    )
    sample_payload = sample_result.payload
    training_points = sample_payload["training_points"]
    if len(training_points) == 0:
        raise ValueError("No selected COSTE-module transcript points were sampled.")
    training_points.to_parquet(args.output_dir / "transcript_point_histoseg_training_points.parquet", index=False)
    timed_rows.append(timed_row(sample_result, "transcript_point_histoseg", len(training_points)))

    fit_result = timed(
        "fit_transcript_point_histoseg",
        fit_point_domains,
        training_points,
        n_domains=args.n_domains,
        neighborhood_k=args.neighborhood_k,
        feature_batch_size=args.feature_batch_size,
        spatial_weight=args.spatial_weight,
        random_state=args.random_state,
    )
    fit_payload = fit_result.payload
    training_domains = fit_payload["training_domains"]
    training_domains.to_parquet(args.output_dir / "transcript_point_histoseg_training_domains.parquet", index=False)
    timed_rows.append(timed_row(fit_result, "transcript_point_histoseg", len(training_domains)))

    assign_result = timed(
        "assign_selected_transcript_points",
        assign_selected_points_by_nearest_training_domain,
        args.transcripts_zarr,
        identity_to_module=module_map["identity_to_module"],
        training_domains=training_domains,
        modules=module_map["modules"],
        qv_min=args.qv_min,
        chunk_rows=args.chunk_rows,
        max_transcripts=args.max_transcripts,
        max_assign_points=args.max_assign_points,
        assignment_neighbors=args.assignment_neighbors,
        preview_points=args.preview_points,
        random_state=args.random_state,
        write_domain_chunks=args.write_domain_chunks,
        output_dir=args.output_dir,
    )
    assign_payload = assign_result.payload
    timed_rows.append(timed_row(assign_result, "transcript_point_histoseg", assign_payload["assigned_points"]))

    assign_payload["summary"].to_csv(args.output_dir / "transcript_point_histoseg_domain_summary.csv", index=False)
    preview = assign_payload["preview_points"]
    preview.to_parquet(args.output_dir / "transcript_point_histoseg_preview_points.parquet", index=False)
    preview.to_csv(args.output_dir / "transcript_point_histoseg_preview_points.csv", index=False)
    plot_point_domains(preview, args.output_dir)

    method_summary = pd.DataFrame(timed_rows)
    method_summary.to_csv(args.output_dir / "method_summary.csv", index=False)
    manifest = {
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "workflow": "transcript_point_histoseg_no_spatial_grid",
        "transcripts_zarr": str(args.transcripts_zarr),
        "gene_summary": str(args.gene_summary),
        "modules": str(args.modules),
        "qv_min": float(args.qv_min),
        "max_transcripts": int(args.max_transcripts),
        "uses_spatial_grid": False,
        "histoseg_unit": "transcript_point",
        "selected_coste_genes": int(module_map["selected_identity_count"]),
        "selected_transcript_count_from_summary": int(module_map["selected_transcript_count_from_summary"]),
        "n_rows_scanned_for_sample": int(sample_payload["n_rows_scanned"]),
        "n_valid_qv_transcripts_seen": int(sample_payload["n_valid_qv_transcripts_seen"]),
        "n_selected_module_transcripts_seen": int(sample_payload["n_selected_module_transcripts_seen"]),
        "n_training_points": int(len(training_points)),
        "n_assigned_selected_transcript_points": int(assign_payload["assigned_points"]),
        "n_preview_points": int(len(preview)),
        "n_domains": int(args.n_domains),
        "n_gene_modules": int(module_map["modules"]["coste_module_id"].nunique()),
        "neighborhood_k": int(args.neighborhood_k),
        "assignment_neighbors": int(args.assignment_neighbors),
        "spatial_weight": float(args.spatial_weight),
        "write_domain_chunks": bool(args.write_domain_chunks),
        "written_domain_chunk_parts": int(assign_payload["written_parts"]),
        "outputs": {
            "method_summary": "method_summary.csv",
            "domain_summary": "transcript_point_histoseg_domain_summary.csv",
            "training_domains": "transcript_point_histoseg_training_domains.parquet",
            "preview_points": "transcript_point_histoseg_preview_points.parquet",
            "domain_figure": "transcript_point_histoseg_domains.png",
        },
    }
    _write_json(args.output_dir / "transcript_point_histoseg_manifest.json", manifest)
    print(json.dumps({"manifest": manifest, "method_summary": method_summary.to_dict(orient="records")}, indent=2), flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Transcript-point HistoSeg using COSTE modules without spatial grid bins.")
    parser.add_argument("--transcripts-zarr", type=Path, default=DEFAULT_TRANSCRIPTS_ZARR)
    parser.add_argument("--gene-summary", type=Path, default=DEFAULT_CELLFREE_COSTE_DIR / "cellfree_transcript_gene_summary.csv")
    parser.add_argument("--modules", type=Path, default=DEFAULT_CELLFREE_COSTE_DIR / "cellfree_coste_gene_modules.csv")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--qv-min", type=float, default=20.0)
    parser.add_argument("--chunk-rows", type=int, default=2_000_000)
    parser.add_argument("--max-transcripts", type=int, default=0, help="0 means all transcript rows; positive values are for smoke tests.")
    parser.add_argument("--max-train-points", type=int, default=750_000)
    parser.add_argument("--max-assign-points", type=int, default=0, help="0 assigns all selected COSTE-module transcript points.")
    parser.add_argument("--preview-points", type=int, default=500_000)
    parser.add_argument("--neighborhood-k", type=int, default=128)
    parser.add_argument("--assignment-neighbors", type=int, default=15)
    parser.add_argument("--feature-batch-size", type=int, default=25_000)
    parser.add_argument("--n-domains", type=int, default=12)
    parser.add_argument("--spatial-weight", type=float, default=0.35)
    parser.add_argument("--random-state", type=int, default=7)
    parser.add_argument("--write-domain-chunks", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
