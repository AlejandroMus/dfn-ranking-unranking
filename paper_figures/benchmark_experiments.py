#!/usr/bin/env python3
"""Reproducible scaling, order-comparison, and enumeration experiments.

Timing and tracemalloc measurements are deliberately separate.  Proposed-method
query timings exclude cached interval-order preprocessing; cold preprocessing is
reported as its own operation.  Every query checks rank(unrank(i)) == i.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import itertools
import json
import math
import os
import platform
import random
import statistics
import subprocess
import sys
import time
import tracemalloc
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

import paper_algorithm as paper  # noqa: E402


ORDERS = ("t-inc", "lex1", "lex2", "xy")
PREPROCESS_TIMING_REPEATS = 10
RAW_FIELDS = (
    "suite",
    "method",
    "order",
    "n",
    "m",
    "seed",
    "trial",
    "measurement",
    "operation",
    "time_ms",
    "peak_kib",
    "correct",
    "total_dfns",
    "python_version",
    "platform",
    "commit_sha",
    "timestamp_utc",
)


def parse_int_list(value: str) -> list[int]:
    values = [int(item.strip()) for item in value.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected a comma-separated integer list")
    return values


def parse_str_list(value: str) -> list[str]:
    values = [item.strip().lower() for item in value.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected a comma-separated list")
    return values


def parse_configs(value: str) -> list[tuple[int, int]]:
    try:
        configs = [tuple(map(int, item.split(":"))) for item in value.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("configs must look like 5:6,7:6") from exc
    if not configs or any(len(config) != 2 for config in configs):
        raise argparse.ArgumentTypeError("configs must look like 5:6,7:6")
    return [(n, m) for n, m in configs]


def commit_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def working_tree_state() -> tuple[bool, str]:
    """Record whether benchmark code differs from the named commit."""
    try:
        status = subprocess.check_output(
            ["git", "status", "--porcelain"], text=True, stderr=subprocess.DEVNULL
        )
        diff = subprocess.check_output(
            ["git", "diff", "--binary", "HEAD"], stderr=subprocess.DEVNULL
        )
        return bool(status.strip()), hashlib.sha256(diff).hexdigest()
    except (OSError, subprocess.CalledProcessError):
        return False, "unknown"


def benchmark_code_sha256() -> str:
    digest = hashlib.sha256()
    for name in ("paper_algorithm.py", "dfn_cuts_rank_unrank.py", Path(__file__).name):
        digest.update(name.encode("utf-8"))
        digest.update(Path(name).read_bytes())
    return digest.hexdigest()


def metadata() -> dict[str, str]:
    dirty, diff_sha256 = working_tree_state()
    return {
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "commit_sha": commit_sha(),
        "working_tree_dirty": str(dirty).lower(),
        "tracked_diff_sha256": diff_sha256,
        "benchmark_code_sha256": benchmark_code_sha256(),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "arithmetic_model": (
            "Python exact integers (math.comb); empirical timings include bigint "
            "cost, while the manuscript asymptotic bound assumes unit-cost arithmetic."
        ),
        "memory_method": "tracemalloc peak allocated memory, measured separately from timing",
        "timing_method": "time.perf_counter_ns; amortized queries exclude preprocessing",
        "preprocess_timing_repeats": PREPROCESS_TIMING_REPEATS,
        "preprocess_timing_statistic": "arithmetic mean",
    }


class RawWriter:
    def __init__(self, path: Path, common: dict[str, object]):
        self.handle = path.open("w", newline="", encoding="utf-8")
        self.writer = csv.DictWriter(self.handle, fieldnames=RAW_FIELDS)
        self.writer.writeheader()
        self.common = common

    def write(self, **values: object) -> None:
        row = {name: "" for name in RAW_FIELDS}
        row.update(self.common)
        row.update(values)
        self.writer.writerow(row)
        self.handle.flush()

    def close(self) -> None:
        self.handle.close()


def _interval_key(interval: tuple[int, int], order: str) -> tuple[int, int]:
    a, b = interval
    if order == "t-inc":
        return a, -b
    if order == "lex1":
        return a, b
    if order == "lex2":
        return b, a
    if order == "xy":
        return a + b, b - a
    raise ValueError(f"unsupported order: {order}")


def _dfn_key(
    levels: Sequence[int],
    n: int,
    m: int,
    order: str,
) -> tuple:
    """Literal descending-cut key from the manuscript definition."""
    key: list[tuple[int, int]] = []
    for threshold in range(m - 1, 0, -1):
        left = 0
        while left <= n and levels[left] < threshold:
            left += 1
        right = n
        while right >= 0 and levels[right] < threshold:
            right -= 1
        if left > right:
            raise ValueError(f"empty alpha-cut at threshold={threshold}: {levels}")
        key.append(_interval_key((left, right), order))
    return tuple(key)


def enumerate_and_sort_dfns(n: int, m: int, order: str) -> list[tuple[int, ...]]:
    """Complete enumeration baseline: generate every DFN, then sort by cuts."""
    if m < 2:
        raise ValueError("m must be at least 2")
    levels_below_core = range(m - 1)
    dfns: list[tuple[int, ...]] = []
    for core_left in range(n + 1):
        for core_right in range(core_left, n + 1):
            left_values = itertools.combinations_with_replacement(
                levels_below_core, core_left
            )
            right_length = n - core_right
            for left in left_values:
                for right_ascending in itertools.combinations_with_replacement(
                    levels_below_core, right_length
                ):
                    right = tuple(reversed(right_ascending))
                    core = (m - 1,) * (core_right - core_left + 1)
                    dfns.append(tuple(left) + core + right)
    dfns.sort(key=lambda levels: _dfn_key(levels, n, m, order))
    return dfns


def _proposed_preprocess(order: str, n: int, m: int) -> None:
    paper.preprocess_interval_order(order, n, m)


def _proposed_unrank(order: str, n: int, m: int, index: int) -> list[int]:
    if order in ("t-inc", "tinc", "t_inc"):
        _, cuts, _ = paper._unrank_tinc_fast(n, m, index)
        return paper.levels_from_alpha_intervals(cuts, n, m)
    sequence = paper._unrank_engine(order, n, m, index)
    return list(paper.seq_levels(sequence, m))


def _proposed_rank(order: str, n: int, m: int, levels: Sequence[int]) -> int:
    return paper._rank_engine(order, n, m, list(levels))


def _measure_call(function: Callable[[], object]) -> tuple[object, float]:
    start = time.perf_counter_ns()
    result = function()
    elapsed_ms = (time.perf_counter_ns() - start) / 1_000_000.0
    return result, elapsed_ms


def _measure_peak(function: Callable[[], object]) -> tuple[object, float]:
    gc.collect()
    tracemalloc.start()
    try:
        result = function()
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return result, peak / 1024.0


def _prepare_method(
    method: str, order: str, n: int, m: int
) -> tuple[object | None, float, float]:
    """Return reusable state, clean cold-start time, and peak Python memory."""
    if method == "proposed":
        elapsed_samples = []
        for _ in range(PREPROCESS_TIMING_REPEATS):
            paper.clear_preprocessing_caches()
            _, elapsed_ms = _measure_call(lambda: _proposed_preprocess(order, n, m))
            elapsed_samples.append(elapsed_ms)
        paper.clear_preprocessing_caches()
        _, peak_kib = _measure_peak(lambda: _proposed_preprocess(order, n, m))
        return None, statistics.mean(elapsed_samples), peak_kib

    if method == "bruteforce":
        elapsed_samples = []
        for _ in range(PREPROCESS_TIMING_REPEATS):
            timed_state, elapsed_ms = _measure_call(
                lambda: enumerate_and_sort_dfns(n, m, order)
            )
            elapsed_samples.append(elapsed_ms)
            del timed_state
            gc.collect()
        state, peak_kib = _measure_peak(lambda: enumerate_and_sort_dfns(n, m, order))
        return state, statistics.mean(elapsed_samples), peak_kib

    raise ValueError(f"unsupported method: {method}")


def _method_operations(
    method: str, state: object | None, order: str, n: int, m: int
) -> tuple[Callable[[int], list[int]], Callable[[Sequence[int]], int]]:
    if method == "proposed":
        return (
            lambda index: _proposed_unrank(order, n, m, index),
            lambda levels: _proposed_rank(order, n, m, levels),
        )
    if method == "bruteforce":
        ordered = state
        assert isinstance(ordered, list)
        return (
            lambda index: list(ordered[index]),
            lambda levels: ordered.index(tuple(levels)),
        )
    raise ValueError(f"unsupported method: {method}")


def run_configuration(
    writer: RawWriter,
    *,
    method: str,
    order: str,
    n: int,
    m: int,
    trials: int,
    warmups: int,
    measure_memory: bool,
    memory_trials: int,
    seed: int,
    max_enumerated: int,
) -> None:
    total = paper.dcru.total_dfns(n, m)
    label = f"method={method} order={order} n={n} m={m} total={total}"
    print(f"[START] {label}", flush=True)

    if method == "bruteforce" and total > max_enumerated:
        writer.write(
            method=method,
            order=order,
            n=n,
            m=m,
            seed=seed,
            trial=-1,
            measurement="skipped",
            operation="preprocess",
            correct=False,
            total_dfns=total,
        )
        print(f"[SKIP] {label}: exceeds --max-enumerated={max_enumerated}", flush=True)
        return

    state, preprocess_ms, preprocess_peak_kib = _prepare_method(method, order, n, m)
    writer.write(
        method=method,
        order=order,
        n=n,
        m=m,
        seed=seed,
        trial=-1,
        measurement="timing",
        operation="preprocess",
        time_ms=f"{preprocess_ms:.9f}",
        correct=True,
        total_dfns=total,
    )
    writer.write(
        method=method,
        order=order,
        n=n,
        m=m,
        seed=seed,
        trial=-1,
        measurement="memory",
        operation="preprocess",
        peak_kib=f"{preprocess_peak_kib:.6f}",
        correct=True,
        total_dfns=total,
    )

    unrank, rank = _method_operations(method, state, order, n, m)
    rng = random.Random(seed)
    warmup_indices = [rng.randrange(total) for _ in range(warmups)]
    trial_indices = [rng.randrange(total) for _ in range(trials)]

    for index in warmup_indices:
        levels = unrank(index)
        observed = rank(levels)
        if observed != index:
            raise RuntimeError(f"warm-up round trip failed: {label} {index=} {observed=}")

    for trial, index in enumerate(trial_indices):
        levels, unrank_ms = _measure_call(lambda index=index: unrank(index))
        observed, rank_ms = _measure_call(lambda levels=levels: rank(levels))
        correct = observed == index
        if not correct:
            raise RuntimeError(f"timed round trip failed: {label} {index=} {observed=}")
        for operation, elapsed_ms in (("unrank", unrank_ms), ("rank", rank_ms)):
            writer.write(
                method=method,
                order=order,
                n=n,
                m=m,
                seed=seed,
                trial=trial,
                measurement="timing",
                operation=operation,
                time_ms=f"{elapsed_ms:.9f}",
                correct=correct,
                total_dfns=total,
            )

    if measure_memory:
        for trial, index in enumerate(trial_indices[:memory_trials]):
            levels, unrank_peak = _measure_peak(lambda index=index: unrank(index))
            observed, rank_peak = _measure_peak(lambda levels=levels: rank(levels))
            correct = observed == index
            if not correct:
                raise RuntimeError(f"memory round trip failed: {label} {index=} {observed=}")
            for operation, peak_kib in (("unrank", unrank_peak), ("rank", rank_peak)):
                writer.write(
                    method=method,
                    order=order,
                    n=n,
                    m=m,
                    seed=seed,
                    trial=trial,
                    measurement="memory",
                    operation=operation,
                    peak_kib=f"{peak_kib:.6f}",
                    correct=correct,
                    total_dfns=total,
                )

    print(f"[DONE] {label}", flush=True)


def _iqr(values: pd.Series) -> float:
    return float(values.quantile(0.75) - values.quantile(0.25))


def build_summary(raw_path: Path, outdir: Path) -> pd.DataFrame:
    raw = pd.read_csv(raw_path)
    timing = raw[(raw["measurement"] == "timing") & raw["time_ms"].notna()].copy()
    timing["time_ms"] = pd.to_numeric(timing["time_ms"])
    keys = ["method", "order", "n", "m", "operation"]
    timing_summary = (
        timing.groupby(keys, dropna=False)["time_ms"]
        .agg(
            trials="count",
            mean_ms="mean",
            median_ms="median",
            std_ms=lambda values: statistics.stdev(values) if len(values) > 1 else 0.0,
            iqr_ms=_iqr,
        )
        .reset_index()
    )

    memory = raw[(raw["measurement"] == "memory") & raw["peak_kib"].notna()].copy()
    memory["peak_kib"] = pd.to_numeric(memory["peak_kib"])
    memory_summary = (
        memory.groupby(keys, dropna=False)["peak_kib"]
        .agg(memory_trials="count", median_peak_kib="median", max_peak_kib="max")
        .reset_index()
    )
    summary = timing_summary.merge(memory_summary, on=keys, how="outer")
    summary.sort_values(keys, inplace=True)
    summary.to_csv(outdir / "summary.csv", index=False, float_format="%.9f")
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="In future versions `DataFrame.to_latex`.*",
            category=FutureWarning,
        )
        latex = summary.to_latex(index=False, float_format="%.3f")
    (outdir / "summary_table.tex").write_text(latex, encoding="utf-8")
    return summary


def _query_summary(summary: pd.DataFrame) -> pd.DataFrame:
    return summary[summary["operation"].isin(("unrank", "rank"))].copy()


def make_plots(suite: str, summary: pd.DataFrame, outdir: Path) -> None:
    query = _query_summary(summary)
    if query.empty:
        return

    if suite in ("m-scaling", "n-scaling"):
        x_name = "m" if suite == "m-scaling" else "n"
        figure, axes = plt.subplots(1, 2 if suite == "n-scaling" else 1, figsize=(11, 4.5))
        axes_list = list(axes) if hasattr(axes, "__len__") else [axes]
        raw_axis = axes_list[0]
        for operation, part in query.groupby("operation"):
            part = part.sort_values(x_name)
            raw_axis.plot(part[x_name], part["median_ms"], marker="o", label=operation)
        raw_axis.set_xlabel(x_name)
        raw_axis.set_ylabel("Median query time (ms)")
        raw_axis.set_title(f"{suite}: raw timing")
        raw_axis.grid(True, alpha=0.3)
        raw_axis.legend()

        if suite == "n-scaling":
            norm_axis = axes_list[1]
            for operation, part in query.groupby("operation"):
                part = part.sort_values("n").copy()
                scale = part["n"] ** 2 * part["m"] * part["n"].map(math.log)
                part["normalized"] = part["median_ms"] / scale
                norm_axis.plot(part["n"], part["normalized"], marker="o", label=operation)
            norm_axis.set_xlabel("n")
            norm_axis.set_ylabel(r"Median ms / $n^2m\log(n)$")
            norm_axis.set_title("Implementation-bound normalization")
            norm_axis.grid(True, alpha=0.3)
            norm_axis.legend()
        figure.tight_layout()
        figure.savefig(outdir / f"figure_{suite.replace('-', '_')}.png", dpi=180)
        plt.close(figure)
        return

    if suite == "orders":
        pivot = query.pivot_table(
            index=["n", "m", "operation"], columns="order", values="median_ms"
        )
        axis = pivot.plot(kind="bar", figsize=(12, 5))
        axis.set_ylabel("Median query time (ms)")
        axis.set_title("Four interval total orders")
        axis.grid(True, axis="y", alpha=0.3)
        plt.tight_layout()
        plt.savefig(outdir / "figure_order_comparison.png", dpi=180)
        plt.close()
        return

    if suite == "baseline":
        config_keys = ["order", "n", "m"]
        query_totals = (
            query.groupby(["method"] + config_keys, as_index=False)["median_ms"]
            .sum()
            .rename(columns={"median_ms": "roundtrip_query_ms"})
        )
        preprocess = summary[summary["operation"] == "preprocess"][
            ["method"] + config_keys + ["median_ms", "max_peak_kib"]
        ].rename(
            columns={
                "median_ms": "preprocess_ms",
                "max_peak_kib": "preprocess_peak_kib",
            }
        )
        totals = query_totals.merge(preprocess, on=["method"] + config_keys)
        totals["cold_roundtrip_ms"] = (
            totals["preprocess_ms"] + totals["roundtrip_query_ms"]
        )
        proposed = totals[totals["method"] == "proposed"]
        brute = totals[totals["method"] == "bruteforce"]
        merged = brute.merge(proposed, on=config_keys, suffixes=("_brute", "_proposed"))
        if merged.empty:
            return
        merged["speedup"] = (
            merged["cold_roundtrip_ms_brute"]
            / merged["cold_roundtrip_ms_proposed"]
        )
        merged["memory_ratio"] = (
            merged["preprocess_peak_kib_brute"]
            / merged["preprocess_peak_kib_proposed"]
        )
        merged.to_csv(outdir / "baseline_ratios.csv", index=False, float_format="%.6f")
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="In future versions `DataFrame.to_latex`.*",
                category=FutureWarning,
            )
            latex = merged.to_latex(index=False, float_format="%.3f")
        (outdir / "baseline_ratios_table.tex").write_text(latex, encoding="utf-8")
        figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        merged = merged.sort_values("n")
        axes[0].plot(merged["n"], merged["speedup"], marker="o")
        axes[1].plot(merged["n"], merged["memory_ratio"], marker="o")
        axes[0].set_title("Cold brute-force / proposed runtime")
        axes[0].set_ylabel("Speedup (preprocess + median round trip)")
        axes[1].set_title("Preprocessing peak-memory ratio")
        axes[1].set_ylabel("Memory ratio")
        for axis in axes:
            axis.set_xlabel("n")
            axis.grid(True, alpha=0.3)
        figure.tight_layout()
        figure.savefig(outdir / "figure_baseline_speedup_memory.png", dpi=180)
        plt.close(figure)


def configurations(args: argparse.Namespace) -> Iterable[tuple[str, str, int, int]]:
    if args.command == "m-scaling":
        for order in args.orders:
            for m in args.m_values:
                yield "proposed", order, args.n, m
    elif args.command == "n-scaling":
        for order in args.orders:
            for n in args.n_values:
                yield "proposed", order, n, args.m
    elif args.command == "orders":
        for order in args.orders:
            for n in args.n_values:
                for m in args.m_values:
                    yield "proposed", order, n, m
    elif args.command == "baseline":
        for method in args.methods:
            for order in args.orders:
                for n, m in args.configs:
                    yield method, order, n, m


def add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--orders", type=parse_str_list, default=["t-inc"])
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--warmups", type=int, default=10)
    parser.add_argument("--measure-memory", action="store_true")
    parser.add_argument("--memory-trials", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260803)
    parser.add_argument("--max-enumerated", type=int, default=200_000)
    parser.add_argument("--outdir", type=Path, required=True)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)

    m_scaling = subparsers.add_parser("m-scaling")
    add_common(m_scaling)
    m_scaling.add_argument("--n", type=int, default=10)
    m_scaling.add_argument("--m-values", type=parse_int_list, required=True)

    n_scaling = subparsers.add_parser("n-scaling")
    add_common(n_scaling)
    n_scaling.add_argument("--m", type=int, required=True)
    n_scaling.add_argument("--n-values", type=parse_int_list, required=True)

    orders = subparsers.add_parser("orders")
    add_common(orders)
    orders.add_argument("--n-values", type=parse_int_list, required=True)
    orders.add_argument("--m-values", type=parse_int_list, required=True)

    baseline = subparsers.add_parser("baseline")
    add_common(baseline)
    baseline.add_argument("--configs", type=parse_configs, required=True)
    baseline.add_argument(
        "--methods", type=parse_str_list, default=["proposed", "bruteforce"]
    )
    return result


def validate_args(args: argparse.Namespace) -> None:
    invalid_orders = sorted(set(args.orders) - set(ORDERS))
    if invalid_orders:
        raise ValueError(f"unsupported orders: {invalid_orders}")
    if args.trials < 1 or args.warmups < 0 or args.memory_trials < 1:
        raise ValueError("trials/memory-trials must be positive and warmups non-negative")
    if hasattr(args, "methods"):
        invalid_methods = sorted(set(args.methods) - {"proposed", "bruteforce"})
        if invalid_methods:
            raise ValueError(f"unsupported methods: {invalid_methods}")


def main() -> None:
    args = parser().parse_args()
    validate_args(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    meta = metadata()
    (args.outdir / "environment.json").write_text(
        json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    serializable_args = {
        key: str(value) if isinstance(value, Path) else value
        for key, value in vars(args).items()
    }
    (args.outdir / "configuration.json").write_text(
        json.dumps(serializable_args, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    raw_path = args.outdir / "raw_trials.csv"
    writer = RawWriter(
        raw_path,
        {
            "suite": args.command,
            "python_version": meta["python_version"],
            "platform": meta["platform"],
            "commit_sha": meta["commit_sha"],
            "timestamp_utc": meta["timestamp_utc"],
        },
    )
    try:
        for config_index, (method, order, n, m) in enumerate(configurations(args)):
            run_configuration(
                writer,
                method=method,
                order=order,
                n=n,
                m=m,
                trials=args.trials,
                warmups=args.warmups,
                measure_memory=args.measure_memory,
                memory_trials=min(args.memory_trials, args.trials),
                seed=args.seed + config_index * 100_003,
                max_enumerated=args.max_enumerated,
            )
    finally:
        writer.close()

    summary = build_summary(raw_path, args.outdir)
    make_plots(args.command, summary, args.outdir)
    print(f"[COMPLETE] suite={args.command} outdir={args.outdir}", flush=True)


if __name__ == "__main__":
    main()
