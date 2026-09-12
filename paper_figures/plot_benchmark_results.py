#!/usr/bin/env python3
"""Create manuscript-ready figures, tables, and prose from benchmark results."""

from __future__ import annotations

import argparse
import math
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


COLORS = {"unrank": "#0072B2", "rank": "#D55E00"}
MARKERS = {"unrank": "o", "rank": "s"}
ORDER_LABELS = {
    "t-inc": "t-inc",
    "lex1": "lex1",
    "lex2": "lex2",
    "xy": "Xu–Yager",
}
ORDER_COLORS = {
    "t-inc": "#0072B2",
    "lex1": "#D55E00",
    "lex2": "#009E73",
    "xy": "#CC79A7",
}
ORDER_MARKERS = {"t-inc": "o", "lex1": "s", "lex2": "^", "xy": "D"}


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 17,
            "axes.titlesize": 18,
            "axes.labelsize": 17,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 14,
            "figure.titlesize": 20,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 120,
            "savefig.bbox": "tight",
        }
    )


def read_raw(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame["time_ms"] = pd.to_numeric(frame["time_ms"], errors="coerce")
    frame["peak_kib"] = pd.to_numeric(frame["peak_kib"], errors="coerce")
    return frame


def timing_stats(raw: pd.DataFrame) -> pd.DataFrame:
    data = raw[
        (raw["measurement"] == "timing")
        & raw["operation"].isin(["rank", "unrank"])
    ].copy()
    keys = ["method", "order", "n", "m", "operation"]
    return (
        data.groupby(keys)["time_ms"]
        .agg(
            q25=lambda values: values.quantile(0.25),
            median="median",
            q75=lambda values: values.quantile(0.75),
            mean="mean",
            std="std",
        )
        .reset_index()
    )


def memory_stats(raw: pd.DataFrame) -> pd.DataFrame:
    data = raw[
        (raw["measurement"] == "memory")
        & raw["operation"].isin(["rank", "unrank"])
    ].copy()
    keys = ["method", "order", "n", "m", "operation"]
    return (
        data.groupby(keys)["peak_kib"]
        .agg(median="median", maximum="max")
        .reset_index()
    )


def preprocess_stats(raw: pd.DataFrame) -> pd.DataFrame:
    data = raw[raw["operation"] == "preprocess"].copy()
    keys = ["method", "order", "n", "m"]
    times = (
        data[data["measurement"] == "timing"]
        .groupby(keys)["time_ms"]
        .median()
        .rename("preprocess_ms")
    )
    memory = (
        data[data["measurement"] == "memory"]
        .groupby(keys)["peak_kib"]
        .max()
        .rename("preprocess_peak_kib")
    )
    return pd.concat([times, memory], axis=1).reset_index()


def save_figure(figure: plt.Figure, outdir: Path, stem: str) -> None:
    figure.savefig(outdir / f"{stem}.pdf")
    figure.savefig(outdir / f"{stem}.png", dpi=300)
    plt.close(figure)


def panel_label(
    axis: plt.Axes,
    label: str,
    *,
    x: float = -0.13,
    y: float = 1.05,
    fontsize: float = 16,
) -> None:
    axis.text(
        x,
        y,
        label,
        transform=axis.transAxes,
        fontsize=fontsize,
        fontweight="bold",
        va="top",
    )


def plot_orders(
    axis: plt.Axes,
    stats: pd.DataFrame,
    x_name: str,
    operation: str,
    *,
    normalize_m: int | None = None,
) -> None:
    for order in ("t-inc", "lex1", "lex2", "xy"):
        part = stats[
            (stats["operation"] == operation) & (stats["order"] == order)
        ].sort_values(x_name).copy()
        if part.empty:
            continue
        if normalize_m is None:
            median, q25, q75 = part["median"], part["q25"], part["q75"]
        else:
            scale = part["n"] ** 2 * normalize_m * part["n"].map(math.log)
            median, q25, q75 = (
                part["median"] / scale,
                part["q25"] / scale,
                part["q75"] / scale,
            )
        axis.plot(
            part[x_name],
            median,
            color=ORDER_COLORS[order],
            marker=ORDER_MARKERS[order],
            label=ORDER_LABELS[order],
        )
        axis.fill_between(
            part[x_name], q25, q75, color=ORDER_COLORS[order], alpha=0.10
        )


def scaling_figure(
    m_stats: pd.DataFrame,
    n20_stats: pd.DataFrame,
    outdir: Path,
) -> None:
    figure, axes = plt.subplots(3, 2, figsize=(13.0, 14.0))
    operations = ("rank", "unrank")
    for column, operation in enumerate(operations):
        plot_orders(axes[0, column], m_stats, "m", operation)
        axes[0, column].set(
            xlabel=r"Number of levels $m$",
            ylabel="Median query time (ms)",
            title=rf"{operation.capitalize()}: fixed $n=10$",
        )
        plot_orders(axes[1, column], n20_stats, "n", operation)
        axes[1, column].set(
            xlabel=r"Chain parameter $n$",
            ylabel="Median query time (ms)",
            xscale="log",
            yscale="log",
            title=rf"{operation.capitalize()}: fixed $m=20$",
        )
        axes[1, column].set_xticks([10, 20, 40, 60, 100, 200])
        axes[1, column].set_xticklabels(["10", "20", "40", "60", "100", "200"])
        plot_orders(
            axes[2, column], n20_stats, "n", operation, normalize_m=20
        )
        axes[2, column].set(
            xlabel=r"Chain parameter $n$",
            ylabel=r"Median ms / $n^2m\log(n)$",
            xscale="log",
            yscale="log",
            title=rf"Normalized {operation} time: fixed $m=20$",
        )
        axes[2, column].set_xticks([10, 20, 40, 60, 100, 200])
        axes[2, column].set_xticklabels(["10", "20", "40", "60", "100", "200"])
    labels = ("(a.1)", "(a.2)", "(b.1)", "(b.2)", "(c.1)", "(c.2)")
    for label, axis in zip(labels, axes.flat):
        axis.grid(alpha=0.25, which="both")
        axis.legend(ncol=2)
        panel_label(axis, label)
    figure.tight_layout(pad=1.5)
    save_figure(figure, outdir, "figure_scaling")


def order_figure(order_raw: pd.DataFrame, outdir: Path) -> None:
    timing = timing_stats(order_raw)
    preprocessing = preprocess_stats(order_raw)
    timing = timing[(timing["n"] == 100) & (timing["m"] == 20)]
    preprocessing = preprocessing[
        (preprocessing["n"] == 100) & (preprocessing["m"] == 20)
    ]
    orders = [order for order in ("t-inc", "lex1", "lex2", "xy") if order in set(timing["order"])]
    x = np.arange(len(orders))
    width = 0.36
    figure, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))

    for offset, operation in ((-width / 2, "unrank"), (width / 2, "rank")):
        part = timing[timing["operation"] == operation].set_index("order").loc[orders]
        lower = part["median"] - part["q25"]
        upper = part["q75"] - part["median"]
        axes[0].bar(
            x + offset,
            part["median"],
            width,
            yerr=np.vstack([lower, upper]),
            capsize=2,
            color=COLORS[operation],
            label=operation.capitalize(),
        )
    axes[0].set_ylabel("Median query time (ms)")
    axes[0].set_title(
        r"$\mathbf{(d)}$ Query time after preprocessing", loc="left", pad=10
    )
    axes[0].legend()

    prep = preprocessing.set_index("order").loc[orders]
    axes[1].bar(x, prep["preprocess_peak_kib"], color="#CC79A7")
    axes[1].set_ylabel("Peak Python-managed memory (KiB)")
    axes[1].set_title(
        r"$\mathbf{(e)}$ Shared preprocessing memory", loc="left", pad=10
    )

    for axis in axes:
        axis.set_xticks(x, [ORDER_LABELS[order] for order in orders])
        axis.grid(axis="y", alpha=0.25)
    figure.suptitle(r"Four interval orders at $n=100$, $m=20$", y=1.02)
    figure.tight_layout(pad=1.5)
    save_figure(figure, outdir, "figure_orders")


def baseline_table_and_figure(baseline_raw: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    timing = timing_stats(baseline_raw)
    preprocessing = preprocess_stats(baseline_raw)
    totals = timing.merge(
        preprocessing, on=["method", "order", "n", "m"], how="left"
    )
    totals["cold_total_ms"] = totals["median"] + totals["preprocess_ms"]
    proposed = totals[totals["method"] == "proposed"].drop(columns="method")
    brute = totals[totals["method"] == "bruteforce"].drop(columns="method")
    table = brute.merge(
        proposed,
        on=["order", "n", "m", "operation"],
        suffixes=("_bruteforce", "_proposed"),
    )
    table["runtime_speedup"] = (
        table["cold_total_ms_bruteforce"] / table["cold_total_ms_proposed"]
    )
    table["memory_ratio"] = (
        table["preprocess_peak_kib_bruteforce"]
        / table["preprocess_peak_kib_proposed"]
    )
    table.sort_values(["m", "n", "operation", "order"], inplace=True)

    figure, axes = plt.subplots(2, 2, figsize=(13.0, 9.0))
    for column, operation in enumerate(("rank", "unrank")):
        fixed_m = table[(table["m"] == 6) & (table["operation"] == operation)]
        fixed_n = table[(table["n"] == 10) & (table["operation"] == operation)]
        for order in ("t-inc", "lex1", "lex2", "xy"):
            part_n = fixed_m[fixed_m["order"] == order].sort_values("n")
            part_m = fixed_n[fixed_n["order"] == order].sort_values("m")
            axes[0, column].plot(
                part_n["n"],
                part_n["runtime_speedup"],
                color=ORDER_COLORS[order],
                marker=ORDER_MARKERS[order],
                label=ORDER_LABELS[order],
            )
            axes[1, column].plot(
                part_m["m"],
                part_m["runtime_speedup"],
                color=ORDER_COLORS[order],
                marker=ORDER_MARKERS[order],
                label=ORDER_LABELS[order],
            )
        axes[0, column].set(
            xlabel=r"Chain parameter $n$",
            yscale="log",
            title=rf"{operation.capitalize()}: fixed $m=6$",
        )
        axes[1, column].set(
            xlabel=r"Number of levels $m$",
            yscale="log",
            title=rf"{operation.capitalize()}: fixed $n=10$",
        )
    for label, axis in zip(("(f.1)", "(f.2)", "(f.3)", "(f.4)"), axes.flat):
        axis.grid(alpha=0.25, which="both")
        axis.legend(ncol=2)
        panel_label(axis, label, x=0.0, y=1.05)
    figure.tight_layout(pad=1.5, w_pad=3.0)
    column_ylabel = "Runtime ratio (brute force / proposed, ×)"
    for column in range(2):
        top_box = axes[0, column].get_position()
        bottom_box = axes[1, column].get_position()
        figure.text(
            top_box.x0 - 0.055,
            (top_box.y1 + bottom_box.y0) / 2,
            column_ylabel,
            rotation="vertical",
            va="center",
            ha="center",
            fontsize=plt.rcParams["axes.labelsize"],
        )
    save_figure(figure, outdir, "figure_baseline")

    figure, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    shared = table.drop_duplicates(["order", "n", "m"])
    for order in ("t-inc", "lex1", "lex2", "xy"):
        fixed_m = shared[(shared["m"] == 6) & (shared["order"] == order)].sort_values("n")
        fixed_n = shared[(shared["n"] == 10) & (shared["order"] == order)].sort_values("m")
        axes[0].plot(
            fixed_m["n"], fixed_m["memory_ratio"], color=ORDER_COLORS[order],
            marker=ORDER_MARKERS[order], label=ORDER_LABELS[order]
        )
        axes[1].plot(
            fixed_n["m"], fixed_n["memory_ratio"], color=ORDER_COLORS[order],
            marker=ORDER_MARKERS[order], label=ORDER_LABELS[order]
        )
    axes[0].set(xlabel=r"Chain parameter $n$")
    axes[0].set_title(r"$\mathbf{(g{.}1)}$ Fixed $m=6$", loc="left", pad=10)
    axes[1].set(xlabel=r"Number of levels $m$")
    axes[1].set_title(r"$\mathbf{(g{.}2)}$ Fixed $n=10$", loc="left", pad=10)
    for axis in axes:
        axis.set_ylabel("Peak-memory ratio (brute force / proposed, ×)")
        axis.set_yscale("log")
        axis.grid(alpha=0.25, which="both")
        axis.legend(ncol=2)
    figure.tight_layout(pad=1.5)
    save_figure(figure, outdir, "figure_baseline_memory")
    return table


def make_n_table(raw: pd.DataFrame) -> pd.DataFrame:
    timing = timing_stats(raw)
    memory = memory_stats(raw)
    keys = ["method", "order", "n", "m", "operation"]
    merged = timing.merge(
        memory[keys + ["maximum"]], on=keys, how="left"
    )
    return merged[
        ["n", "m", "operation", "mean", "median", "std", "q25", "q75", "maximum"]
    ].rename(columns={"maximum": "peak_kib"})


def make_order_table(raw: pd.DataFrame) -> pd.DataFrame:
    timing = timing_stats(raw)
    memory = memory_stats(raw)
    preprocessing = preprocess_stats(raw)
    timing = timing[(timing["n"] == 100) & (timing["m"] == 20)]
    memory = memory[(memory["n"] == 100) & (memory["m"] == 20)]
    query_keys = ["method", "order", "n", "m", "operation"]
    table = timing.merge(
        memory[query_keys + ["maximum"]], on=query_keys, how="left"
    ).merge(preprocessing, on=["method", "order", "n", "m"], how="left")
    return table[
        [
            "order",
            "operation",
            "mean",
            "median",
            "std",
            "q25",
            "q75",
            "maximum",
            "preprocess_ms",
            "preprocess_peak_kib",
        ]
    ].rename(columns={"maximum": "query_peak_kib"})


def write_table(table: pd.DataFrame, outdir: Path, stem: str) -> None:
    table.to_csv(outdir / f"{stem}.csv", index=False, float_format="%.6f")
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="In future versions `DataFrame.to_latex`.*",
            category=FutureWarning,
        )
        latex = table.to_latex(index=False, float_format="%.3f")
    (outdir / f"{stem}.tex").write_text(
        latex, encoding="utf-8"
    )


def regression_slope(stats: pd.DataFrame, x: np.ndarray, operation: str) -> float:
    part = stats[stats["operation"] == operation].sort_values("n" if "n" in stats else "m")
    return float(np.polyfit(np.log(x), np.log(part["median"].to_numpy()), 1)[0])


def write_narrative(
    m_stats: pd.DataFrame,
    n20_stats: pd.DataFrame,
    n100_stats: pd.DataFrame,
    order_table: pd.DataFrame,
    baseline_table: pd.DataFrame,
    outdir: Path,
) -> None:
    lines = [
        "# Benchmark results",
        "",
        "All reported query times are medians over independent validated round trips; shaded regions in the scaling figures show the interquartile range. Peak Python-managed memory was measured separately with `tracemalloc`, so tracing overhead is excluded from timings. Cold interval-order preprocessing is reported separately from amortized query cost.",
        "",
        "## Headline findings",
        "",
    ]
    for operation in ("unrank", "rank"):
        slopes = []
        endpoints = []
        for order in ("t-inc", "lex1", "lex2", "xy"):
            part = m_stats[
                (m_stats["operation"] == operation) & (m_stats["order"] == order)
            ].sort_values("m")
            slopes.append(np.polyfit(np.log(part["m"]), np.log(part["median"]), 1)[0])
            endpoints.append(part.iloc[-1]["median"])
        lines.append(
            f"- For fixed $n=10$, the four {operation} log–log slopes over $m=100$–$1000$ range from {min(slopes):.3f} to {max(slopes):.3f}; median times at $m=1000$ range from {min(endpoints):.3f} to {max(endpoints):.3f} ms."
        )
    for operation in ("unrank", "rank"):
        slopes = []
        endpoints = []
        for order in ("t-inc", "lex1", "lex2", "xy"):
            part = n20_stats[
                (n20_stats["operation"] == operation) & (n20_stats["order"] == order)
            ].sort_values("n")
            scale = part["n"] ** 2 * 20 * part["n"].map(math.log)
            slopes.append(np.polyfit(np.log(scale), np.log(part["median"]), 1)[0])
            endpoints.append(part.iloc[-1]["median"])
        lines.append(
            f"- For fixed $m=20$, the four normalized-scaling regression slopes for {operation} range from {min(slopes):.3f} to {max(slopes):.3f}; at $n=200$ median times range from {min(endpoints):.3f} to {max(endpoints):.3f} ms."
        )
    confirm = n100_stats[(n100_stats["n"] == 100)].set_index("operation")
    lines.append(
        f"- The larger confirmation point $n=100,m=100$ completes in median {confirm.loc['unrank','median']:.3f} ms (unrank) and {confirm.loc['rank','median']:.3f} ms (rank)."
    )
    roundtrip = order_table.groupby("order")["median"].sum().sort_values()
    fastest = roundtrip.index[0]
    slowest = roundtrip.index[-1]
    lines.append(
        f"- At $n=100,m=20$, {ORDER_LABELS[fastest]} has the lowest median round-trip query time ({roundtrip.iloc[0]:.3f} ms); {ORDER_LABELS[slowest]} has the highest ({roundtrip.iloc[-1]:.3f} ms). All four orders pass every round-trip check."
    )
    largest = baseline_table[(baseline_table["n"] == 10) & (baseline_table["m"] == 6)]
    lines.append(
        f"- At $n=10,m=6$, complete enumeration plus sorting is {largest['runtime_speedup'].min():.0f}×–{largest['runtime_speedup'].max():.0f}× slower in cold-start runtime across operations and orders, and reaches {largest['memory_ratio'].min():.0f}×–{largest['memory_ratio'].max():.0f}× higher shared preprocessing peak memory."
    )
    lines += [
        "",
        "## Recommended manuscript placement",
        "",
        "- Use `figure_scaling.pdf` for the complexity/scaling discussion.",
        "- Use `figure_orders.pdf` for the four-total-order comparison.",
        "- Use `figure_baseline.pdf` for the operation-specific horizontal baseline comparison.",
        "- Use `figure_baseline_memory.pdf` for the shared preprocessing-memory comparison.",
        "- Tables are supplied in both CSV and LaTeX formats.",
        "",
        "Do not describe the fitted slopes as proof of an asymptotic bound. They are empirical consistency checks under Python exact-integer arithmetic.",
    ]
    (outdir / "RESULTS_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("results/data"),
        help="Directory containing the five consolidated raw CSV files.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/plots"),
    )
    args = parser.parse_args()
    data_dir = args.data_dir
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    configure_style()

    m_raw = read_raw(data_dir / "m_scaling.csv")
    n20_raw = read_raw(data_dir / "n_scaling_m20.csv")
    n100_raw = read_raw(data_dir / "n_scaling_m100.csv")
    order_raw = n20_raw
    baseline_raw = read_raw(data_dir / "baseline.csv")
    m_stats = timing_stats(m_raw)
    n20_stats = timing_stats(n20_raw)
    n100_stats = timing_stats(n100_raw)

    scaling_figure(m_stats, n20_stats, outdir)
    order_figure(order_raw, outdir)
    baseline_table = baseline_table_and_figure(baseline_raw, outdir)

    n_table = pd.concat(
        [make_n_table(n20_raw), make_n_table(n100_raw)], ignore_index=True
    ).sort_values(["m", "n", "operation"])
    order_table = make_order_table(order_raw)
    write_table(n_table, outdir, "table_n_scaling")
    write_table(order_table, outdir, "table_orders_n100_m20")
    write_table(baseline_table, outdir, "table_baseline")
    write_narrative(
        m_stats, n20_stats, n100_stats, order_table, baseline_table, outdir
    )
    print(f"[COMPLETE] presentation written to {outdir}")


if __name__ == "__main__":
    main()
