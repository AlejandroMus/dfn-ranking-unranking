# -*- coding: utf-8 -*-
"""
DFN CUTS — Section-6–consistent t-inc (by α-cuts) + external module backend + full benchmark.
Clean CLI (no demo helpers), index base control, and CSV/plot generation.

- Global interval order is selected via --interval_order {lex1,lex2,xy,t-inc}.
- Orders lex1/lex2/xy are implemented in a Python module given by --engine/--module.
- The t-inc order is implemented internally by the Section-6 α-cut algorithm.
"""

import sys
import os
import time
import random
import argparse
import math
import importlib
import re
from dataclasses import dataclass
from typing import List, Tuple, Optional
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# ---------------------------------------------------------------------
# External module loader (lex1 / lex2 / xy / custom)
# ---------------------------------------------------------------------
try:
    import dfn_cuts_rank_unrank as _dfn_default
except Exception:
    # Fallback: add current and /mnt/data to sys.path
    for extra in [os.getcwd(), "/mnt/data"]:
        if extra not in sys.path:
            sys.path.append(extra)
    import dfn_cuts_rank_unrank as _dfn_default  # may still raise

# Active backend module used for module-based orders (lex1/lex2/xy/custom)
dcru = _dfn_default


def set_module(module_name: str) -> None:
    """
    Select the backend module implementing total_dfns / rank_dfn_cuts / unrank_dfn_cuts.

    - If module_name is empty or "dfn_cuts_rank_unrank", use the default companion module.
    - Otherwise, dynamically import the module with importlib.import_module(module_name).
    """
    global dcru
    if not module_name or module_name == "dfn_cuts_rank_unrank":
        dcru = _dfn_default
    else:
        dcru = importlib.import_module(module_name)


# ---------------------------------------------------------------------
# Basic utilities
# ---------------------------------------------------------------------
def ensure_outdir(path: str) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def seq_levels(seq: List[float], m: int) -> Tuple[int, ...]:
    """
    Quantize a membership sequence μ(i) in [0,1] into integer levels 0..m-1.
    """
    y_s = m - 1
    return tuple(int(round(v * y_s)) for v in seq)


def mu_form_from_levels(levels: Tuple[int, ...], m: int) -> str:
    """
    Pretty-print a DFN given as integer levels 0..m-1 into μ-form "μ/position".
    """
    y_s = m - 1 if m > 1 else 1
    items = [f"{(lv / y_s):.1f}/{x}" for x, lv in enumerate(levels)]
    return "{ " + ", ".join(items) + " }"


def parse_mu(mu_str: Optional[str], n: int, m: int) -> List[float]:
    """
    Parse a comma-separated μ-string into a list of floats of length n+1.
    """
    if mu_str is None:
        return [1, 1, 1, 0.4, 0.2, 0.2] if (n, m) == (5, 6) else [0.0] * (n + 1)
    vals = [float(x.strip()) for x in mu_str.split(",")]
    if len(vals) != n + 1:
        raise ValueError(f"--mu must have length n+1={n+1}")
    return vals


def parse_levels_arg(x: Optional[str]) -> Optional[List[int]]:
    """
    Parse levels or subindices from CLI.

    Accepts formats like:
      "0,1,3,5,3,1"   or   "[0, 1, 3, 5, 3, 1]"

    Returns:
        list[int] or None
    """
    if x is None:
        return None
    nums = re.findall(r"-?\d+", x)
    if not nums:
        raise ValueError(f"Cannot parse integer list from: {x!r}")
    return [int(t) for t in nums]


# ---------------------------------------------------------------------
# Comparator mapping for module-based engine
# ---------------------------------------------------------------------
def _effective_comp(order: str) -> str:
    """
    Map the global interval order (CLI) to the name understood by the backend module.

    - "lex1"   -> "lex1"
    - "lex2"   -> "lex2"
    - "xy"     -> "xy"     (Xu–Yager)
    - "engine" -> "engine" (custom order via engine.py)
    - "t-inc", "tinc", "t_inc" -> "tinc"  (if the module provides such a comparator)

    Any unknown value defaults to "lex1".
    """
    c = (order or "lex1").lower()
    if c in ("lex1", "lex2", "xy", "engine"):
        return c
    if c in ("tinc", "t-inc", "t_inc"):
        return "tinc"
    return "lex1"


def subindices_to_levels(subidx: List[int], n: int, m: int) -> List[int]:
    """
    Convert a list of subindices in {1,...,m} to internal levels 0..m-1.
    """
    if len(subidx) != n + 1:
        raise ValueError(f"subindex must have length n+1 = {n+1}")
    if any(s < 1 or s > m for s in subidx):
        raise ValueError(f"each subindex must be in 1..{m}")
    return [s - 1 for s in subidx]


# ---------------------------------------------------------------------
# DFN validators (discrete α-cuts and unimodality)
# ---------------------------------------------------------------------
def validate_dfn_levels(
    levels: List[int],
    n: int,
    m: int,
    require_normal: bool = True,
    check_unimodal: bool = False,
) -> List[str]:
    """
    Validate a Discrete Fuzzy Number (DFN) given as integer levels on L_n = {0,...,n}
    with m discrete membership levels {0,...,m-1}.

    Checks:
      1) Shape & range: length n+1, integers, each in [0, m-1].
      2) Normality (optional): max(levels) == m-1.
      3) α-cut convexity and nesting:
         - For each threshold thr = j-1 (j from m down to 2), the set { i : levels[i] >= thr }
           must be non-empty and contiguous (no holes).
         - α-cuts must be nested as j decreases (intervals expand outwards or stay equal).
      4) Unimodality (optional): non-decreasing up to the first maximum and
         non-increasing afterwards.

    Returns:
      A list of error messages. Empty list means "valid".
    """
    errs: List[str] = []

    # 1) type/length/range
    if len(levels) != n + 1:
        errs.append(f"subindex must have length n+1 = {n+1}")
    if any((not isinstance(l, (int, np.integer))) for l in levels):
        errs.append("subindex must be a list of integers")
    if any(l < 0 or l > m - 1 for l in levels):
        errs.append(f"each level must be in 0..{m-1}")

    # 2) normality (optional)
    if require_normal and (len(levels) == n + 1) and max(levels) != m - 1:
        errs.append("non-normal DFN: max(subindex) != m-1")

    # 3) α-cut convexity & nesting
    prev: Optional[Tuple[int, int]] = None  # previous inner α-cut interval [L,R]
    for j in range(m, 1, -1):  # j = m, m-1, ..., 2
        thr = j - 1
        idx = [i for i, lv in enumerate(levels) if lv >= thr]
        if not idx:
            errs.append(f"empty α-cut at j={j} (threshold={thr})")
            continue

        # contiguity (convexity): no holes in the α-cut indices
        if idx[-1] - idx[0] + 1 != len(idx):
            errs.append(f"non-convex α-cut (holes) at j={j}: indices={idx}")

        # nesting with previous α-cut: as threshold decreases, α-cuts must expand or stay equal
        if prev is not None:
            L, R = idx[0], idx[-1]
            pL, pR = prev
            if not (L <= pL and R >= pR):
                errs.append(
                    f"non-nested alpha-cuts at j={j} vs j={j+1}: "
                    f"[{L},{R}] does not contain [{pL},{pR}]"
                )
        prev = (idx[0], idx[-1])

    # 4) unimodality (optional)
    if check_unimodal and (len(levels) == n + 1):
        pmax = max(levels)
        i0 = levels.index(pmax)  # first occurrence of the maximum
        left_ok = all(levels[i] <= levels[i + 1] for i in range(0, i0))
        right_ok = all(levels[i] >= levels[i + 1] for i in range(i0, len(levels) - 1))
        if not (left_ok and right_ok):
            errs.append("non-unimodal DFN (not monotone up then monotone down)")

    return errs


def validate_dfn_mu(mu: List[float], n: int, m: int, **kw) -> List[str]:
    """
    Validate a DFN given as a membership vector μ (length n+1, values in [0,1]).
    Internally discretizes μ into integer levels via seq_levels and delegates to
    validate_dfn_levels.
    """
    errs: List[str] = []
    if len(mu) != n + 1:
        errs.append(f"μ must have length n+1 = {n+1}")
    if any((v < 0 or v > 1) for v in mu):
        errs.append("some μ[i] is outside [0,1]")

    try:
        levels = list(seq_levels(mu, m))
    except Exception as e:
        errs.append(f"discretization failed in seq_levels: {e}")
        return errs

    errs += validate_dfn_levels(levels, n, m, **kw)
    return errs


def ensure_valid_levels(levels: List[int], n: int, m: int, **kw) -> List[int]:
    """
    Raise ValueError if levels do not represent a valid DFN.
    """
    errs = validate_dfn_levels(levels, n, m, **kw)
    if errs:
        raise ValueError("Invalid DFN (subindex):\n- " + "\n- ".join(errs))
    return levels


def ensure_valid_mu(mu: List[float], n: int, m: int, **kw) -> List[float]:
    """
    Raise ValueError if μ does not represent a valid DFN.
    """
    errs = validate_dfn_mu(mu, n, m, **kw)
    if errs:
        raise ValueError("Invalid DFN (mu):\n- " + "\n- ".join(errs))
    return mu


# ---------------------------------------------------------------------
# Internal t-inc interval ordering (Section 6)
# ---------------------------------------------------------------------
def intervals_tinc_sorted(n: int) -> List[Tuple[int, int]]:
    """
    t-inc order on intervals [a,b] with a <= b in {0,...,n}:
      - a increases  0..n
      - for each a,  b decreases n..a
    """
    order: List[Tuple[int, int]] = []
    for a in range(0, n + 1):
        for b in range(n, a - 1, -1):
            order.append((a, b))
    return order


def sdfn_count_at_level(n: int, j: int, a: int, b: int) -> int:
    """
    |SDFN(a,b,j)| = C(a+T,T) * C((n-b)+T,T) with T = j-2
    (stars-and-bars counts for left/right extensions).
    """
    T = max(0, j - 2)
    return math.comb(a + T, T) * math.comb((n - b) + T, T)


# ---------------------------------------------------------------------
# Section-6 α-cut helpers
# ---------------------------------------------------------------------
def alpha_intervals_from_levels(
    levels: List[int], n: int, m: int
) -> List[Tuple[int, int]]:
    """
    For j = m..2, return [a_j,b_j] where A_{y_j} = {x : level(x) >= j-1} = [a_j,b_j].
    """
    res: List[Tuple[int, int]] = []
    for j in range(m, 1, -1):
        thr = j - 1
        indices = [x for x, lv in enumerate(levels) if lv >= thr]
        if not indices:
            raise ValueError(f"Empty α-cut at level j={j}; not a valid CUTS DFN.")
        res.append((indices[0], indices[-1]))
    return res  # length m-1, entries for j = m,...,2


def levels_from_alpha_intervals(
    alpha_list: List[Tuple[int, int]], n: int, m: int
) -> List[int]:
    """
    From nested intervals [a_j,b_j] for j=m..2, rebuild integer levels (0..m-1).
    """
    lv = [0] * (n + 1)
    for idx, (a, b) in enumerate(alpha_list):
        j = m - idx
        for x in range(a, b + 1):
            lv[x] = max(lv[x], j - 1)
    return lv


# ---------------------------------------------------------------------
# UNRANK (t-inc) — by cuts (algorithm + optional log)
# ---------------------------------------------------------------------
def unrank_tinc_by_cuts(
    n: int,
    m: int,
    i0: int,
    index_base: int = 1,
    show_both: bool = False,
) -> Tuple[List[float], List[Tuple[int, int]], str]:
    """
    Section-6 UNRANK under the t-inc order (index i0 is 0-based internally).

    Returns:
        (seq, alpha_list, log_text) where:
        - seq        : membership vector μ(i) in [0,1],
        - alpha_list : list of intervals [a_j,b_j] for j = m..2,
        - log_text   : detailed human-readable trace.
    """
    lines: List[str] = []
    total = dcru.total_dfns(n, m)
    i_disp = i0 if index_base == 0 else i0 + 1
    min_disp = 0 if index_base == 0 else 1
    max_disp = (total - 1) if index_base == 0 else total

    lines.append("=== t-inc UNRANK (by α-cuts) ===")
    lines.append(f"Parameters: n={n}, m={m}, order=t-inc")
    lines.append(
        f"Index i ({'0-based' if index_base == 0 else '1-based'}) = {i_disp} "
        f"of {min_disp}..{max_disp}"
    )
    if show_both:
        lines.append(f"(also: 0-based={i0} / 1-based={i0 + 1})")

    # Level j = m: choose the core interval
    order = intervals_tinc_sorted(n)
    acc = 0
    a_m = b_m = None
    lines.append(f"\nLevel j={m} candidates: (a asc, b desc)")
    lines.append(" idx | [a,b] | |SDFN(a,b,j)| | accumulated")
    for idx, (a, b) in enumerate(order, start=1):
        cnt = sdfn_count_at_level(n, m, a, b)
        acc2 = acc + cnt
        mark = ""
        if a_m is None and acc2 > i0:
            a_m, b_m = a, b
            mark = "  <-- pick"
        lines.append(
            f"{idx:>4} | [{a},{b}] | {cnt:>12}    | {acc2:>11}{mark}"
        )
        acc = acc2
        if a_m is not None and acc2 > i0:
            break
    before = acc - sdfn_count_at_level(n, m, a_m, b_m)  # type: ignore[arg-type]
    resid = i0 - before
    lines.append(
        f"\nChosen core: [{a_m},{b_m}] ; subtract previous = {before} → residual i0 = {resid}"
    )
    alpha: List[Tuple[int, int]] = [(a_m, b_m)]  # type: ignore[list-item]

    # Lower levels j = m-1,...,2
    prev_a, prev_b = a_m, b_m  # type: ignore[assignment]
    for j in range(m - 1, 1, -1):
        acc = 0
        lines.append(f"\nLevel j={j} candidates containing [{prev_a},{prev_b}]:")
        lines.append(" idx | [a,b] | |SDFN(a,b,j)| | accumulated")
        idx_sel = None
        chosen: Optional[Tuple[int, int]] = None
        candidates = [(a, b) for (a, b) in order if a <= prev_a and b >= prev_b]
        for idx, (a, b) in enumerate(candidates, start=1):
            cnt = sdfn_count_at_level(n, j, a, b)
            acc2 = acc + cnt
            mark = ""
            if idx_sel is None and acc2 > resid:
                idx_sel = idx
                chosen = (a, b)
                mark = "  <-- pick"
            lines.append(
                f"{idx:>4} | [{a},{b}] | {cnt:>12}    | {acc2:>11}{mark}"
            )
            acc = acc2
            if chosen is not None and acc2 > resid:
                break
        before2 = acc - sdfn_count_at_level(
            n, j, chosen[0], chosen[1]  # type: ignore[index]
        )
        resid = resid - before2
        lines.append(
            f"Chosen j={j}: [{chosen[0]},{chosen[1]}] ; subtract previous = {before2} → residual i0 = {resid}"  # type: ignore[index]
        )
        alpha.append(chosen)  # type: ignore[arg-type]
        prev_a, prev_b = chosen  # type: ignore[assignment]

    # Rebuild μ from α-cuts
    lv = levels_from_alpha_intervals(alpha, n, m)
    seq = [x / (m - 1) for x in lv]
    subidx_paper = tuple(lv_i + 1 for lv_i in lv)
    lines.append("\nFinal DFN subindex (1..m): " + str(subidx_paper))
    lines.append("μ-form: " + mu_form_from_levels(tuple(lv), m))

    return seq, alpha, "\n".join(lines)


# ---------------------------------------------------------------------
# RANK (t-inc) — by cuts (algorithm + optional log)
# ---------------------------------------------------------------------
def rank_tinc_by_cuts_from_levels(
    levels: List[int],
    n: int,
    m: int,
    index_base: int = 1,
    show_both: bool = False,
) -> Tuple[int, str]:
    """
    Section-6 RANK under the t-inc order.

    Args:
        levels     : integer levels 0..m-1 (length n+1).
        index_base : output base {0,1}.
        show_both  : if True, both 0-based and 1-based indices appear in the log.

    Returns:
        (i0, log_text) where i0 is the internal 0-based index.
    """
    ensure_valid_levels(levels, n, m, require_normal=True, check_unimodal=False)

    lines: List[str] = []
    lines.append("=== t-inc RANK (by α-cuts) ===")
    lines.append(f"Parameters: n={n}, m={m}, order=t-inc")
    lines.append("Input μ: " + mu_form_from_levels(tuple(levels), m))

    order = intervals_tinc_sorted(n)
    alpha = alpha_intervals_from_levels(levels, n, m)

    i0 = 0
    prev_a: Optional[int] = None
    prev_b: Optional[int] = None
    for idx_level, (a_sel, b_sel) in enumerate(alpha):
        j = m - idx_level
        if prev_a is None:
            candidates = order
            lines.append(f"\nLevel j={j} candidates:")
        else:
            candidates = [(a, b) for (a, b) in order if a <= prev_a and b >= prev_b]  # type: ignore[arg-type]
            lines.append(
                f"\nLevel j={j} candidates containing [{prev_a},{prev_b}]:"
            )
        lines.append(" idx | [a,b] | |SDFN(a,b,j)| | cumulative sum added to i")
        cum = 0
        for k, (a, b) in enumerate(candidates, start=1):
            cnt = sdfn_count_at_level(n, j, a, b)
            if (a, b) == (a_sel, b_sel):
                lines.append(
                    f"{k:>4} | [{a},{b}] | {cnt:>12}    | (stop here)"
                )
                break
            else:
                i0 += cnt
                cum += cnt
                lines.append(
                    f"{k:>4} | [{a},{b}] | {cnt:>12}    | +{cum}"
                )
        lines.append(f"Chosen j={j}: [{a_sel},{b_sel}] → partial i0 = {i0}")
        prev_a, prev_b = a_sel, b_sel

    i_disp = i0 if index_base == 0 else i0 + 1
    base_label = "0-based" if index_base == 0 else "1-based"
    lines.append(f"\nResult: pos^-1(A) = i ({base_label}) = {i_disp}")
    if show_both:
        lines.append(f"(also: 0-based={i0} / 1-based={i0 + 1})")

    return i0, "\n".join(lines)


# ---------------------------------------------------------------------
# Benchmark utilities
# ---------------------------------------------------------------------
@dataclass
class BenchRow:
    m: int
    mean_unrank_ms: float
    std_unrank_ms: float
    mean_rank_ms: float
    std_rank_ms: float


def _now_ms() -> float:
    return time.time() * 1000.0


def _unrank_engine(order: str, n: int, m: int, i0: int) -> List[float]:
    """
    Compute μ for index i0 (0-based) under a given global order.

    - order in {"lex1","lex2","xy","engine"} → use module-based backend dcru.
    - order ~ "t-inc"                        → use internal t-inc engine.
    """
    c = (order or "t-inc").lower()
    if c in ("t-inc", "tinc", "t_inc"):
        # Internal t-inc by α-cuts (we ignore logs in the benchmark)
        seq, _, _ = unrank_tinc_by_cuts(n, m, i0, index_base=0, show_both=False)
        return seq
    else:
        comp = _effective_comp(c)
        return dcru.unrank_dfn_cuts(n, m, i0, comp)


def _rank_engine(order: str, n: int, m: int, levels: List[int]) -> int:
    """
    Compute the 0-based index i0 from levels under a given global order.

    - order in {"lex1","lex2","xy","engine"} → use module-based backend dcru.
    - order ~ "t-inc"                        → use internal t-inc engine.
    """
    c = (order or "t-inc").lower()
    if c in ("t-inc", "tinc", "t_inc"):
        i0, _ = rank_tinc_by_cuts_from_levels(
            levels, n, m, index_base=0, show_both=False
        )
        return i0
    else:
        ensure_valid_levels(levels, n, m, require_normal=True, check_unimodal=False)
        seq = [lv / float(m - 1) for lv in levels]
        comp = _effective_comp(c)
        return dcru.rank_dfn_cuts(n, m, seq, comp)


def benchmark(
    outdir: Path,
    n: int = 10,
    m_from: int = 10,
    m_to: int = 100,
    m_step: int = 10,
    trials: int = 100,
    interval_order: str = "t-inc",
) -> pd.DataFrame:
    """
    Run timing benchmarks for a fixed n and varying m.

    - interval_order in {"lex1","lex2","xy","engine"} → module-based backend.
    - interval_order == "t-inc"                       → internal t-inc engine.
    """
    rows: List[dict] = []
    order = interval_order

    for m in range(m_from, m_to + 1, m_step):
        unrank_ms: List[float] = []
        rank_ms: List[float] = []
        total = dcru.total_dfns(n, m)
        rng = random.Random(12345 + m)

        for _ in range(trials):
            i0 = rng.randrange(total)

            # UNRANK
            t0 = _now_ms()
            mu = _unrank_engine(order, n, m, i0)
            t1 = _now_ms()
            unrank_ms.append(t1 - t0)

            # RANK
            levels = list(seq_levels(mu, m))
            t2 = _now_ms()
            idx = _rank_engine(order, n, m, levels)
            t3 = _now_ms()
            rank_ms.append(t3 - t2)

            if idx != i0:
                raise RuntimeError(f"rank/unrank are not inverse for i0={i0}, idx={idx}")

        rows.append(
            {
                "m": m,
                "mean_unrank_ms": float(np.mean(unrank_ms)),
                "std_unrank_ms": float(np.std(unrank_ms, ddof=1)),
                "mean_rank_ms": float(np.mean(rank_ms)),
                "std_rank_ms": float(np.std(rank_ms, ddof=1)),
            }
        )

    df = pd.DataFrame(rows)

    # Save CSV with both metrics
    df_csv = df.rename(
        columns={
            "mean_unrank_ms": "unrank t̄m (ms)",
            "std_unrank_ms": "unrank sm (ms)",
            "mean_rank_ms": "rank t̄m (ms)",
            "std_rank_ms": "rank sm (ms)",
        }
    )
    (outdir / "bench_times_unrank_rank.csv").write_text(
        df_csv.to_csv(index=False, float_format="%.3f"), encoding="utf-8"
    )

    # (a) Average execution time vs m
    plt.figure()
    plt.plot(df["m"], df["mean_unrank_ms"], marker="o", label="unranking")
    plt.plot(df["m"], df["mean_rank_ms"], marker="s", label="ranking")
    plt.xlabel(r"$m$")
    plt.ylabel(r"Average execution time $\bar{t}_m$ (ms)")
    plt.title(
        r"Average execution time as a function of $m$ (for fixed $n$) "
        f"\n(n={n}, trials={trials}, order={interval_order})"
    )
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outdir / "fig_times_unrank_rank.png", dpi=160)
    plt.close()

    # (b) Log–log figure and slope fitting
    x = df["m"].to_numpy(dtype=float)
    y_u = df["mean_unrank_ms"].to_numpy(dtype=float)
    y_r = df["mean_rank_ms"].to_numpy(dtype=float)
    mask_u = (x > 0) & (y_u > 0) & np.isfinite(y_u)
    mask_r = (x > 0) & (y_r > 0) & np.isfinite(y_r)

    logx_u = np.log(x[mask_u])
    logy_u = np.log(y_u[mask_u])
    slope_u, intercept_u = np.polyfit(logx_u, logy_u, 1)

    logx_r = np.log(x[mask_r])
    logy_r = np.log(y_r[mask_r])
    slope_r, intercept_r = np.polyfit(logx_r, logy_r, 1)

    plt.figure()
    plt.plot(
        logx_u,
        logy_u,
        marker="o",
        linestyle="",
        label=f"unranking (slope={slope_u:.3f})",
    )
    plt.plot(
        logx_r,
        logy_r,
        marker="s",
        linestyle="",
        label=f"ranking (slope={slope_r:.3f})",
    )
    plt.plot(logx_u, intercept_u + slope_u * logx_u, label="fit unranking")
    plt.plot(logx_r, intercept_r + slope_r * logx_r, label="fit ranking")
    plt.xlabel(r"$\log m$")
    plt.ylabel(r"$\log \bar{t}_m$")
    plt.title(
        r"Log–log plot of average runtime vs.~$m$ showing fitted slope $\approx 1$"
    )
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outdir / "fig_loglog_unrank_rank.png", dpi=160)
    plt.close()

    latex_table = df_csv.to_latex(
        index=False,
        escape=False,
        float_format="%.3f",
        caption=(
            r"Average execution time as a function of $m$ (for fixed $n=10$). "
            rf"Results over $K={trials}$ trials; order={interval_order}."
        ),
        label="tab:table1",
        column_format="r r r r r",
    )
    (outdir / "table1_results.tex").write_text(latex_table, encoding="utf-8")
    print(f"[Table 1 → {outdir/'table1_results.tex'}]")

    return df


# ---------------------------------------------------------------------
# Plotter from CSV
# ---------------------------------------------------------------------
def make_plots_from_csv(
    csv_path: str,
    outdir: str = "out",
    n: int = 10,
    trials: int = 500,
    order_label: str = "t-inc",
    fig_name_times: str = "fig_times_unrank_rank.png",
    fig_name_loglog: str = "fig_loglog_unrank_rank.png",
):
    """
    Load a CSV with timing information and regenerate:

      (a) Average execution time vs m
      (b) Log–log plot with fitted slopes

    The function automatically detects whether there is one or two series
    (unranking, ranking) in the CSV.
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(csv_path)

    cols = {c.lower().strip(): c for c in df.columns}
    col_m = cols.get("m", "m")

    # candidates for a single series
    one_series_candidates = [
        "mean_ms",
        "t̄m (ms)",
        "tm (ms)",
        "mean_time_ms",
        "average_ms",
    ]

    # candidates for two-series CSV
    unr_mean = next(
        (
            cols[k]
            for k in [
                "mean_unrank_ms",
                "unrank t̄m (ms)",
                "unrank mean (ms)",
                "unrank mean_ms",
            ]
            if k in cols
        ),
        None,
    )
    unr_std = next(
        (
            cols[k]
            for k in [
                "std_unrank_ms",
                "unrank sm (ms)",
                "unrank std (ms)",
                "unrank std_ms",
            ]
            if k in cols
        ),
        None,
    )
    rnk_mean = next(
        (
            cols[k]
            for k in [
                "mean_rank_ms",
                "rank t̄m (ms)",
                "rank mean (ms)",
                "rank mean_ms",
            ]
            if k in cols
        ),
        None,
    )
    rnk_std = next(
        (
            cols[k]
            for k in [
                "std_rank_ms",
                "rank sm (ms)",
                "rank std (ms)",
                "rank std_ms",
            ]
            if k in cols
        ),
        None,
    )

    single_mean = None
    if not unr_mean and not rnk_mean:
        for k in one_series_candidates:
            if k in cols:
                single_mean = cols[k]
                break
        if not single_mean:
            for c in df.columns:
                if "t" in c and "ms" in c and "̄" in c:
                    single_mean = c
                    break
        if not single_mean:
            raise ValueError(
                "Could not find columns with timing data in the CSV."
            )

    # (a) time vs m
    plt.figure()
    x = df[col_m].to_numpy(dtype=float)

    if unr_mean or rnk_mean:
        # two series
        if unr_mean:
            y_u = df[unr_mean].to_numpy(dtype=float)
            plt.plot(x, y_u, marker="o", label="unranking")
        if rnk_mean:
            y_r = df[rnk_mean].to_numpy(dtype=float)
            plt.plot(x, y_r, marker="s", label="ranking")
        plt.legend()
    else:
        # one series
        y = df[single_mean].to_numpy(dtype=float)
        plt.plot(x, y, marker="o", label="unranking (single series)")
        plt.legend()

    plt.xlabel(r"$m$")
    plt.ylabel(r"Average execution time $\bar{t}_m$ (ms)")
    plt.title(
        r"Average execution time as a function of $m$ (for fixed $n$) "
        f"(n={n}, trials={trials}, order={order_label})"
    )
    plt.grid(True)
    plt.tight_layout()
    path_times = out / fig_name_times
    plt.savefig(path_times, dpi=160)
    plt.close()

    # (b) log–log + fits
    plt.figure()

    def _fit_and_plot(series: np.ndarray, marker: str, label_prefix: str):
        mask = (
            np.isfinite(x)
            & (x > 0)
            & np.isfinite(series)
            & (series > 0)
        )
        logx = np.log(x[mask])
        logy = np.log(series[mask])
        slope, intercept = np.polyfit(logx, logy, 1)
        plt.plot(
            logx,
            logy,
            marker=marker,
            linestyle="",
            label=f"{label_prefix} (slope={slope:.3f})",
        )
        plt.plot(
            logx,
            intercept + slope * logx,
            label=f"fit {label_prefix}",
        )
        return slope, intercept

    lines_info = []
    if unr_mean or rnk_mean:
        if unr_mean:
            y_u = df[unr_mean].to_numpy(dtype=float)
            su, iu = _fit_and_plot(y_u, "o", "unranking")
            lines_info.append(("unranking", su, iu))
        if rnk_mean:
            y_r = df[rnk_mean].to_numpy(dtype=float)
            sr, ir = _fit_and_plot(y_r, "s", "ranking")
            lines_info.append(("ranking", sr, ir))
    else:
        y = df[single_mean].to_numpy(dtype=float)
        s, b = _fit_and_plot(y, "o", "unranking")
        lines_info.append(("unranking", s, b))

    plt.xlabel(r"$\log m$")
    plt.ylabel(r"$\log \bar{t}_m$")
    plt.title(
        r"Log–log plot of average runtime vs.~$m$ showing fitted slope $\approx 1$"
    )
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    path_loglog = out / fig_name_loglog
    plt.savefig(path_loglog, dpi=160)
    plt.close()

    info_txt = "\n".join(
        [
            f"{name} slope = {s:.6f} ; intercept = {b:.6f}"
            for (name, s, b) in lines_info
        ]
    )
    (out / "bench_fit_fromcsv.txt").write_text(info_txt, encoding="utf-8")

    print(f"[OK] Plots from CSV →\n  {path_times}\n  {path_loglog}")


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=(
            "Instructions to use the algorithm"
        )
    )
    ap.add_argument(
        "--unrank",
        action="store_true",
        help="Run UNRANK (index -> DFN) under the chosen interval order.",
    )
    ap.add_argument(
        "--rank",
        action="store_true",
        help="Run RANK (DFN -> index) under the chosen interval order.",
    )
    ap.add_argument("--n", type=int, default=5, help="Length of L_n (we work on {0,...,n}).")
    ap.add_argument("--m", type=int, default=6, help="Number of membership levels |Y_m| = m.")
    ap.add_argument(
        "--i",
        type=int,
        default=50,
        help="Index for unrank, interpreted in the base given by --index_base.",
    )
    ap.add_argument(
        "--mu",
        type=str,
        default=None,
        help='Comma-separated μ for rank, e.g. "1,1,1,0.4,0.2,0.2".',
    )
    ap.add_argument(
        "--subindex",
        type=str,
        default=None,
        help='Comma-separated integer subindices (1..m) for rank, e.g. "6,6,6,2,1,1".',
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="If given, show detailed t-inc α-cut log (only for interval_order=t-inc).",
    )
    ap.add_argument(
        "--bench",
        action="store_true",
        help="Run timing benchmarks (results in CSV + plots + LaTeX table).",
    )
    ap.add_argument(
        "--interval_order",
        type=str,
        default="t-inc",
        choices=["lex1", "lex2", "xy", "t-inc", "engine"],
        help=(
            "Global interval order: lex1, lex2, xy (Xu–Yager), "
            "t-inc (internal), or engine (custom order from engine.py)."
        ),
    )
    # module backend (engine); alias --module to match the README text
    ap.add_argument(
        "--engine",
        "--module",
        dest="engine",
        type=str,
        default="dfn_cuts_rank_unrank",
        help=(
            "Python module implementing rank/unrank for orders lex1/lex2/xy. "
            "Default: dfn_cuts_rank_unrank."
        ),
    )
    ap.add_argument(
        "--m_from", type=int, default=10, help="Initial m value for the benchmark."
    )
    ap.add_argument(
        "--m_to", type=int, default=100, help="Final m value for the benchmark."
    )
    ap.add_argument(
        "--m_step", type=int, default=10, help="Step size for m in the benchmark."
    )
    ap.add_argument(
        "--trials",
        type=int,
        default=100,
        help="Number of random rank/unrank queries per (n,m) configuration.",
    )
    ap.add_argument(
        "--outdir",
        type=str,
        default="out",
        help="Output directory for logs, CSVs, plots, and LaTeX tables.",
    )
    ap.add_argument(
        "--index_base",
        type=int,
        choices=[0, 1],
        default=1,
        help="Index base in the CLI: 0 (internal) or 1 (matches the paper).",
    )
    ap.add_argument(
        "--show_both_indices",
        action="store_true",
        help="When printing detailed t-inc logs, also show both 0-based and 1-based indices.",
    )
    ap.add_argument(
        "--plots_from_csv",
        type=str,
        default=None,
        help="Path to a CSV to regenerate plots without re-running the benchmark.",
    )
    ap.add_argument(
        "--plots_order_label",
        type=str,
        default="t-inc",
        help="Order label that will appear in plot titles (for plots_from_csv).",
    )

    args = ap.parse_args()
    outdir = ensure_outdir(args.outdir)

    # Select backend module for lex1/lex2/xy orders
    set_module(args.engine)

    did = False

    # ---------------- UNRANK ----------------
    if args.unrank:
        did = True
        total = dcru.total_dfns(args.n, args.m)
        i0 = args.i if args.index_base == 0 else args.i - 1
        if not (0 <= i0 < total):
            raise ValueError(
                f"Index out of range after base conversion: i0={i0}, valid 0..{total-1}"
            )

        order = args.interval_order.lower()

        if order in ("t-inc", "tinc", "t_inc"):
            # Internal t-inc UNRANK
            seq, alpha, log = unrank_tinc_by_cuts(
                args.n,
                args.m,
                i0,
                index_base=args.index_base,
                show_both=args.show_both_indices,
            )
            if args.verbose:
                print(log)
        else:
            # Module-based order
            comp = _effective_comp(order)
            seq = dcru.unrank_dfn_cuts(args.n, args.m, i0, comp)

        levels = tuple(seq_levels(seq, args.m))
        subidx = [lv + 1 for lv in levels]

        print("DFN (subindex 1..m):", subidx)
        print("DFN (μ-form):", mu_form_from_levels(levels, args.m))

    # ---------------- RANK ----------------
    if args.rank:
        did = True
        order = args.interval_order.lower()

        # Input DFN: either from subindices or from μ
        if args.subindex is not None:
            subidx_list = parse_levels_arg(args.subindex)
            levels = subindices_to_levels(subidx_list, args.n, args.m)
        else:
            mu = parse_mu(args.mu, args.n, args.m)
            seq = mu
            levels = list(seq_levels(seq, args.m))

        if order in ("t-inc", "tinc", "t_inc"):
            # Internal t-inc RANK
            i0, log = rank_tinc_by_cuts_from_levels(
                levels,
                args.n,
                args.m,
                index_base=args.index_base,
                show_both=args.show_both_indices,
            )
            if args.verbose:
                print(log)
        else:
            ensure_valid_levels(
                levels, args.n, args.m, require_normal=True, check_unimodal=False
            )
            seq = [lv / float(args.m - 1) for lv in levels]
            comp = _effective_comp(order)
            i0 = dcru.rank_dfn_cuts(args.n, args.m, seq, comp)

        i_disp = i0 if args.index_base == 0 else i0 + 1
        base_label = "0-based" if args.index_base == 0 else "1-based"
        print(f"Index i ({base_label}) = {i_disp}")

    # ---------------- BENCHMARK ----------------
    if args.bench:
        did = True
        df = benchmark(
            outdir,
            n=args.n,
            m_from=args.m_from,
            m_to=args.m_to,
            m_step=args.m_step,
            trials=args.trials,
            interval_order=args.interval_order,
        )
        print("\n== BENCHMARK RESULTS ==")
        print(df.to_string(index=False))

    # ---------------- PLOTS FROM CSV ----------------
    if args.plots_from_csv:
        did = True
        make_plots_from_csv(
            csv_path=args.plots_from_csv,
            outdir=args.outdir,
            n=args.n,
            trials=args.trials,
            order_label=args.plots_order_label,
        )

    if not did:
        print("Nothing selected. Use --unrank, --rank, --bench, or --plots_from_csv.")


if __name__ == "__main__":
    main()
