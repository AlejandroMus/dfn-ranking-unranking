# -*- coding: utf-8 -*-
"""
DFN CUTS — Section-6–consistent t-inc (by α-cuts) + module path + full benchmark + clean index base.
v7
----
- Adds --index_base {0,1} (default 1 = paper) and --show_both_indices (off by default).
- Demos interpret --i in the chosen base and print only that base (no mixed lines).
- Optional dual display only if --show_both_indices is provided.
- Keeps v6 functionality: module/tinc engines, full bench, detailed logs.
"""

import sys, os, time, random, argparse, math
from dataclasses import dataclass
from typing import List, Tuple, Iterable, Optional
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ------------- Module path (lex1/lex2/scopia) -------------
try:
    import dfn_cuts_rank_unrank as dcru
except Exception:
    for extra in [os.getcwd(), "/mnt/data"]:
        if extra not in sys.path:
            sys.path.append(extra)
    import dfn_cuts_rank_unrank as dcru  # may still raise

# ----------------------- Utilities -----------------------
def ensure_outdir(path: str) -> Path:
    p = Path(path); p.mkdir(parents=True, exist_ok=True); return p

def seq_levels(seq: List[float], m: int) -> Tuple[int, ...]:
    y_s = m - 1
    return tuple(int(round(v * y_s)) for v in seq)

def mu_form_from_levels(levels: Tuple[int, ...], m: int) -> str:
    y_s = m - 1 if m > 1 else 1
    items = [f"{(lv/y_s):.1f}/{x}" for x, lv in enumerate(levels)]
    return "{ " + ", ".join(items) + " }"

def parse_mu(mu_str: Optional[str], n: int, m: int) -> List[float]:
    if mu_str is None:
        return [1,1,1,0.4,0.2,0.2] if (n,m)==(5,6) else [0.0]*(n+1)
    vals = [float(x.strip()) for x in mu_str.split(",")]
    if len(vals) != n+1:
        raise ValueError(f"--mu must have length n+1={n+1}")
    return vals

# Comparator mapping for module engine
def _effective_comp(comp: str) -> str:
    c = (comp or 'lex1').lower()
    if c == 'scopia': return 'lex2'
    if c in ('lex1','lex2','tinc'): return c
    return 'lex1'

# ------------------ t-inc interval ordering ------------------
def intervals_tinc_sorted(n: int) -> List[Tuple[int,int]]:
    """t-inc: a asc, b desc within each a, for all a<=b in [0..n]."""
    order = []
    for a in range(0, n+1):
        for b in range(n, a-1, -1):
            order.append((a,b))
    return order

def sdfn_count_at_level(n: int, j: int, a: int, b: int) -> int:
    """|SDFN(a,b,j)| = C(a+T,T) * C((n-b)+T,T) with T=j-2."""
    T = max(0, j - 2)
    return math.comb(a + T, T) * math.comb((n - b) + T, T)

# ---------------- Section-6 α-cut helpers ----------------
def alpha_intervals_from_levels(levels: List[int], n: int, m: int) -> List[Tuple[int,int]]:
    """
    For j=m..2, return [a_j,b_j] where A_{y_j} = {x : level(x) >= j-1} = [a_j,b_j].
    """
    res = []
    for j in range(m, 1, -1):
        thr = j - 1
        indices = [x for x, lv in enumerate(levels) if lv >= thr]
        if not indices:
            raise ValueError(f"Empty α-cut at level j={j}; not a valid CUTS DFN.")
        res.append((indices[0], indices[-1]))
    return res  # length m-1 entries: j=m,..,2

def levels_from_alpha_intervals(alpha_list: List[Tuple[int,int]], n: int, m: int) -> List[int]:
    """
    From nested intervals [a_j,b_j] for j=m..2, rebuild integer levels (0..m-1).
    """
    lv = [0]*(n+1)
    for idx, (a,b) in enumerate(alpha_list):
        j = m - idx
        for x in range(a, b+1):
            lv[x] = max(lv[x], j-1)
    return lv

# --------------- UNRANK (t-inc) — by cuts ---------------
def unrank_tinc_by_cuts(n: int, m: int, i0: int, index_base: int = 1, show_both: bool = False) -> Tuple[List[float], List[Tuple[int,int]], str]:
    """
    True Section-6 unrank (internal index is 0-based i0). Logs in the chosen base.
    Returns: (seq, alpha_list[j=m..2], log_text)
    """
    lines = []
    total = dcru.total_dfns(n, m)
    i_disp = i0 if index_base == 0 else i0 + 1
    min_disp = 0 if index_base == 0 else 1
    max_disp = (total - 1) if index_base == 0 else total
    lines.append("=== Academic Example — t-inc UNRANK (by cuts) ===")
    lines.append(f"Parameters: n={n}, m={m}, comp=tinc")
    lines.append(f"Index i ({'0-based' if index_base==0 else '1-based'}) = {i_disp}  of {min_disp}..{max_disp}")
    if show_both:
        lines.append(f"(also: 0-based={i0} / 1-based={i0+1})")

    # j = m
    order = intervals_tinc_sorted(n)
    acc = 0
    a_m = b_m = None
    lines.append(f"\nLevel j={m} candidates: (a asc, b desc)")
    lines.append(" idx | [a,b] | |SDFN(a,b,j)| | accumulated")
    for idx, (a,b) in enumerate(order, start=1):
        cnt = sdfn_count_at_level(n, m, a, b); acc2 = acc + cnt
        mark = ""
        if a_m is None and acc2 > i0:
            a_m, b_m = a, b
            mark = "  <-- pick"
        lines.append(f"{idx:>4} | [{a},{b}] | {cnt:>12}    | {acc2:>11}{mark}")
        acc = acc2
        if a_m is not None and acc2 > i0:
            break
    before = acc - sdfn_count_at_level(n, m, a_m, b_m)
    resid = i0 - before
    lines.append(f"\nChosen core: [{a_m},{b_m}] ; subtract previous = {before} → residual i0 = {resid}")
    alpha = [(a_m, b_m)]

    # lower levels
    prev_a, prev_b = a_m, b_m
    for j in range(m-1, 1, -1):
        acc = 0
        lines.append(f"\nLevel j={j} candidates containing [{prev_a},{prev_b}]:")
        lines.append(" idx | [a,b] | |SDFN(a,b,j)| | accumulated")
        idx_sel = None; chosen = None
        candidates = [(a,b) for (a,b) in order if a <= prev_a and b >= prev_b]
        for idx, (a,b) in enumerate(candidates, start=1):
            cnt = sdfn_count_at_level(n, j, a, b); acc2 = acc + cnt
            mark = ""
            if idx_sel is None and acc2 > resid:
                idx_sel = idx; chosen = (a,b); mark = "  <-- pick"
            lines.append(f"{idx:>4} | [{a},{b}] | {cnt:>12}    | {acc2:>11}{mark}")
            acc = acc2
            if chosen is not None and acc2 > resid:
                break
        before2 = acc - sdfn_count_at_level(n, j, chosen[0], chosen[1])
        resid = resid - before2
        lines.append(f"Chosen j={j}: [{chosen[0]},{chosen[1]}] ; subtract previous = {before2} → residual i0 = {resid}")
        alpha.append(chosen)
        prev_a, prev_b = chosen

    # rebuild μ
    lv = levels_from_alpha_intervals(alpha, n, m)
    seq = [x/(m-1) for x in lv]
    lines.append("\nFinal DFN levels: " + str(tuple(lv)))
    lines.append("μ-form: " + mu_form_from_levels(tuple(lv), m))
    return seq, alpha, "\n".join(lines)

# --------------- RANK (t-inc) — by cuts ---------------
def rank_tinc_by_cuts_from_levels(levels: List[int], n: int, m: int, index_base: int = 1, show_both: bool = False) -> Tuple[int,str]:
    """
    True Section-6 rank. Returns i0 (0-based), but logs in the chosen base.
    """
    lines = []
    lines.append("=== Inverse Academic Example — t-inc RANK (by cuts) ===")
    lines.append(f"Parameters: n={n}, m={m}, comp=tinc")
    lines.append("Input μ: " + mu_form_from_levels(tuple(levels), m))

    order = intervals_tinc_sorted(n)
    alpha = alpha_intervals_from_levels(levels, n, m)

    i0 = 0
    prev_a, prev_b = None, None
    for idx_level, (a_sel, b_sel) in enumerate(alpha):
        j = m - idx_level
        if prev_a is None:
            candidates = order
            lines.append(f"\nLevel j={j} candidates:")
        else:
            candidates = [(a,b) for (a,b) in order if a <= prev_a and b >= prev_b]
            lines.append(f"\nLevel j={j} candidates containing [{prev_a},{prev_b}]:")
        lines.append(" idx | [a,b] | |SDFN(a,b,j)| | cumulative sum added to i")
        cum = 0
        for k, (a,b) in enumerate(candidates, start=1):
            cnt = sdfn_count_at_level(n, j, a, b)
            if (a,b) == (a_sel,b_sel):
                lines.append(f"{k:>4} | [{a},{b}] | {cnt:>12}    | (stop here)")
                break
            else:
                i0 += cnt; cum += cnt
                lines.append(f"{k:>4} | [{a},{b}] | {cnt:>12}    | +{cum}")
        lines.append(f"Chosen j={j}: [{a_sel},{b_sel}] → partial i0 = {i0}")
        prev_a, prev_b = a_sel, b_sel

    i_disp = i0 if index_base == 0 else i0 + 1
    base_label = "0-based" if index_base == 0 else "1-based"
    lines.append(f"\nResult: pos^-1(A) = i ({base_label}) = {i_disp}")
    if show_both:
        lines.append(f"(also: 0-based={i0} / 1-based={i0+1})")
    return i0, "\n".join(lines)

# ---------------- Benchmarks (two plots) ----------------
@dataclass
class BenchRow: m: int; mean_ms: float
def _now_ms() -> float: return time.time() * 1000.0

def _unrank_engine(engine: str, n: int, m: int, i0: int, comp: str) -> List[float]:
    if engine == "module":
        return dcru.unrank_dfn_cuts(n, m, i0, _effective_comp(comp))
    elif engine == "tinc":
        return unrank_tinc_by_cuts(n, m, i0)[0]
    else:
        raise ValueError("Unknown engine")

def _rank_engine(engine: str, n: int, m: int, levels: List[int], comp: str) -> int:
    """Returns the (0-based) index after ranking the levels'."""
    if engine == "module":
        # dcru.rank_dfn_cuts expects μ or levels depending on your lib chosen;
        # If it expects μ in [0,1], must be converted from levels:
        seq = [lv/float(m-1) for lv in levels]
        return dcru.rank_dfn_cuts(n, m, seq, _effective_comp(comp))
    elif engine == "tinc":
        i0, _ = rank_tinc_by_cuts_from_levels(levels, n, m, index_base=0)
        return i0
    else:
        raise ValueError("Unknown engine")


def benchmark(outdir: Path, n: int = 10, m_from: int = 10, m_to: int = 100, m_step: int = 10,
              trials: int = 100, comp: str = "lex1", engine: str = "module") -> pd.DataFrame:
    rows: List[BenchRow] = []
    rows = []
    for m in range(m_from, m_to + 1, m_step):
        unrank_ms = []
        rank_ms = []
        total = dcru.total_dfns(n, m)
        rng = random.Random(12345 + m)

        for _ in range(trials):
            # Random index selection
            i0 = rng.randrange(total)

            # 1) Time UNRANK
            t0 = _now_ms()
            seq = _unrank_engine(engine, n, m, i0, comp)  # devuelve μ en [0,1]
            t1 = _now_ms()
            unrank_ms.append(t1 - t0)

            # 2) Time RANK
            levels = list(seq_levels(seq, m))  # levels 0..m-1
            t2 = _now_ms()
            _ = _rank_engine(engine, n, m, levels, comp)
            t3 = _now_ms()
            rank_ms.append(t3 - t2)

        rows.append({
            "m": m,
            "mean_unrank_ms": float(np.mean(unrank_ms)),
            "std_unrank_ms": float(np.std(unrank_ms, ddof=1)),
            "mean_rank_ms": float(np.mean(rank_ms)),
            "std_rank_ms": float(np.std(rank_ms, ddof=1)),
        })

    df = pd.DataFrame(rows)
    # Save CSV with both metrics
    df_csv = df.rename(columns={
        "mean_unrank_ms": "unrank t̄m (ms)", "std_unrank_ms": "unrank sm (ms)",
        "mean_rank_ms": "rank t̄m (ms)", "std_rank_ms": "rank sm (ms)"
    })
    (outdir / "bench_times_unrank_rank.csv").write_text(
        df_csv.to_csv(index=False, float_format="%.3f"), encoding="utf-8"
    )

    # -------- (a) Average execution time vs m (both curves) ----------
    plt.figure()
    plt.plot(df["m"], df["mean_unrank_ms"], marker="o", label="unranking")
    plt.plot(df["m"], df["mean_rank_ms"], marker="s", label="ranking")
    plt.xlabel(r"$m$")
    plt.ylabel(r"Average execution time $\bar{t}_m$ (ms)")
    plt.title(r"Average execution time as a function of $m$ (for fixed $n=10$) "
              f"\n(n={n}, trials={trials}, order={comp})")
    plt.legend()
    plt.grid(True); plt.tight_layout()
    plt.savefig(outdir / "fig_times_unrank_rank.png", dpi=160); plt.close()

    # -------- (b) Log–log figure and slope fitting (for UNRANK) -------
    x = df["m"].to_numpy(dtype=float)
    y_u = df["mean_unrank_ms"].to_numpy(dtype=float)
    y_r = df["mean_rank_ms"].to_numpy(dtype=float)
    mask_u = (x > 0) & (y_u > 0) & np.isfinite(y_u)
    mask_r = (x > 0) & (y_r > 0) & np.isfinite(y_r)

    logx_u = np.log(x[mask_u]); logy_u = np.log(y_u[mask_u])
    slope_u, intercept_u = np.polyfit(logx_u, logy_u, 1)

    logx_r = np.log(x[mask_r]); logy_r = np.log(y_r[mask_r])
    slope_r, intercept_r = np.polyfit(logx_r, logy_r, 1)

    plt.figure()

    plt.plot(logx_u, logy_u, marker="o", linestyle="", label=f"unranking (slope={slope_u:.3f})")
    plt.plot(logx_r, logy_r, marker="s", linestyle="", label=f"ranking (slope={slope_r:.3f})")
    # fitting lines
    plt.plot(logx_u, intercept_u + slope_u*logx_u, label="fit unranking")
    plt.plot(logx_r, intercept_r + slope_r*logx_r, label="fit ranking")
    plt.xlabel(r"$\log m$")
    plt.ylabel(r"$\log \bar{t}_m$")
    plt.title(r"Log–log plot of average runtime vs.~$m$ showing fitted slope $\approx 1$")
    plt.legend(); plt.grid(True); plt.tight_layout()
    plt.savefig(outdir / "fig_loglog_unrank_rank.png", dpi=160); plt.close()

    latex_table = df_csv.to_latex(
    index=False,
    escape=False,
    float_format="%.3f",
    caption=r"Average execution time as a function of $m$ (for fixed $n=10$). "
            r"Results over $K=500$ trials; order=t-inc.",
    label="tab:table1",
    column_format="r r r r r"
    )
    (outdir / "table1_results.tex").write_text(latex_table, encoding="utf-8")
    print(f"[Tabla 1 → {outdir/'table1_results.tex'}]")

    return df

# ---------------- Demos ----------------
def demo_unrank_tinc(n=5,m=6,i_disp=50,index_base=1,show_both=False,outdir=Path("out")) -> str:
    # Convert displayed/base index to internal 0-based
    total = dcru.total_dfns(n, m)
    i0 = i_disp if index_base == 0 else i_disp - 1
    if not (0 <= i0 < total):
        raise ValueError(f"Index out of range after base conversion: got i0={i0}, valid 0..{total-1}")
    seq, alpha, log = unrank_tinc_by_cuts(n,m,i0,index_base=index_base,show_both=show_both)
    (outdir/"example_tinc_unrank.txt").write_text(log, encoding="utf-8")
    return log

def demo_rank_tinc(mu=None, levels=None, n=5, m=6, index_base=1, show_both=False, outdir=Path("out")) -> Tuple[int,str]:
    if levels is None:
        if mu is None:
            seq = [1,1,1,0.4,0.2,0.2] if (n,m)==(5,6) else [0]*(n+1)
        else:
            seq = [float(x.strip()) for x in mu.split(",")]
        levels = list(seq_levels(seq, m))
    i0, log = rank_tinc_by_cuts_from_levels(levels, n, m, index_base=index_base, show_both=show_both)
    (outdir/"example_tinc_rank.txt").write_text(log, encoding="utf-8")
    return i0, log


# ---------------- Plotter ----------------
def make_plots_from_csv(
    csv_path: str,
    outdir: str = "out",
    n: int = 10,
    trials: int = 500,
    order_label: str = "t-inc",
    fig_name_times: str = "fig_times_unrank_rank.png",
    fig_name_loglog: str = "fig_loglog_unrank_rank.png.png",
):
    """
    Load a CSV withtimes and plot:
      (a) Average execution time vs m
      (b) Log–log with fitted slope
    Automatically detects if there is ranking and/or unranking (Automatically detects number of lines to plot)
    """
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(csv_path)

    # Name normalization
    cols = {c.lower().strip(): c for c in df.columns}
    # Columns candidates
    col_m = cols.get("m", "m")

    # Only one serie to plot
    one_series_candidates = [
        "mean_ms", "t̄m (ms)", "tm (ms)", "mean_time_ms", "average_ms"
    ]
    # Two series to plot
    unr_mean = next((cols[k] for k in [
        "mean_unrank_ms", "unrank t̄m (ms)", "unrank mean (ms)", "unrank mean_ms"
    ] if k in cols), None)
    unr_std  = next((cols[k] for k in [
        "std_unrank_ms", "unrank sm (ms)", "unrank std (ms)", "unrank std_ms"
    ] if k in cols), None)
    rnk_mean = next((cols[k] for k in [
        "mean_rank_ms", "rank t̄m (ms)", "rank mean (ms)", "rank mean_ms"
    ] if k in cols), None)
    rnk_std  = next((cols[k] for k in [
        "std_rank_ms", "rank sm (ms)", "rank std (ms)", "rank std_ms"
    ] if k in cols), None)

    # If only one series
    single_mean = None
    if not unr_mean and not rnk_mean:
        for k in one_series_candidates:
            if k in cols:
                single_mean = cols[k]; break
        if not single_mean:
            for c in df.columns:
                if "t" in c and "ms" in c and "̄" in c:
                    single_mean = c; break
        if not single_mean:
            raise ValueError("No encuentro columnas de medias en el CSV. "
                             "Esperaba algo como 'mean_ms' o 'mean_unrank_ms' / 'mean_rank_ms'.")

    # -------- (a) Plot time vs m --------
    plt.figure()
    x = df[col_m].to_numpy(dtype=float)

    if unr_mean or rnk_mean:
        # two series
        y_u = df[unr_mean].to_numpy(dtype=float) if unr_mean else None
        y_r = df[rnk_mean].to_numpy(dtype=float) if rnk_mean else None
        if y_u is not None:
            plt.plot(x, y_u, marker="o", label="unranking")
        if y_r is not None:
            plt.plot(x, y_r, marker="s", label="ranking")
        plt.legend()
    else:
        # one serie
        y = df[single_mean].to_numpy(dtype=float)
        plt.plot(x, y, marker="o", label="unranking (single series)")
        plt.legend()

    plt.xlabel(r"$m$")
    plt.ylabel(r"Average execution time $\bar{t}_m$ (ms)")
    plt.title(r"Average execution time as a function of $m$ (for fixed $n=10$) "
              f"(n={n}, trials={trials}, order={order_label})")
    plt.grid(True); plt.tight_layout()
    path_times = out / fig_name_times
    plt.savefig(path_times, dpi=160); plt.close()

    # -------- (b) Plot log–log + fits --------
    plt.figure()

    def _fit_and_plot(series, marker, label_prefix):
        mask = np.isfinite(x) & (x > 0) & np.isfinite(series) & (series > 0)
        logx = np.log(x[mask]); logy = np.log(series[mask])
        slope, intercept = np.polyfit(logx, logy, 1)
        # points
        plt.plot(logx, logy, marker=marker, linestyle="", label=f"{label_prefix} (slope={slope:.3f})")
        # fitting line
        plt.plot(logx, intercept + slope*logx, label=f"fit {label_prefix}")
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
    plt.title(r"Log–log plot of average runtime vs.~$m$ showing fitted slope $\approx 1$")
    plt.legend(); plt.grid(True); plt.tight_layout()
    path_loglog = out / fig_name_loglog
    plt.savefig(path_loglog, dpi=160); plt.close()

    # Save csv with fits
    info_txt = "\n".join([f"{name} slope = {s:.6f} ; intercept = {b:.6f}" for (name,s,b) in lines_info])
    (out / "bench_fit_fromcsv.txt").write_text(info_txt, encoding="utf-8")

    print(f"[OK] Figuras desde CSV →\n  {path_times}\n  {path_loglog}")


# ---------------- CLI ----------------
def main():
    ap = argparse.ArgumentParser(description="Section-6–consistent t-inc rank/unrank + module engine + full bench CLI")
    ap.add_argument("--demo_unrank", action="store_true", help="Run t-inc UNRANK demo (i -> DFN)")
    ap.add_argument("--demo_rank", action="store_true", help="Run t-inc RANK demo (DFN -> i)")
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--m", type=int, default=6)
    ap.add_argument("--i", type=int, default=50, help="Index for unrank demo, interpreted in --index_base")
    ap.add_argument("--mu", type=str, default=None, help="Comma-separated μ for rank demo")
    ap.add_argument("--levels", type=str, default=None, help="Comma-separated integer levels (0..m-1) for rank demo")
    ap.add_argument("--bench", action="store_true", help="Run timing benchmarks")
    ap.add_argument("--engine", type=str, default="module", choices=["module","tinc"],
                    help="Benchmark engine: 'module' (dfn_cuts_rank_unrank) or 'tinc' (Section-6 by-cuts)")
    ap.add_argument("--comp", type=str, default="lex1", choices=["lex1","lex2","scopia","tinc"],
                    help="Module comparator (ignored for engine=tinc). 'scopia' aliases 'lex2'.")
    ap.add_argument("--m_from", type=int, default=10, help="initial m for the benchmark")
    ap.add_argument("--m_to", type=int, default=100, help="final m for the benchmark")
    ap.add_argument("--m_step", type=int, default=10, help="m step")
    ap.add_argument("--trials", type=int, default=100, help="random queries per point")
    ap.add_argument("--outdir", type=str, default="out", help="Output directory")
    ap.add_argument("--index_base", type=int, choices=[0,1], default=1, help="Index base for demos/logs (1 matches the paper)")
    ap.add_argument("--show_both_indices", action="store_true", help="Also display the other index base")
    ap.add_argument("--plots_from_csv", type=str, default=None,
                    help="Ruta a CSV para regenerar figuras sin cronometrar")
    ap.add_argument("--plots_order_label", type=str, default="t-inc",
                    help="Etiqueta del 'order' que aparecerá en el título (p.ej. t-inc/lex1)")

    args = ap.parse_args()

    outdir = ensure_outdir(args.outdir)

    did = False
    if args.demo_unrank:
        did = True
        log = demo_unrank_tinc(n=args.n, m=args.m, i_disp=args.i,
                               index_base=args.index_base, show_both=args.show_both_indices, outdir=outdir)
        print(log)
    if args.demo_rank:
        did = True
        i0, log = demo_rank_tinc(mu=args.mu, levels=args.levels, n=args.n, m=args.m,
                                 index_base=args.index_base, show_both=args.show_both_indices, outdir=outdir)
        print(log)

    if args.bench:
        did = True
        df = benchmark(outdir, n=args.n, m_from=args.m_from, m_to=args.m_to,
                       m_step=args.m_step, trials=args.trials, comp=args.comp, engine=args.engine)
        print("\n== BENCH ==")
        print(df.to_string(index=False))

    if args.plots_from_csv:
        did = True
        make_plots_from_csv(
            csv_path=args.plots_from_csv,
            outdir=outdir,
            n=args.n,
            trials=args.trials,
            order_label=args.plots_order_label,
        )


    if not did:
        print("Nothing selected. Try --demo_unrank, --demo_rank, or --bench.")

if __name__ == "__main__":
    main()
