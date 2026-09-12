
# coding: utf-8
from math import comb
from functools import cmp_to_key, lru_cache

try:
    import engine as _user_engine  # engine.py en el mismo directorio
except ImportError:
    _user_engine = None


def _interval_cmp_lex1(a, b, adj):
    La, Ra = a
    Lb, Rb = b
    A = (La, adj - Ra - 1)
    B = (Lb, adj - Rb - 1)
    return (A > B) - (A < B)

def _interval_cmp_lex2(a, b, adj):
    """Order [a,b] by (b,a), both coordinates ascending."""
    La, Ra = a
    Lb, Rb = b
    A = (adj - Ra - 1, La)
    B = (adj - Rb - 1, Lb)
    return (A > B) - (A < B)

def _interval_cmp_xy(a, b, adj):
    """
    Xu–Yager order on intervals [a,b] embedded in L_n.
    Given (L,R) we interpret the core interval as [L, n-R] = [L, adj-R-1].

    [a,b] <=_XY [c,d]  si:
      (a + b < c + d)  o
      (a + b == c + d y b - a <= d - c)
    """
    La, Ra = a
    Lb, Rb = b
    # Convert (L,R) -> true endpoints [a,b] in {0,...,n}
    a1, b1 = La, adj - Ra - 1
    a2, b2 = Lb, adj - Rb - 1

    s1 = a1 + b1
    s2 = a2 + b2
    if s1 < s2:
        return -1
    if s1 > s2:
        return 1

    # tie-break: shorter interval first
    len1 = b1 - a1
    len2 = b2 - a2
    if len1 < len2:
        return -1
    if len1 > len2:
        return 1
    return 0

def _interval_cmp_engine(a, b, adj):
    """
    Comparator that leaves the intervalar order to engine._engine(I, J).

    engine._engine(I, J) should return True iff [a,b] <= [c,d] in the defined order
    """
    if _user_engine is None or not hasattr(_user_engine, "_engine"):
        raise RuntimeError(
            "interval_order='engine' requires a engine.py with "
            "some function _engine(I, J) defined."
        )

    La, Ra = a
    Lb, Rb = b

    a1, b1 = La, adj - Ra - 1
    a2, b2 = Lb, adj - Rb - 1

    # [a1,b1] <= [a2,b2]?
    le_ab = _user_engine._engine((a1, b1), (a2, b2))
    # [a2,b2] <= [a1,b1]?
    le_ba = _user_engine._engine((a2, b2), (a1, b1))


    if le_ab and not le_ba:
        return -1
    if le_ba and not le_ab:
        return 1
    return 0



@lru_cache(maxsize=None)
def _order_groups(n, m, comparator_name):
    """List all (L,R,C) blocks in ∆-up order by the core interval I_1 (top α-cut)."""
    adj = n + 1
    groups = []
    for L in range(adj):
        for R in range(adj - L):
            C = adj - L - R
            if C < 1:
                continue
            size = comb(L + m - 2, m - 2) * comb(R + m - 2, m - 2) if m >= 2 else 1
            groups.append((L, R, C, size))

    if comparator_name == "lex1":
        def cmp(a, b):
            return _interval_cmp_lex1((a[0], a[1]), (b[0], b[1]), adj)
    elif comparator_name == "lex2":
        def cmp(a, b):
            return _interval_cmp_lex2((a[0], a[1]), (b[0], b[1]), adj)
    elif comparator_name == "xy":
        def cmp(a, b):
            return _interval_cmp_xy((a[0], a[1]), (b[0], b[1]), adj)
    elif comparator_name == "engine":
        def cmp(a, b):
            return _interval_cmp_engine((a[0], a[1]), (b[0], b[1]), adj)
    else:
        raise ValueError(f"Unknown comparator_name={comparator_name!r} in _order_groups")


    return tuple(sorted(groups, key=cmp_to_key(cmp)))


def preprocess_order(n, m, comparator_name):
    """Build and cache order-dependent blocks before amortized query timing."""
    return _order_groups(n, m, comparator_name)


def clear_order_cache():
    """Clear cached interval groups for reproducible cold-start measurements."""
    _order_groups.cache_clear()


def _ordered_containing_intervals(n, left, right, comparator_name):
    """Yield intervals containing [left,right] in the requested order."""
    if comparator_name == "lex1":
        return ((a, b) for a in range(left + 1) for b in range(right, n + 1))
    if comparator_name == "lex2":
        return ((a, b) for b in range(right, n + 1) for a in range(left + 1))
    if comparator_name == "xy":
        def xy_intervals():
            for endpoint_sum in range(right, left + n + 1):
                max_a = min(left, endpoint_sum - right)
                min_a = max(0, endpoint_sum - n)
                for a in range(max_a, min_a - 1, -1):
                    yield a, endpoint_sum - a
        return xy_intervals()
    if comparator_name == "engine":
        adj = n + 1
        intervals = [
            (a, b) for a in range(left + 1) for b in range(right, n + 1)
        ]
        return iter(sorted(
            intervals,
            key=cmp_to_key(lambda first, second: _interval_cmp_engine(
                (first[0], n - first[1]),
                (second[0], n - second[1]),
                adj,
            )),
        ))
    raise ValueError(f"Unknown comparator_name={comparator_name!r}")


def _nested_unrank_in_group(n, L, R, m, rem, comparator_name):
    """Choose all lower cuts in the selected interval order."""
    y_s = m - 1
    if y_s <= 1:
        return [0] * y_s, [0] * y_s
    Lext = [0] * y_s
    Rext = [0] * y_s
    ell_prev = r_prev = 0
    core_right = n - R
    for t in range(y_s - 1, 0, -1):
        remaining = t - 1
        previous_left = L - ell_prev
        previous_right = core_right + r_prev
        candidates = _ordered_containing_intervals(
            n, previous_left, previous_right, comparator_name
        )
        for a, b in candidates:
            count = comb(a + remaining, remaining) * comb(
                (n - b) + remaining, remaining
            )
            if rem >= count:
                rem -= count
                continue
            ell, rr = L - a, b - core_right
            Lext[t], Rext[t] = ell, rr
            ell_prev, r_prev = ell, rr
            break
        else:
            raise RuntimeError("rem desfasado dentro del bloque.")
    return Lext, Rext

def _rebuild_seq_from_ext(n, L, R, C, m, Lext, Rext):
    """Build DFN (sequence of membership values) from chosen extensions per level."""
    y_s = m - 1
    adj = n + 1
    seq = [0.0]*adj
    s = L
    # core
    for i in range(s, s + C):
        seq[i] = 1.0
    # left
    for i in range(0, s):
        d = s - i
        val = 0.0
        for t in range(1, y_s):
            if Lext[t] >= d:
                val = max(val, t / y_s)
        seq[i] = val
    # right
    for i in range(s + C, adj):
        d = i - (s + C - 1)
        val = 0.0
        for t in range(1, y_s):
            if Rext[t] >= d:
                val = max(val, t / y_s)
        seq[i] = val
    return seq

def unrank_dfn_cuts(n, m, idx, comparator_name="lex1"):
    """Global unranking → DFN for index idx using the nested-cuts method only (no heap)."""
    if comparator_name not in ("lex1", "lex2", "xy", "engine"):
        raise ValueError("comparator_name debe ser 'lex1', 'lex2', 'xy' or 'engine'")
    groups_sorted = _order_groups(n, m, comparator_name)
    rem = idx
    chosen = None
    for (L, R, C, size) in groups_sorted:
        if rem < size:
            chosen = (L, R, C, size); break
        rem -= size
    if chosen is None:
        raise IndexError("idx fuera de rango")
    L, R, C, _ = chosen
    Lext, Rext = _nested_unrank_in_group(
        n, L, R, m, rem, comparator_name
    )
    return _rebuild_seq_from_ext(n, L, R, C, m, Lext, Rext)

def _extract_ext_from_seq(seq, n, m):
    """Given a DFN sequence (values in {0,1/(m-1),...,1}), extract L,R,C and the extensions Lext,Rext.
       Propagates level t to all u<=t (nested cuts)."""
    y_s = m - 1
    adj = n + 1
    ones = [i for i, v in enumerate(seq) if abs(v - 1.0) < 1e-12]
    if not ones:
        raise ValueError("DFN should have a core")
    s = min(ones); e = max(ones)
    C = e - s + 1
    L = s
    R = adj - (s + C)
    Lext = [0]*y_s
    Rext = [0]*y_s
    for i in range(s):
        d = s - i
        t = int(round(seq[i]*(y_s)))
        for u in range(1, t+1):
            Lext[u] = max(Lext[u], d)
    for i in range(e+1, adj):
        d = i - e
        t = int(round(seq[i]*(y_s)))
        for u in range(1, t+1):
            Rext[u] = max(Rext[u], d)
    return L, R, C, Lext, Rext

def _rank_in_group_from_ext(n, L, R, m, Lext, Rext, comparator_name):
    """Compute the intra-block rank using the selected order at every cut."""
    y_s = m - 1
    if y_s <= 1:
        return 0
    rem = 0
    ell_prev = r_prev = 0
    core_right = n - R
    for t in range(y_s - 1, 0, -1):
        remaining = t - 1
        target_ell = Lext[t]
        target_rr = Rext[t]
        previous_left = L - ell_prev
        previous_right = core_right + r_prev
        candidates = _ordered_containing_intervals(
            n, previous_left, previous_right, comparator_name
        )
        target = (L - target_ell, core_right + target_rr)
        for a, b in candidates:
            if (a, b) == target:
                ell_prev, r_prev = target_ell, target_rr
                break
            rem += comb(a + remaining, remaining) * comb(
                (n - b) + remaining, remaining
            )
        else:
            raise ValueError(f"cut {target} is absent from the interval order")
    return rem

def rank_dfn_cuts(n, m, seq, comparator_name="lex1"):
    """Global rank (idx) of DFN 'seq' under ∆-up + comparator_name
       (lex1, lex2, xy o engine) usando nested cuts.
    """
    groups = _order_groups(n, m, comparator_name)
    # extract its group and extensions
    L, R, C, Lext, Rext = _extract_ext_from_seq(seq, n, m)
    # add sizes of previous groups
    idx = 0
    for (Lg, Rg, Cg, size) in groups:
        if (Lg, Rg, Cg) == (L, R, C):
            break
        idx += size
    # intra-block offset
    idx += _rank_in_group_from_ext(
        n, L, R, m, Lext, Rext, comparator_name
    )
    return idx

def total_dfns(n, m):
    """Total number of DFNs in AL_n×Y_m (support interval)"""
    adj = n + 1
    total = 0
    for L in range(adj):
        for R in range(adj - L):
            C = adj - L - R
            if C < 1:
                continue
            size = comb(L + m - 2, m - 2) * comb(R + m - 2, m - 2) if m >= 2 else 1
            total += size
    return total
