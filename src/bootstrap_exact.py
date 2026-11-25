"""Pure-stdlib exact-Q arithmetic & discrete bootstrap percolation engine for Paper 32:
The Tower Transition in Majority Bootstrap Percolation
(Vyacheslav Koval & Yuval Peres, arXiv:2611.33087, Nov 24, 2026).

ZERO floating-point arithmetic anywhere.
Author: Lluis Eriksson
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction
from typing import Dict, List, Sequence, Set, Tuple

Coord = Tuple[int, ...]


def fmt_rat(q: Fraction) -> str:
    if q.denominator == 1:
        return str(q.numerator)
    return f"{q.numerator}/{q.denominator}"


# =============================================================================
# 1. Exact Tower Function & Tower-Curve Bijection (Theorem 1.1 & Cor. 1.2)
# =============================================================================

TOWER_TABLE = [1, 2, 4, 16, 65536]  # T_0..T_4; T_5 = 2^65536


def tower_height_H_exact(n: int) -> int:
    """Exact log*_2(n) = min {k >= 0 : n <= T_k} for n >= 1 (exact up to T_4 = 65536; >65536 -> 5)."""
    if n < 1:
        raise ValueError("Require n >= 1")
    for k, tk in enumerate(TOWER_TABLE):
        if n <= tk:
            return k
    # For n > 65536, T_5 = 2^65536 has 19729 decimal digits
    if n.bit_length() <= 65537 and n <= (1 << 65536):
        return 5
    raise ValueError("n exceeds T_5")


def tower_curve_p_from_rho(rho: Fraction) -> Fraction:
    """p_rho = (1 - rho) / (2 - rho) for rho in [0, 1]."""
    return (Fraction(1) - rho) / (Fraction(2) - rho)


def tower_curve_rho_from_p(p: Fraction) -> Fraction:
    """rho_p = (1 - 2p) / (1 - p) for p in [0, 1/2]."""
    return (Fraction(1) - 2 * p) / (Fraction(1) - p)


# =============================================================================
# 2. Exact Finite Bootstrap Percolation & Eq. (2.3) Barrier Engine
# =============================================================================

def graph_neighbors(coord: Coord, n: int, topology: str) -> List[Coord]:
    """Nearest neighbours on torus Z_n^d (n >= 3) or box [n]^d (0..n-1)."""
    d = len(coord)
    nbrs: List[Coord] = []
    for axis in range(d):
        val = coord[axis]
        if topology == "torus":
            for delta in (-1, 1):
                nxt = list(coord)
                nxt[axis] = (val + delta) % n
                nbrs.append(tuple(nxt))
        elif topology == "box":
            if val > 0:
                nxt = list(coord)
                nxt[axis] = val - 1
                nbrs.append(tuple(nxt))
            if val + 1 < n:
                nxt = list(coord)
                nxt[axis] = val + 1
                nbrs.append(tuple(nxt))
        else:
            raise ValueError(f"Unknown topology: {topology!r}")
    return nbrs


def run_bootstrap_percolation_exact(
    d: int, n: int, topology: str, seeds: Sequence[Coord], threshold_r: int | None = None
) -> Tuple[Set[Coord], int]:
    """Run exact deterministic r-neighbour bootstrap percolation (default r = d).
    Returns (final_infected_set A_inf, num_rounds)."""
    r = d if threshold_r is None else threshold_r
    all_verts = list(itertools.product(range(n), repeat=d))
    infected: Set[Coord] = set(seeds)
    rounds = 0
    while True:
        newly: List[Coord] = []
        for v in all_verts:
            if v in infected:
                continue
            cnt = sum(1 for nb in graph_neighbors(v, n, topology) if nb in infected)
            if cnt >= r:
                newly.append(v)
        if not newly:
            break
        infected.update(newly)
        rounds += 1
    return infected, rounds


def compute_eq23_barrier_exact(
    d: int, n: int, topology: str, init_I: Set[Coord], h: int, g: int
) -> Tuple[List[Set[Coord]], Set[Coord], bool]:
    """Exact Eq. (2.3) layers E_0(V)..E_h(V) and barrier B(V) = I U E_0 U ... U E_{h-1}.
    Returns (layers [E_0..E_h], B(V), is_closed_Eh_subset_B)."""
    G = h * g
    all_verts = list(itertools.product(range(n), repeat=d))
    E0 = {
        v for v in all_verts
        if v not in init_I and sum(1 for nb in graph_neighbors(v, n, topology) if nb in init_I) >= d - G
    }
    layers: List[Set[Coord]] = [E0]
    for _ in range(1, h + 1):
        prev = layers[-1]
        Ei = {
            v for v in all_verts
            if sum(1 for nb in graph_neighbors(v, n, topology) if nb in prev) > g
        }
        layers.append(Ei)
    B = set(init_I)
    for idx in range(h):
        B.update(layers[idx])
    closed = layers[h].issubset(B)
    return layers, B, closed


# =============================================================================
# 3. Exact Rational Binomial Tail w_0 (Lemma 2.1 Eq. 2.4)
# =============================================================================

def exact_binom_tail_q(m: int, p: Fraction, thresh: int) -> Fraction:
    """Exact rational Pr[Bin(m, p) >= thresh] = sum_{k=thresh}^m C(m, k) p^k (1-p)^{m-k}."""
    if thresh <= 0:
        return Fraction(1)
    if thresh > m:
        return Fraction(0)
    q = Fraction(1) - p
    total = Fraction(0)
    for k in range(thresh, m + 1):
        total += Fraction(math.comb(m, k), 1) * (p ** k) * (q ** (m - k))
    return total


# =============================================================================
# 4. Exact Lemma 2.3 Separated-Rectangle l_inf Distance Check
# =============================================================================

def rect_linf_distance(rect1: Sequence[Tuple[int, int]], rect2: Sequence[Tuple[int, int]]) -> int:
    """Exact l_inf distance between two axis-aligned boxes [lo_i, hi_i] in Z^d."""
    best = 0
    for (a1, b1), (a2, b2) in zip(rect1, rect2):
        if b1 < a2:
            gap = a2 - b1
        elif b2 < a1:
            gap = a1 - b2
        else:
            gap = 0
        if gap > best:
            best = gap
    return best
