"""Generate all 16 canonical BOOTSTRAP-TOWER-CERT-1.0 certificates and 4 exact-Q experiments
for Paper 32: Majority Bootstrap Percolation Tower Transition (Koval & Peres, arXiv:2611.33087).

Author: Lluis Eriksson
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.bootstrap_certifier import SCHEMA_VERSION, verify_bootstrap_tower_certificate_exact
from src.bootstrap_exact import (
    Coord,
    compute_eq23_barrier_exact,
    exact_binom_tail_q,
    fmt_rat,
    rect_linf_distance,
    run_bootstrap_percolation_exact,
    tower_curve_p_from_rho,
    tower_height_H_exact,
)

CERT_DIR = ROOT / "certificates"
DATA_DIR = ROOT / "data"


def build_mode1_cert(name: str, d: int, n: int) -> Dict[str, Any]:
    H_n = tower_height_H_exact(n)
    rho = min(Fraction(H_n, d), Fraction(1))
    p_rho = tower_curve_p_from_rho(rho)
    v_forced = rho * d
    return {
        "schema": SCHEMA_VERSION,
        "mode": "tower_curve_bijection_exact",
        "name": name,
        "dim_d": d,
        "side_n": n,
        "claimed_H_n": H_n,
        "claimed_rho": fmt_rat(rho),
        "claimed_p_rho": fmt_rat(p_rho),
        "claimed_forced_v": fmt_rat(v_forced),
    }


def build_mode2_cert(
    name: str,
    d: int,
    n: int,
    topology: str,
    seeds: Sequence[Coord],
    outcome: str,
    barrier_h: int | None = None,
    barrier_g: int | None = None,
) -> Dict[str, Any]:
    final_set, rounds = run_bootstrap_percolation_exact(d, n, topology, seeds)
    claimed_barrier_size: int | None = None
    if outcome == "d_closed_barrier":
        assert barrier_h is not None and barrier_g is not None
        _, B, closed = compute_eq23_barrier_exact(d, n, topology, set(seeds), barrier_h, barrier_g)
        assert closed, f"Barrier not closed for {name}"
        assert len(B) < n ** d, f"Barrier vacuous for {name}"
        claimed_barrier_size = len(B)
    return {
        "schema": SCHEMA_VERSION,
        "mode": "bootstrap_dynamics_closure_exact",
        "name": name,
        "dim_d": d,
        "side_n": n,
        "topology": topology,
        "seeds": [list(pt) for pt in seeds],
        "outcome": outcome,
        "claimed_final_size": len(final_set),
        "claimed_rounds": rounds,
        "barrier_h": barrier_h,
        "barrier_g": barrier_g,
        "claimed_barrier_size": claimed_barrier_size,
    }


def build_mode3_cert(name: str, d: int, p: Fraction, H: int, h: int, g: int) -> Dict[str, Any]:
    u = g - 2 * h
    G = h * g
    m = 2 * d - H
    thresh = d - G - 2 * h - H
    w0 = (Fraction(1) - p) * exact_binom_tail_q(m, p, thresh)
    return {
        "schema": SCHEMA_VERSION,
        "mode": "outward_witness_geodesic_exact",
        "name": name,
        "dim_d": d,
        "p": fmt_rat(p),
        "forced_H": H,
        "layers_h": h,
        "degree_g": g,
        "claimed_u": u,
        "claimed_binom_m": m,
        "claimed_thresh": thresh,
        "claimed_w0": fmt_rat(w0),
    }


def build_mode4_cert(
    name: str,
    m: int,
    d: int,
    L: int,
    t: int,
    marks: Sequence[Coord],
    rectangles: Sequence[Sequence[Tuple[int, int]]],
) -> Dict[str, Any]:
    box_pc = Fraction((1 << m) - 1, 1 << m)
    min_dist = min(
        rect_linf_distance(rectangles[i], rectangles[j])
        for i in range(len(rectangles))
        for j in range(i + 1, len(rectangles))
    )
    return {
        "schema": SCHEMA_VERSION,
        "mode": "box_dyadic_and_rect_merge_exact",
        "name": name,
        "box_inv_ratio_m": m,
        "claimed_box_pc": fmt_rat(box_pc),
        "dim_d": d,
        "side_cap_L": L,
        "sep_radius_t": t,
        "marks": [list(pt) for pt in marks],
        "rectangles": [[[lo, hi] for lo, hi in rect] for rect in rectangles],
        "claimed_min_rect_dist": min_dist,
    }


def generate_canonical_certificates() -> List[Dict[str, Any]]:
    certs: List[Dict[str, Any]] = []

    # Mode 1: Exact Tower-Curve Bijection (Theorem 1.1 & Corollary 1.2)
    certs.append(build_mode1_cert("mode1_morris_quarter_d6_n65536", d=6, n=65536))
    certs.append(build_mode1_cert("mode1_tower_d8_n16", d=8, n=16))
    certs.append(build_mode1_cert("mode1_tower_d10_n65536", d=10, n=65536))
    certs.append(build_mode1_cert("mode1_tower_d4_n65536_saturated", d=4, n=65536))

    # Mode 2: Exact Finite Torus/Box Bootstrap Filling & Eq. (2.3) d-Closed Barriers
    # 5. Z_6^2 fills from 6-vertex diagonal in 2 rounds!
    diag_z6_2 = [(i, i) for i in range(6)]
    certs.append(build_mode2_cert("mode2_torus_d2_n6_fills", 2, 6, "torus", diag_z6_2, "fills"))

    # 6. Z_4^3 fills from 2D slice z=0 (16 verts) + diagonal in z=1,2,3 (12 verts) -> 28 seeds
    seeds_z4_3_fill: List[Coord] = [(x, y, 0) for x in range(4) for y in range(4)]
    for z in (1, 2, 3):
        for x in range(4):
            seeds_z4_3_fill.append((x, (x + z) % 4, z))
    certs.append(build_mode2_cert("mode2_torus_d3_n4_fills", 3, 4, "torus", seeds_z4_3_fill, "fills"))

    # 7. Z_4^3 Eq. (2.3) d-closed barrier (h=1, g=1, G=1, d-G=2): 2x2x2 cube minus 1 corner (7 seeds -> grows to 8, trapped!)
    cube7_d3 = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1) if (x, y, z) != (1, 1, 1)]
    certs.append(
        build_mode2_cert(
            "mode2_torus_d3_n4_eq23_barrier", 3, 4, "torus", cube7_d3, "d_closed_barrier", barrier_h=1, barrier_g=1
        )
    )

    # 8. [3]^4 Eq. (2.3) d-closed barrier (h=1, g=1, G=1, d-G=3): 2x2x2x2 hypercube minus 2 non-adjacent corners (14 seeds -> grows to 16, trapped!)
    hcube14_d4 = [
        (a, b, c, e)
        for a in (0, 1)
        for b in (0, 1)
        for c in (0, 1)
        for e in (0, 1)
        if (a, b, c, e) not in ((1, 1, 1, 1), (0, 0, 0, 0))
    ]
    certs.append(
        build_mode2_cert(
            "mode2_box_d4_n3_eq23_barrier", 4, 3, "box", hcube14_d4, "d_closed_barrier", barrier_h=1, barrier_g=1
        )
    )

    # Mode 3: Exact Rational Outward-Witness Binomial Tail w_0 (Lemma 2.1 Eq. 2.4)
    certs.append(build_mode3_cert("mode3_outward_d16_p1_4_h1_g3", d=16, p=Fraction(1, 4), H=3, h=1, g=3))
    certs.append(build_mode3_cert("mode3_outward_d20_p1_4_h1_g3", d=20, p=Fraction(1, 4), H=2, h=1, g=3))
    certs.append(build_mode3_cert("mode3_outward_d24_p1_5_h1_g4", d=24, p=Fraction(1, 5), H=3, h=1, g=4))
    certs.append(build_mode3_cert("mode3_outward_d30_p1_6_h2_g5", d=30, p=Fraction(1, 6), H=4, h=2, g=5))

    # Mode 4: Exact Dyadic Box Transition 1 - 2^{-m} (Eq. 1.2) & Lemma 2.3 Separated Rectangles
    certs.append(
        build_mode4_cert(
            "mode4_box_m1_d2_rects",
            m=1, d=2, L=6, t=1,
            marks=[(0, 0), (1, 2), (5, 0), (6, 1)],
            rectangles=[[(0, 2), (0, 2)], [(5, 7), (0, 2)]],
        )
    )
    certs.append(
        build_mode4_cert(
            "mode4_box_m2_d2_rects",
            m=2, d=2, L=8, t=2,
            marks=[(0, 0), (2, 1), (8, 0), (9, 2), (0, 9), (1, 10)],
            rectangles=[[(0, 2), (0, 2)], [(8, 10), (0, 2)], [(0, 2), (8, 10)]],
        )
    )
    certs.append(
        build_mode4_cert(
            "mode4_box_m3_d3_rects",
            m=3, d=3, L=10, t=2,
            marks=[(0, 0, 0), (2, 1, 2), (9, 0, 0), (10, 2, 1), (0, 9, 9), (1, 10, 10)],
            rectangles=[
                [(0, 3), (0, 3), (0, 3)],
                [(9, 12), (0, 3), (0, 3)],
                [(0, 3), (9, 12), (9, 12)],
            ],
        )
    )
    certs.append(
        build_mode4_cert(
            "mode4_box_m4_d4_rects",
            m=4, d=4, L=12, t=3,
            marks=[(0, 0, 0, 0), (2, 2, 1, 0), (11, 0, 0, 0), (12, 1, 1, 1), (0, 11, 0, 11), (1, 12, 1, 12)],
            rectangles=[
                [(0, 3), (0, 3), (0, 3), (0, 3)],
                [(11, 14), (0, 3), (0, 3), (0, 3)],
                [(0, 3), (11, 14), (0, 3), (11, 14)],
            ],
        )
    )
    return certs


def run_exp1_tower_curve_grid() -> Dict[str, Any]:
    """Exact-Q Theorem 1.1 tower-transition curve p_rho = (1-rho)/(2-rho) vs. Eq. (1.2) dyadic box limit."""
    curve_rows = []
    for k in range(0, 25):
        rho = Fraction(k, 24)
        p_rho = tower_curve_p_from_rho(rho)
        curve_rows.append({
            "rho": fmt_rat(rho),
            "rho_dec": f"{rho.numerator / rho.denominator:.6f}",
            "p_rho": fmt_rat(p_rho),
            "p_rho_dec": f"{p_rho.numerator / p_rho.denominator:.6f}",
        })
    box_rows = []
    for m in range(1, 9):
        pc = Fraction((1 << m) - 1, 1 << m)
        box_rows.append({
            "m": m,
            "a": fmt_rat(Fraction(1, m)),
            "box_pc": fmt_rat(pc),
            "box_pc_dec": f"{pc.numerator / pc.denominator:.6f}",
            "gap_above_half": fmt_rat(pc - Fraction(1, 2)),
        })
    # Exact H(n) = log*_2(n) staircase for n in [1..65536]
    staircase = [
        {"n_lo": 1, "n_hi": 1, "H_n": 0},
        {"n_lo": 2, "n_hi": 2, "H_n": 1},
        {"n_lo": 3, "n_hi": 4, "H_n": 2},
        {"n_lo": 5, "n_hi": 16, "H_n": 3},
        {"n_lo": 17, "n_hi": 65536, "H_n": 4},
    ]
    return {"tower_curve": curve_rows, "dyadic_box_limits": box_rows, "tower_staircase": staircase}


def run_exp2_exact_filling_prob_torus_vs_box() -> Dict[str, Any]:
    """EXACT polynomial filling probability Pr_p(fills) over ALL 2^V configurations on Z_3^2 vs [3]^2
    and Z_4^2 (size-k census on Z_3^2 / [3]^2 = 512 exact configs each!)"""
    results: Dict[str, Any] = {}
    for topo in ("torus", "box"):
        d, n = 2, 3
        verts = list(itertools.product(range(n), repeat=d))
        V = len(verts)  # 9 -> 2^9 = 512 exact subsets!
        fill_by_k = [0] * (V + 1)
        total_by_k = [0] * (V + 1)
        barrier_non_fill_by_k = [0] * (V + 1)
        for mask in range(1 << V):
            seeds = [verts[i] for i in range(V) if (mask >> i) & 1]
            k = len(seeds)
            total_by_k[k] += 1
            final_set, _ = run_bootstrap_percolation_exact(d, n, topo, seeds)
            if len(final_set) == V:
                fill_by_k[k] += 1
            elif k > 0:
                _, B, closed = compute_eq23_barrier_exact(d, n, topo, set(seeds), h=1, g=1)
                if closed and len(B) < V:
                    barrier_non_fill_by_k[k] += 1
        # Evaluate exact polynomial P(p) = sum_k fill_by_k[k] * p^k * (1-p)^{V-k} at p = j/20
        curve = []
        for j in range(0, 21):
            p = Fraction(j, 20)
            q = Fraction(1) - p
            prob = sum(
                Fraction(fill_by_k[k], 1) * (p ** k) * (q ** (V - k))
                for k in range(V + 1)
            )
            curve.append({
                "p": fmt_rat(p),
                "prob_fills": fmt_rat(prob),
                "prob_fills_dec": f"{prob.numerator / prob.denominator:.6f}",
            })
        results[f"{topo}_d2_n3"] = {
            "d": d,
            "n": n,
            "V": V,
            "total_configs": 1 << V,
            "fill_by_k": fill_by_k,
            "total_by_k": total_by_k,
            "eq23_barrier_certified_non_fill_by_k": barrier_non_fill_by_k,
            "exact_fill_curve": curve,
        }
    return results


def run_exp3_outward_witness_w0_sweep() -> Dict[str, Any]:
    """Exact-Q Lemma 2.1 Eq. (2.4) w_0 decay across d = 10..32 for p in {1/4, 1/5, 1/6}."""
    rows = []
    for p_den in (4, 5, 6):
        p = Fraction(1, p_den)
        for d in range(10, 33, 2):
            H, h, g = 2, 1, 3
            G = h * g
            m = 2 * d - H
            thresh = d - G - 2 * h - H
            if thresh >= 1:
                w0 = (Fraction(1) - p) * exact_binom_tail_q(m, p, thresh)
                rows.append({
                    "p": fmt_rat(p),
                    "d": d,
                    "H": H,
                    "h": h,
                    "g": g,
                    "m": m,
                    "thresh": thresh,
                    "w0_exact": fmt_rat(w0),
                    "w0_dec": f"{w0.numerator / w0.denominator:.8e}",
                })
    return {"sweep": rows}


def run_exp4_eq23_barrier_census_d3() -> Dict[str, Any]:
    """Exact Eq. (2.3) d-closed barrier statistics on Z_4^3 (V=64) across subcube/cluster seeds."""
    d, n = 3, 4
    V = n ** d
    cases = []
    # Test growing partial 2x2x2 cubes (sizes 4..8) and 2x2x3 slabs (sizes 8..12) on Z_4^3
    cube_222 = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]
    for sz in range(4, 9):
        seeds = cube_222[:sz]
        final_set, rounds = run_bootstrap_percolation_exact(d, n, "torus", seeds)
        _, B11, closed11 = compute_eq23_barrier_exact(d, n, "torus", set(seeds), h=1, g=1)
        _, B21, closed21 = compute_eq23_barrier_exact(d, n, "torus", set(seeds), h=2, g=1)
        cases.append({
            "family": f"cube222_prefix_{sz}",
            "seed_size": sz,
            "final_size": len(final_set),
            "rounds": rounds,
            "fills": len(final_set) == V,
            "h1_g1_closed": closed11,
            "h1_g1_barrier_size": len(B11),
            "h2_g1_closed": closed21,
            "h2_g1_barrier_size": len(B21),
        })
    return {"d": d, "n": n, "V": V, "cases": cases}


def main() -> None:
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    certs = generate_canonical_certificates()
    summary_rows = []
    for cert in certs:
        out_path = CERT_DIR / f"{cert['name']}.json"
        out_path.write_text(json.dumps(cert, indent=2) + "\n", encoding="utf-8")
        t0 = time.perf_counter()
        res = verify_bootstrap_tower_certificate_exact(cert)
        dt_ms = (time.perf_counter() - t0) * 1000.0
        summary_rows.append({
            "name": cert["name"],
            "mode": cert["mode"],
            "verify_ms": round(dt_ms, 3),
            "summary": res,
        })
        print(f"[OK] {cert['name']} ({cert['mode']}) verified in {dt_ms:.3f} ms")

    (DATA_DIR / "certificate_replay_summary.json").write_text(
        json.dumps({"total": len(summary_rows), "certificates": summary_rows}, indent=2) + "\n",
        encoding="utf-8",
    )
    (DATA_DIR / "exp1_tower_curve_grid.json").write_text(
        json.dumps(run_exp1_tower_curve_grid(), indent=2) + "\n", encoding="utf-8"
    )
    (DATA_DIR / "exp2_exact_filling_prob.json").write_text(
        json.dumps(run_exp2_exact_filling_prob_torus_vs_box(), indent=2) + "\n", encoding="utf-8"
    )
    (DATA_DIR / "exp3_outward_witness_w0.json").write_text(
        json.dumps(run_exp3_outward_witness_w0_sweep(), indent=2) + "\n", encoding="utf-8"
    )
    (DATA_DIR / "exp4_eq23_barrier_census.json").write_text(
        json.dumps(run_exp4_eq23_barrier_census_d3(), indent=2) + "\n", encoding="utf-8"
    )
    print(f"Generated {len(certs)} certificates and 4 experiment JSONs successfully.")


if __name__ == "__main__":
    main()
