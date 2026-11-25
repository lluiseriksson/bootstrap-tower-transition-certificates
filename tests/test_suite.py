"""Comprehensive pytest suite + 58 live adversarial forgery attacks for Paper 32:
BOOTSTRAP-TOWER-CERT-1.0 (Koval & Peres, arXiv:2611.33087).

Author: Lluis Eriksson
"""
from __future__ import annotations

import copy
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple

import pytest

from src.bootstrap_certifier import (
    BootstrapTowerCertificateError,
    verify_bootstrap_tower_certificate_exact,
)
from src.bootstrap_exact import (
    compute_eq23_barrier_exact,
    exact_binom_tail_q,
    run_bootstrap_percolation_exact,
    tower_curve_p_from_rho,
    tower_curve_rho_from_p,
    tower_height_H_exact,
)

ROOT = Path(__file__).resolve().parent.parent
CERT_DIR = ROOT / "certificates"


def _load_cert(name: str) -> Dict[str, Any]:
    return json.loads((CERT_DIR / f"{name}.json").read_text(encoding="utf-8"))


def test_all_16_canonical_certificates_pass() -> None:
    cert_files = sorted(CERT_DIR.glob("*.json"))
    assert len(cert_files) == 16
    for path in cert_files:
        cert = json.loads(path.read_text(encoding="utf-8"))
        res = verify_bootstrap_tower_certificate_exact(cert)
        assert res["verified"] is True


def test_exact_mathematical_invariants() -> None:
    # 1. Exact tower height H(n) = log*_2(n) boundary transitions
    assert tower_height_H_exact(1) == 0
    assert tower_height_H_exact(2) == 1
    assert tower_height_H_exact(3) == 2
    assert tower_height_H_exact(4) == 2
    assert tower_height_H_exact(5) == 3
    assert tower_height_H_exact(16) == 3
    assert tower_height_H_exact(17) == 4
    assert tower_height_H_exact(65536) == 4

    # 2. Exact tower bijection & Morris benchmark p=1/4 <-> rho=2/3
    assert tower_curve_p_from_rho(Fraction(2, 3)) == Fraction(1, 4)
    assert tower_curve_rho_from_p(Fraction(1, 4)) == Fraction(2, 3)
    for num in range(0, 11):
        rho = Fraction(num, 10)
        p = tower_curve_p_from_rho(rho)
        assert tower_curve_rho_from_p(p) == rho
        for d in (4, 7, 12):
            v = rho * d
            assert v + (2 * d - v) * p == Fraction(d, 1)

    # 3. Exact binomial tail symmetry at p = 1/2: Pr[Bin(2k+1, 1/2) >= k+1] == 1/2
    for k in range(1, 8):
        assert exact_binom_tail_q(2 * k + 1, Fraction(1, 2), k + 1) == Fraction(1, 2)

    # 4. Exact Eq. (2.3) barrier containment on Z_4^3
    seeds = {(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1) if (x, y, z) != (1, 1, 1)}
    final_set, rounds = run_bootstrap_percolation_exact(3, 4, "torus", list(seeds))
    assert len(final_set) == 8 and rounds == 1
    _, B, closed = compute_eq23_barrier_exact(3, 4, "torus", seeds, h=1, g=1)
    assert closed and final_set.issubset(B) and len(B) < 64


def test_adversarial_red_team_forgeries() -> None:
    m1 = _load_cert("mode1_morris_quarter_d6_n65536")
    m2_fill = _load_cert("mode2_torus_d2_n6_fills")
    m2_bar = _load_cert("mode2_torus_d3_n4_eq23_barrier")
    m3 = _load_cert("mode3_outward_d16_p1_4_h1_g3")
    m4 = _load_cert("mode4_box_m2_d2_rects")

    def mut(base: Dict[str, Any], fn: Callable[[Dict[str, Any]], None]) -> Dict[str, Any]:
        c = copy.deepcopy(base)
        fn(c)
        return c

    attacks: List[Tuple[str, Any]] = [
        # --- Global / Schema / Type Confusion (12) ---
        ("root_not_dict", [1, 2, 3]),
        ("wrong_schema", mut(m1, lambda c: c.__setitem__("schema", "BOOTSTRAP-TOWER-CERT-0.9"))),
        ("unknown_mode", mut(m1, lambda c: c.__setitem__("mode", "fake_mode"))),
        ("empty_name", mut(m1, lambda c: c.__setitem__("name", "   "))),
        ("float_injected_m1", mut(m1, lambda c: c.__setitem__("claimed_p_rho", 0.25))),
        ("bool_for_int_d", mut(m1, lambda c: c.__setitem__("dim_d", True))),
        ("extra_key_m1", mut(m1, lambda c: c.__setitem__("extra", 1))),
        ("missing_key_m1", mut(m1, lambda c: c.pop("claimed_H_n"))),
        ("unreduced_frac_m1", mut(m1, lambda c: c.__setitem__("claimed_p_rho", "2/8"))),
        ("den_one_frac_m1", mut(m1, lambda c: c.__setitem__("claimed_forced_v", "4/1"))),
        ("leading_zero_int_str", mut(m1, lambda c: c.__setitem__("claimed_forced_v", "04"))),
        ("whitespace_frac_str", mut(m1, lambda c: c.__setitem__("claimed_p_rho", " 1/4 "))),

        # --- Mode 1: Tower-Curve Bijection Forgeries (10) ---
        ("m1_d_too_small", mut(m1, lambda c: c.__setitem__("dim_d", 1))),
        ("m1_d_too_large", mut(m1, lambda c: c.__setitem__("dim_d", 1001))),
        ("m1_n_too_small", mut(m1, lambda c: c.__setitem__("side_n", 2))),
        ("m1_n_too_large", mut(m1, lambda c: c.__setitem__("side_n", 65537))),
        ("m1_wrong_H_n", mut(m1, lambda c: c.__setitem__("claimed_H_n", 3))),
        ("m1_wrong_rho", mut(m1, lambda c: c.__setitem__("claimed_rho", "1/2"))),
        ("m1_wrong_p_rho", mut(m1, lambda c: c.__setitem__("claimed_p_rho", "1/3"))),
        ("m1_wrong_forced_v", mut(m1, lambda c: c.__setitem__("claimed_forced_v", "3"))),
        ("m1_neg_rho", mut(m1, lambda c: c.__setitem__("claimed_rho", "-2/3"))),
        ("m1_float_n", mut(m1, lambda c: c.__setitem__("side_n", 65536.0))),

        # --- Mode 2: Dynamics & Eq. (2.3) Barrier Forgeries (14) ---
        ("m2_bad_topology", mut(m2_fill, lambda c: c.__setitem__("topology", "sphere"))),
        ("m2_d_out_of_range", mut(m2_fill, lambda c: c.__setitem__("dim_d", 5))),
        ("m2_torus_n_2", mut(m2_fill, lambda c: c.__setitem__("side_n", 2))),
        ("m2_graph_too_big", mut(m2_fill, lambda c: c.__setitem__("side_n", 17))),
        ("m2_empty_seeds", mut(m2_fill, lambda c: c.__setitem__("seeds", []))),
        ("m2_duplicate_seed", mut(m2_fill, lambda c: c["seeds"].append(c["seeds"][0]))),
        ("m2_seed_out_of_range", mut(m2_fill, lambda c: c["seeds"].__setitem__(0, [6, 0]))),
        ("m2_seed_wrong_dim", mut(m2_fill, lambda c: c["seeds"].__setitem__(0, [0, 0, 0]))),
        ("m2_bad_outcome", mut(m2_fill, lambda c: c.__setitem__("outcome", "unknown"))),
        ("m2_wrong_final_size", mut(m2_fill, lambda c: c.__setitem__("claimed_final_size", 35))),
        ("m2_wrong_rounds", mut(m2_fill, lambda c: c.__setitem__("claimed_rounds", c["claimed_rounds"] + 1))),
        ("m2_fills_with_nonnull_barrier", mut(m2_fill, lambda c: c.__setitem__("barrier_h", 1))),
        ("m2_false_fill_claim", mut(m2_bar, lambda c: (c.__setitem__("outcome", "fills"), c.__setitem__("barrier_h", None), c.__setitem__("barrier_g", None), c.__setitem__("claimed_barrier_size", None)))),
        ("m2_barrier_hg_ge_d", mut(m2_bar, lambda c: c.__setitem__("barrier_g", 3))),
        ("m2_barrier_wrong_size", mut(m2_bar, lambda c: c.__setitem__("claimed_barrier_size", 10))),
        ("m2_barrier_not_closed_on_filling_seeds", mut(m2_fill, lambda c: (c.__setitem__("outcome", "d_closed_barrier"), c.__setitem__("barrier_h", 1), c.__setitem__("barrier_g", 1), c.__setitem__("claimed_barrier_size", 20)))),

        # --- Mode 3: Outward-Witness Binomial Tail Forgeries (10) ---
        ("m3_d_too_small", mut(m3, lambda c: c.__setitem__("dim_d", 3))),
        ("m3_p_ge_half", mut(m3, lambda c: c.__setitem__("p", "1/2"))),
        ("m3_p_zero", mut(m3, lambda c: c.__setitem__("p", "0"))),
        ("m3_p_huge_den", mut(m3, lambda c: c.__setitem__("p", "1/1001"))),
        ("m3_g_le_2h", mut(m3, lambda c: c.__setitem__("degree_g", 2))),
        ("m3_H_ge_d", mut(m3, lambda c: c.__setitem__("forced_H", 16))),
        ("m3_wrong_u", mut(m3, lambda c: c.__setitem__("claimed_u", 2))),
        ("m3_wrong_m", mut(m3, lambda c: c.__setitem__("claimed_binom_m", 28))),
        ("m3_wrong_thresh", mut(m3, lambda c: c.__setitem__("claimed_thresh", 7))),
        ("m3_wrong_w0", mut(m3, lambda c: c.__setitem__("claimed_w0", "1/10"))),

        # --- Mode 4: Dyadic Box & Separated-Rectangle Forgeries (12) ---
        ("m4_m_zero", mut(m4, lambda c: c.__setitem__("box_inv_ratio_m", 0))),
        ("m4_wrong_box_pc", mut(m4, lambda c: c.__setitem__("claimed_box_pc", "1/2"))),
        ("m4_rect_side_ge_L", mut(m4, lambda c: c["rectangles"][0].__setitem__(0, [0, 7]))),
        ("m4_rect_inverted_iv", mut(m4, lambda c: c["rectangles"][0].__setitem__(0, [3, 1]))),
        ("m4_uncovered_mark", mut(m4, lambda c: c["marks"].append([5, 5]))),
        ("m4_duplicate_mark", mut(m4, lambda c: c["marks"].append(c["marks"][0]))),
        ("m4_empty_rectangle", mut(m4, lambda c: c["rectangles"].append([[20, 21], [20, 21]]))),
        ("m4_violates_2t_sep", mut(m4, lambda c: c["rectangles"][1].__setitem__(0, [6, 8]))),
        ("m4_wrong_min_dist", mut(m4, lambda c: c.__setitem__("claimed_min_rect_dist", 5))),
        ("m4_too_few_rects", mut(m4, lambda c: c.__setitem__("rectangles", [c["rectangles"][0]]))),
        ("m4_mark_out_of_cap", mut(m4, lambda c: c["marks"].__setitem__(0, [2000000, 0]))),
        ("m4_missing_key", mut(m4, lambda c: c.pop("sep_radius_t"))),
    ]

    assert len(attacks) == 60
    rejected = 0
    for attack_name, payload in attacks:
        with pytest.raises(BootstrapTowerCertificateError, match=".+"):
            verify_bootstrap_tower_certificate_exact(payload)
        rejected += 1
    assert rejected == 60
