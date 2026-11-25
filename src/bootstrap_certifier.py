"""Pure-stdlib zero-float verifier for BOOTSTRAP-TOWER-CERT-1.0 certificates.
Supports four exact-Q verification modes aligned with Vyacheslav Koval & Yuval Peres
(arXiv:2611.33087, Nov 24, 2026) "The tower transition in majority bootstrap percolation":
  1. tower_curve_bijection_exact        (Thm 1.1 & Cor 1.2: exact p_rho <-> rho_p, support balance v+(2d-v)p=d, H(n)=log*_2 n)
  2. bootstrap_dynamics_closure_exact   (Sec 1 & Eq. 2.3: exact Z_n^d / [n]^d filling OR d-closed barrier E_h subset B(V) -> A_inf subset B(V))
  3. outward_witness_geodesic_exact     (Lemma 2.1 Eq. 2.4: exact rational binomial tail w_0, g > 2h, geodesic overlap 2(i+1))
  4. box_dyadic_and_rect_merge_exact    (Eq. 1.2 & Lemma 2.3: exact dyadic box limit 1-2^{-m} + separated-rectangle l_inf > 2t audit)

Author: Lluis Eriksson
"""
from __future__ import annotations

import re
from fractions import Fraction
from typing import Any, Dict, List, Set, Tuple

from src.bootstrap_exact import (
    Coord,
    compute_eq23_barrier_exact,
    exact_binom_tail_q,
    rect_linf_distance,
    run_bootstrap_percolation_exact,
    tower_curve_p_from_rho,
    tower_curve_rho_from_p,
    tower_height_H_exact,
)

SCHEMA_VERSION = "BOOTSTRAP-TOWER-CERT-1.0"
_CANON_INT_RE = re.compile(r"\A(?:0|-?[1-9][0-9]*)\Z")
_CANON_FRAC_RE = re.compile(r"\A-?[1-9][0-9]*/[1-9][0-9]*\Z")


class BootstrapTowerCertificateError(ValueError):
    """Raised whenever any BOOTSTRAP-TOWER-CERT-1.0 verification check fails."""


def _assert_strict_json_types(obj: Any, path: str = "root") -> None:
    t = type(obj)
    if t is int or t is str or obj is None:
        return
    if t is dict:
        for k, v in obj.items():
            if type(k) is not str:
                raise BootstrapTowerCertificateError(f"Non-str dict key at {path}: {type(k).__name__}")
            _assert_strict_json_types(v, f"{path}.{k}")
        return
    if t is list:
        for i, v in enumerate(obj):
            _assert_strict_json_types(v, f"{path}[{i}]")
        return
    raise BootstrapTowerCertificateError(f"Forbidden type {t.__name__} at {path}")


def _parse_int(val: Any, name: str) -> int:
    if type(val) is not int:
        raise BootstrapTowerCertificateError(f"{name} must be strict int, got {type(val).__name__}")
    return val


def _parse_rat(val: Any, name: str) -> Fraction:
    if type(val) is int:
        return Fraction(val, 1)
    if type(val) is str:
        if _CANON_INT_RE.match(val):
            return Fraction(int(val), 1)
        if _CANON_FRAC_RE.match(val):
            num_s, den_s = val.split("/")
            num, den = int(num_s), int(den_s)
            if den < 2:
                raise BootstrapTowerCertificateError(f"{name} denominator must be >= 2: {val!r}")
            q = Fraction(num, den)
            if q.numerator != num or q.denominator != den:
                raise BootstrapTowerCertificateError(f"{name} unreduced fraction: {val!r}")
            return q
        raise BootstrapTowerCertificateError(f"{name} non-canonical rational string: {val!r}")
    raise BootstrapTowerCertificateError(f"{name} invalid type {type(val).__name__}")


def _verify_mode1_tower_bijection(cert: Dict[str, Any]) -> Dict[str, Any]:
    req = {
        "schema", "mode", "name", "dim_d", "side_n",
        "claimed_H_n", "claimed_rho", "claimed_p_rho", "claimed_forced_v",
    }
    if set(cert.keys()) != req:
        raise BootstrapTowerCertificateError(f"Mode 1 key mismatch: {sorted(set(cert.keys()) ^ req)}")
    d = _parse_int(cert["dim_d"], "dim_d")
    n = _parse_int(cert["side_n"], "side_n")
    if d < 2 or d > 1000 or n < 3 or n > 65536:
        raise BootstrapTowerCertificateError("Require 2 <= dim_d <= 1000 and 3 <= side_n <= 65536")

    H_n = tower_height_H_exact(n)
    rho = min(Fraction(H_n, d), Fraction(1))
    p_rho = tower_curve_p_from_rho(rho)
    # Verify exact inverse rho_p == rho when rho < 1 (and p_rho == 0 -> rho_p == 1 when rho == 1)
    if tower_curve_rho_from_p(p_rho) != rho:
        raise BootstrapTowerCertificateError("Internal tower curve inverse check failed")
    # Verify page 2 exact expected support balance: v + (2d - v)*p_rho == d at v = rho * d
    v_forced = rho * d
    if v_forced + (2 * d - v_forced) * p_rho != Fraction(d, 1):
        raise BootstrapTowerCertificateError("Page 2 expected support balance v + (2d - v)*p == d failed")

    if _parse_int(cert["claimed_H_n"], "claimed_H_n") != H_n:
        raise BootstrapTowerCertificateError(f"claimed_H_n mismatch: exact={H_n}")
    if _parse_rat(cert["claimed_rho"], "claimed_rho") != rho:
        raise BootstrapTowerCertificateError(f"claimed_rho mismatch: exact={rho}")
    if _parse_rat(cert["claimed_p_rho"], "claimed_p_rho") != p_rho:
        raise BootstrapTowerCertificateError(f"claimed_p_rho mismatch: exact={p_rho}")
    if _parse_rat(cert["claimed_forced_v"], "claimed_forced_v") != v_forced:
        raise BootstrapTowerCertificateError(f"claimed_forced_v mismatch: exact={v_forced}")

    return {
        "verified": True,
        "mode": cert["mode"],
        "d": d,
        "n": n,
        "H_n": H_n,
        "rho": str(rho),
        "p_rho": str(p_rho),
    }


def _verify_mode2_dynamics_closure(cert: Dict[str, Any]) -> Dict[str, Any]:
    req = {
        "schema", "mode", "name", "dim_d", "side_n", "topology", "seeds",
        "outcome", "claimed_final_size", "claimed_rounds",
        "barrier_h", "barrier_g", "claimed_barrier_size",
    }
    if set(cert.keys()) != req:
        raise BootstrapTowerCertificateError(f"Mode 2 key mismatch: {sorted(set(cert.keys()) ^ req)}")
    d = _parse_int(cert["dim_d"], "dim_d")
    n = _parse_int(cert["side_n"], "side_n")
    topo = cert["topology"]
    if topo not in ("torus", "box"):
        raise BootstrapTowerCertificateError(f"Invalid topology: {topo!r}")
    if not (2 <= d <= 4) or (topo == "torus" and n < 3) or n < 2 or (n ** d) > 256:
        raise BootstrapTowerCertificateError("Invalid (d, n, topology) or graph size > 256")

    seeds_raw = cert["seeds"]
    if type(seeds_raw) is not list or len(seeds_raw) == 0 or len(seeds_raw) > n ** d:
        raise BootstrapTowerCertificateError("seeds must be non-empty list")
    seeds: Set[Coord] = set()
    for i, s in enumerate(seeds_raw):
        if type(s) is not list or len(s) != d:
            raise BootstrapTowerCertificateError(f"seeds[{i}] must be length-{d} coord")
        pt = tuple(_parse_int(c, f"seeds[{i}][{j}]") for j, c in enumerate(s))
        if any(not (0 <= c < n) for c in pt) or pt in seeds:
            raise BootstrapTowerCertificateError(f"seeds[{i}] out of range or duplicate")
        seeds.add(pt)

    outcome = cert["outcome"]
    if outcome not in ("fills", "d_closed_barrier"):
        raise BootstrapTowerCertificateError(f"Invalid outcome: {outcome!r}")

    final_set, rounds = run_bootstrap_percolation_exact(d, n, topo, list(seeds))
    if _parse_int(cert["claimed_final_size"], "claimed_final_size") != len(final_set):
        raise BootstrapTowerCertificateError(f"claimed_final_size mismatch: exact={len(final_set)}")
    if _parse_int(cert["claimed_rounds"], "claimed_rounds") != rounds:
        raise BootstrapTowerCertificateError(f"claimed_rounds mismatch: exact={rounds}")

    if outcome == "fills":
        if cert["barrier_h"] is not None or cert["barrier_g"] is not None or cert["claimed_barrier_size"] is not None:
            raise BootstrapTowerCertificateError("fills outcome requires null barrier fields")
        if len(final_set) != n ** d:
            raise BootstrapTowerCertificateError(f"Graph did NOT fill: {len(final_set)} < {n**d}")
        return {
            "verified": True,
            "mode": cert["mode"],
            "outcome": outcome,
            "d": d,
            "n": n,
            "final_size": len(final_set),
            "rounds": rounds,
        }

    # outcome == "d_closed_barrier": verify Eq. (2.3) closure E_h subset B(V) AND A_inf subset B(V) < n^d!
    h = _parse_int(cert["barrier_h"], "barrier_h")
    g = _parse_int(cert["barrier_g"], "barrier_g")
    if h < 1 or g < 1 or h * g >= d:
        raise BootstrapTowerCertificateError("Require h >= 1, g >= 1, and G = h*g <= d - 1")
    _, B, closed = compute_eq23_barrier_exact(d, n, topo, seeds, h, g)
    if not closed:
        raise BootstrapTowerCertificateError("Eq. (2.3) closure E_h(V) subset B(V) failed")
    if len(B) >= n ** d:
        raise BootstrapTowerCertificateError("Barrier B(V) covers entire graph (vacuous barrier)")
    if not final_set.issubset(B):
        raise BootstrapTowerCertificateError("Impossible: A_inf not contained in closed barrier B(V)")
    if _parse_int(cert["claimed_barrier_size"], "claimed_barrier_size") != len(B):
        raise BootstrapTowerCertificateError(f"claimed_barrier_size mismatch: exact={len(B)}")
    return {
        "verified": True,
        "mode": cert["mode"],
        "outcome": outcome,
        "d": d,
        "n": n,
        "final_size": len(final_set),
        "barrier_size": len(B),
        "total_verts": n ** d,
    }


def _verify_mode3_outward_witness(cert: Dict[str, Any]) -> Dict[str, Any]:
    req = {
        "schema", "mode", "name", "dim_d", "p", "forced_H",
        "layers_h", "degree_g", "claimed_u", "claimed_binom_m",
        "claimed_thresh", "claimed_w0",
    }
    if set(cert.keys()) != req:
        raise BootstrapTowerCertificateError(f"Mode 3 key mismatch: {sorted(set(cert.keys()) ^ req)}")
    d = _parse_int(cert["dim_d"], "dim_d")
    p = _parse_rat(cert["p"], "p")
    H = _parse_int(cert["forced_H"], "forced_H")
    h = _parse_int(cert["layers_h"], "layers_h")
    g = _parse_int(cert["degree_g"], "degree_g")
    if not (4 <= d <= 32) or not (Fraction(0) < p < Fraction(1, 2)) or p.denominator > 1000:
        raise BootstrapTowerCertificateError("Require 4 <= d <= 32, 0 < p < 1/2, and p.denominator <= 1000")
    if not (0 <= H < d) or h < 1 or g <= 2 * h:
        raise BootstrapTowerCertificateError("Require 0 <= H < d, h >= 1, and g > 2*h (Lemma 2.1)")

    u = g - 2 * h
    G = h * g
    m = 2 * d - H
    thresh = d - G - 2 * h - H
    if thresh < 1 or thresh > m:
        raise BootstrapTowerCertificateError(f"Invalid binomial threshold {thresh} for m={m}")
    w0 = (Fraction(1) - p) * exact_binom_tail_q(m, p, thresh)
    if not (Fraction(0) < w0 < Fraction(1)):
        raise BootstrapTowerCertificateError(f"w0={w0} not strictly in (0, 1)")

    if _parse_int(cert["claimed_u"], "claimed_u") != u:
        raise BootstrapTowerCertificateError(f"claimed_u mismatch: exact={u}")
    if _parse_int(cert["claimed_binom_m"], "claimed_binom_m") != m:
        raise BootstrapTowerCertificateError(f"claimed_binom_m mismatch: exact={m}")
    if _parse_int(cert["claimed_thresh"], "claimed_thresh") != thresh:
        raise BootstrapTowerCertificateError(f"claimed_thresh mismatch: exact={thresh}")
    if _parse_rat(cert["claimed_w0"], "claimed_w0") != w0:
        raise BootstrapTowerCertificateError(f"claimed_w0 mismatch: exact={w0}")

    return {
        "verified": True,
        "mode": cert["mode"],
        "d": d,
        "p": str(p),
        "u": u,
        "w0_lt_half": w0 < Fraction(1, 2),
        "w0_num_bits": w0.numerator.bit_length(),
    }


def _verify_mode4_box_and_rect_merge(cert: Dict[str, Any]) -> Dict[str, Any]:
    req = {
        "schema", "mode", "name", "box_inv_ratio_m", "claimed_box_pc",
        "dim_d", "side_cap_L", "sep_radius_t",
        "marks", "rectangles", "claimed_min_rect_dist",
    }
    if set(cert.keys()) != req:
        raise BootstrapTowerCertificateError(f"Mode 4 key mismatch: {sorted(set(cert.keys()) ^ req)}")
    m = _parse_int(cert["box_inv_ratio_m"], "box_inv_ratio_m")
    if not (1 <= m <= 12):
        raise BootstrapTowerCertificateError("Require 1 <= box_inv_ratio_m <= 12")
    box_pc = Fraction((1 << m) - 1, 1 << m)  # exact 1 - 2^{-m} (Eq. 1.2 at a = 1/m!)
    if _parse_rat(cert["claimed_box_pc"], "claimed_box_pc") != box_pc:
        raise BootstrapTowerCertificateError(f"claimed_box_pc mismatch: exact={box_pc}")

    d = _parse_int(cert["dim_d"], "dim_d")
    L = _parse_int(cert["side_cap_L"], "side_cap_L")
    t = _parse_int(cert["sep_radius_t"], "sep_radius_t")
    if not (2 <= d <= 4) or not (2 <= L <= 1000) or not (1 <= t <= 500):
        raise BootstrapTowerCertificateError("Invalid d, side_cap_L, or sep_radius_t")

    marks_raw = cert["marks"]
    rects_raw = cert["rectangles"]
    if (
        type(marks_raw) is not list or not (2 <= len(marks_raw) <= 256)
        or type(rects_raw) is not list or not (2 <= len(rects_raw) <= 64)
    ):
        raise BootstrapTowerCertificateError("Require 2..256 marks and 2..64 separated rectangles")

    marks: List[Coord] = []
    seen_marks: Set[Coord] = set()
    for i, mk in enumerate(marks_raw):
        if type(mk) is not list or len(mk) != d:
            raise BootstrapTowerCertificateError(f"marks[{i}] must have length {d}")
        pt = tuple(_parse_int(c, f"marks[{i}][{j}]") for j, c in enumerate(mk))
        if any(abs(c) > 1000000 for c in pt) or pt in seen_marks:
            raise BootstrapTowerCertificateError(f"marks[{i}] duplicate or out of coordinate cap")
        seen_marks.add(pt)
        marks.append(pt)

    rects: List[List[Tuple[int, int]]] = []
    for i, rc in enumerate(rects_raw):
        if type(rc) is not list or len(rc) != d:
            raise BootstrapTowerCertificateError(f"rectangles[{i}] must have {d} [lo, hi] intervals")
        box: List[Tuple[int, int]] = []
        for j, iv in enumerate(rc):
            if type(iv) is not list or len(iv) != 2:
                raise BootstrapTowerCertificateError(f"rect[{i}][{j}] must be [lo, hi]")
            lo, hi = _parse_int(iv[0], f"r[{i}][{j}][0]"), _parse_int(iv[1], f"r[{i}][{j}][1]")
            if abs(lo) > 1000000 or abs(hi) > 1000000 or lo > hi or (hi - lo + 1) >= L:
                raise BootstrapTowerCertificateError(f"rect[{i}][{j}] invalid interval or side >= L={L}")
            box.append((lo, hi))
        rects.append(box)

    # 1. Verify every mark lies in at least one rectangle AND every rectangle contains >= 1 mark
    rect_hit = [0] * len(rects)
    for idx_m, mk in enumerate(marks):
        covered = False
        for idx_r, box in enumerate(rects):
            if all(lo <= mk[j] <= hi for j, (lo, hi) in enumerate(box)):
                covered = True
                rect_hit[idx_r] += 1
        if not covered:
            raise BootstrapTowerCertificateError(f"marks[{idx_m}]={mk} not covered by any rectangle")
    if any(h == 0 for h in rect_hit):
        raise BootstrapTowerCertificateError("Empty rectangle containing no marks")

    # 2. Verify pairwise l_inf separation > 2t (Lemma 2.3!)
    min_dist = None
    for i in range(len(rects)):
        for j in range(i + 1, len(rects)):
            dist = rect_linf_distance(rects[i], rects[j])
            if dist <= 2 * t:
                raise BootstrapTowerCertificateError(f"Rectangles {i} and {j} violate l_inf > 2t={2*t}: dist={dist}")
            if min_dist is None or dist < min_dist:
                min_dist = dist
    assert min_dist is not None
    if _parse_int(cert["claimed_min_rect_dist"], "claimed_min_rect_dist") != min_dist:
        raise BootstrapTowerCertificateError(f"claimed_min_rect_dist mismatch: exact={min_dist}")

    return {
        "verified": True,
        "mode": cert["mode"],
        "box_pc": str(box_pc),
        "num_rects": len(rects),
        "min_rect_dist": min_dist,
    }


def verify_bootstrap_tower_certificate_exact(cert: Dict[str, Any]) -> Dict[str, Any]:
    if type(cert) is not dict:
        raise BootstrapTowerCertificateError("Certificate root must be a dict")
    _assert_strict_json_types(cert)
    if cert.get("schema") != SCHEMA_VERSION:
        raise BootstrapTowerCertificateError(f"Invalid schema: {cert.get('schema')!r}")
    name = cert.get("name")
    if type(name) is not str or not name.strip():
        raise BootstrapTowerCertificateError(f"Certificate name must be non-empty str, got {name!r}")
    mode = cert.get("mode")
    if mode == "tower_curve_bijection_exact":
        return _verify_mode1_tower_bijection(cert)
    if mode == "bootstrap_dynamics_closure_exact":
        return _verify_mode2_dynamics_closure(cert)
    if mode == "outward_witness_geodesic_exact":
        return _verify_mode3_outward_witness(cert)
    if mode == "box_dyadic_and_rect_merge_exact":
        return _verify_mode4_box_and_rect_merge(cert)
    raise BootstrapTowerCertificateError(f"Unknown mode: {mode!r}")
