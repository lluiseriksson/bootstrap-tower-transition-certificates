"""Generate all three 300-DPI publication figures for Paper 32:
BOOTSTRAP-TOWER-CERT-1.0 (Koval & Peres, arXiv:2611.33087).

Author: Lluis Eriksson
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.grid": True,
    "grid.alpha": 0.28,
})


def make_fig1() -> None:
    d1 = json.loads((DATA_DIR / "exp1_tower_curve_grid.json").read_text(encoding="utf-8"))
    rhos = [float(r["rho_dec"]) for r in d1["tower_curve"]]
    ps = [float(r["p_rho_dec"]) for r in d1["tower_curve"]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.6), dpi=300)

    # Panel (a): Theorem 1.1 exact tower curve p_rho = (1 - rho) / (2 - rho)
    ax1.plot(rhos, ps, color="#1f4e79", lw=2.2, label=r"Torus curve $p_\rho = (1-\rho)/(2-\rho)$")
    ax1.fill_between(rhos, ps, 0.56, color="#d9e8f5", alpha=0.45, label=r"Supercritical ($p > p_\rho$, fills)")
    ax1.fill_between(rhos, 0.0, ps, color="#fbe4d5", alpha=0.45, label=r"Subcritical ($p < p_\rho$, trapped)")
    # Highlight Morris (2/3, 1/4) and certified Mode 1 points
    pts = [
        (3 / 8, 5 / 13, r"$d=8, n=16$"),
        (2 / 5, 3 / 8, r"$d=10, n=65536$"),
        (2 / 3, 1 / 4, r"Morris $(2/3, 1/4)$"),
        (1.0, 0.0, r"Saturated $(1, 0)$"),
    ]
    for rx, py, lbl in pts:
        color = "#c00000" if "Morris" in lbl else "#2e75b6"
        ax1.scatter([rx], [py], s=48, color=color, zorder=5, edgecolor="white", linewidth=0.8)
    ax1.annotate(
        r"Morris $(\rho=2/3, p=1/4)$",
        xy=(2 / 3, 1 / 4),
        xytext=(0.38, 0.14),
        arrowprops=dict(arrowstyle="->", color="#c00000", lw=1.1),
        fontsize=8,
        color="#c00000",
        fontweight="bold",
    )
    ax1.set_xlim(-0.02, 1.04)
    ax1.set_ylim(-0.02, 0.56)
    ax1.set_xlabel(r"Tower-height ratio $\rho = \min(H(n)/d, 1)$")
    ax1.set_ylabel(r"Critical density $p_\rho$")
    ax1.set_title(r"(a) Exact Torus Tower Transition Curve (Thm. 1.1)")
    ax1.legend(loc="upper right", framealpha=0.92, fontsize=7.4)

    # Panel (b): Sharp Box-vs.-Torus Gap (Eq. 1.2)
    ms = [r["m"] for r in d1["dyadic_box_limits"]]
    box_pcs = [float(r["box_pc_dec"]) for r in d1["dyadic_box_limits"]]
    bars = ax2.bar(ms, box_pcs, width=0.62, color="#2b6cb0", edgecolor="#1a365d", alpha=0.88, label=r"Box limit $1 - 2^{-m}$ ($n/d = 1/m$)")
    ax2.axhline(0.5, color="#c00000", ls="--", lw=1.5, label=r"Torus upper envelope $p \leq 1/2$")
    for b, m_val, y_val in zip(bars, ms, box_pcs):
        if m_val <= 4:
            ax2.text(
                b.get_x() + b.get_width() / 2.0,
                y_val + 0.022,
                f"{(1 << m_val) - 1}/{1 << m_val}",
                ha="center",
                va="bottom",
                fontsize=7.5,
                color="#1a365d",
                fontweight="bold",
            )
    ax2.set_xticks(ms)
    ax2.set_ylim(0.0, 1.12)
    ax2.set_xlabel(r"Inverse box ratio $m = d/n$ ($n/d = 1/m$)")
    ax2.set_ylabel(r"Asymptotic critical density $p_c([n]^d)$")
    ax2.set_title(r"(b) Exact Dyadic Box Limit vs. Torus Envelope (Eq. 1.2)")
    ax2.legend(loc="lower right", framealpha=0.92, fontsize=7.6)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig1_tower_curve_and_box_gap.png", dpi=300)
    plt.close(fig)


def make_fig2() -> None:
    d2 = json.loads((DATA_DIR / "exp2_exact_filling_prob.json").read_text(encoding="utf-8"))
    d4 = json.loads((DATA_DIR / "exp4_eq23_barrier_census.json").read_text(encoding="utf-8"))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.6), dpi=300)

    # Panel (a): Exact polynomial filling probability on Z_3^2 vs [3]^2 (512 configs each)
    t_curve = d2["torus_d2_n3"]["exact_fill_curve"]
    b_curve = d2["box_d2_n3"]["exact_fill_curve"]
    ps = [i / 20.0 for i in range(len(t_curve))]
    p_tor = [float(r["prob_fills_dec"]) for r in t_curve]
    p_box = [float(r["prob_fills_dec"]) for r in b_curve]

    ax1.plot(ps, p_tor, "o-", ms=3.5, lw=2.0, color="#1f4e79", label=r"Torus $\mathbb{Z}_3^2$ (exact $2^9=512$ configs)")
    ax1.plot(ps, p_box, "s--", ms=3.5, lw=2.0, color="#c55a11", label=r"Box $[3]^2$ (exact $2^9=512$ configs)")
    ax1.fill_between(ps, p_box, p_tor, color="#ffd966", alpha=0.35, label="Exact boundary deficit")
    ax1.set_xlim(0.0, 1.0)
    ax1.set_ylim(-0.02, 1.05)
    ax1.set_xlabel(r"Initial infection probability $p$")
    ax1.set_ylabel(r"Exact filling probability $\Pr_p(A_\infty = V)$")
    ax1.set_title(r"(a) Exact Filling Polynomials: $\mathbb{Z}_3^2$ vs. $[3]^2$")
    ax1.legend(loc="lower right", framealpha=0.92, fontsize=7.5)

    # Panel (b): Eq. (2.3) d-closed barriers on Z_4^3 (V = 64)
    cases = d4["cases"]
    x = np.arange(len(cases))
    w = 0.22
    seed_sz = [c["seed_size"] for c in cases]
    fin_sz = [c["final_size"] for c in cases]
    b11_sz = [c["h1_g1_barrier_size"] for c in cases]
    b21_sz = [c["h2_g1_barrier_size"] for c in cases]

    ax2.bar(x - 1.5 * w, seed_sz, width=w, color="#7f7f7f", label=r"Initial seeds $|I|$")
    ax2.bar(x - 0.5 * w, fin_sz, width=w, color="#2e75b6", label=r"Final $|A_\infty|$")
    ax2.bar(x + 0.5 * w, b11_sz, width=w, color="#548235", label=r"$|B(V)|$ ($h{=}1,g{=}1$, 4/5 closed)")
    ax2.bar(x + 1.5 * w, b21_sz, width=w, color="#c55a11", label=r"$|B(V)|$ ($h{=}2,g{=}1$, unclosed)")
    ax2.axhline(64, color="#c00000", ls=":", lw=1.4, label=r"Full torus $|V| = 64$")
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"$|I|={c['seed_size']}$" for c in cases])
    ax2.set_ylim(0, 86)
    ax2.set_xlabel(r"Partial $2\times 2\times 2$ subcube seeds in $\mathbb{Z}_4^3$")
    ax2.set_ylabel(r"Vertex count (out of $64$)")
    ax2.set_title(r"(b) Exact Eq. (2.3) $d$-Closed Barriers on $\mathbb{Z}_4^3$")
    ax2.legend(loc="upper center", ncol=2, framealpha=0.95, fontsize=6.7)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig2_exact_filling_and_eq23_barriers.png", dpi=300)
    plt.close(fig)


def make_fig3() -> None:
    d3 = json.loads((DATA_DIR / "exp3_outward_witness_w0.json").read_text(encoding="utf-8"))
    rep = json.loads((DATA_DIR / "certificate_replay_summary.json").read_text(encoding="utf-8"))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.6), dpi=300)

    # Panel (a): Lemma 2.1 exact-Q w_0 decay across d = 10..32
    palette = {"1/4": ("#1f4e79", "o-"), "1/5": ("#548235", "s-"), "1/6": ("#c55a11", "^-")}
    for p_str, (col, fmt) in palette.items():
        sub = [r for r in d3["sweep"] if r["p"] == p_str]
        ds = [r["d"] for r in sub]
        ws = [float(r["w0_dec"]) for r in sub]
        ax1.semilogy(ds, ws, fmt, ms=4.2, lw=1.8, color=col, label=f"$p = {p_str}$ ($H=2,h=1,g=3$)")
    ax1.axhline(0.5, color="#888888", ls="--", lw=1.1, label=r"$w_0 = 1/2$")
    ax1.set_xlim(9.5, 32.8)
    ax1.set_ylim(3e-6, 2.0)
    ax1.set_xlabel(r"Dimension $d$")
    ax1.set_ylabel(r"Exact rational $w_0$ (log scale)")
    ax1.set_title(r"(a) Exact-$\mathbb{Q}$ Outward-Witness Tail $w_0$ (Lemma 2.1)")
    ax1.legend(loc="lower left", framealpha=0.92, fontsize=7.5)

    # Panel (b): Replay latency across all 16 canonical certificates
    certs = rep["certificates"]
    mode_colors = {
        "tower_curve_bijection_exact": "#1f4e79",
        "bootstrap_dynamics_closure_exact": "#2e75b6",
        "outward_witness_geodesic_exact": "#548235",
        "box_dyadic_and_rect_merge_exact": "#c55a11",
    }
    mode_short = {
        "tower_curve_bijection_exact": "Mode 1 (Tower)",
        "bootstrap_dynamics_closure_exact": "Mode 2 (Closure)",
        "outward_witness_geodesic_exact": "Mode 3 (Witness)",
        "box_dyadic_and_rect_merge_exact": "Mode 4 (Box/Rect)",
    }
    ys = np.arange(len(certs))
    lat = [max(c["verify_ms"], 0.005) for c in certs]
    cols = [mode_colors[c["mode"]] for c in certs]
    ax2.barh(ys, lat, color=cols, height=0.72, alpha=0.9)
    ax2.set_xscale("log")
    ax2.set_yticks(ys)
    ax2.set_yticklabels([f"C{i+1:02d}" for i in range(len(certs))], fontsize=7.2)
    ax2.invert_yaxis()
    ax2.set_xlim(0.004, max(lat) * 4.5)
    ax2.set_xlabel(r"Zero-float exact replay latency (ms, log scale)")
    ax2.set_title(r"(b) Replay Latency Across All 16 Certificates")
    # Custom legend
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=col, label=mode_short[m])
        for m, col in mode_colors.items()
    ]
    ax2.legend(handles=handles, loc="lower right", framealpha=0.92, fontsize=7.0)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig3_outward_witness_and_replay_latency.png", dpi=300)
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    make_fig1()
    make_fig2()
    make_fig3()
    print("Generated all 3 publication figures successfully.")


if __name__ == "__main__":
    main()
