# Machine-Checked Exact-Rational Tower Transition, $d$-Closed Barrier, and Outward-Witness Certificates for Majority Bootstrap Percolation

**Author:** Lluis Eriksson (prepared with transparent AI assistance via Google Antigravity)
**Target Reference:** Vyacheslav Koval & Yuval Peres, *The tower transition in majority bootstrap percolation*, `arXiv:2611.33087` (November 24, 2026).
**Licenses:** Manuscript (`paper/main.tex`, `paper/main.pdf`) under `CC-BY-4.0`; code, Lean 4 proofs, and exact certificates under `MIT`.

## Highlights
- **Zero-Dependency Lean 4 Kernel (`lean/BootstrapTowerTransition.lean`, Lean 4.35.0):** 18 constructive theorems, `0` `sorry`, `0` warnings, `0` `Classical.choice`, `0` `Lean.ofReduceBool`, `3` zero-axiom theorems (`tower_values_exact`, `tower_curve_morris_quarter_exact`, `box_transition_dyadic_exact`).
- **Pure-Stdlib Zero-Float Verifier (`BOOTSTRAP-TOWER-CERT-1.0`, `src/bootstrap_certifier.py`):** Four exact-$\mathbb{Q}$ verification modes covering Theorem 1.1 / Corollary 1.2 tower-curve bijections, exact finite $\mathbb{Z}_n^d$/$[n]^d$ bootstrap filling and Eq. (2.3) $d$-closed barriers, Lemma 2.1 Eq. (2.4) exact rational binomial tails $w_0$, and Eq. (1.2) dyadic box limits $1-2^{-m}$ paired with Lemma 2.3 $\ell_\infty > 2t$ separated rectangles.
- **Canonical Atlas & Adversarial Audit:** 16/16 canonical certificates replayed in `0.036–0.446 ms`, exhaustive $2^9=512$ exact filling polynomials on $\mathbb{Z}_3^2$ vs. $[3]^2$, Eq. (2.3) barrier census on $\mathbb{Z}_4^3$, exact-$\mathbb{Q}$ $w_0$ decay sweeps ($d=10..32$), and **60/60 live adversarial forgery attacks rejected**.

## Quick Reproduction
```bash
# 1. Build & audit Lean 4 kernel
cd lean && lake build && lake env lean AuditAxioms.lean && cd ..

# 2. Regenerate all 16 certificates, 4 exact-Q experiments, and 3 figures
python scripts/run_all_experiments.py
python scripts/make_figures.py

# 3. Run full end-to-end verification (Lean + 16 certs + 60/60 attacks + SHA-256 manifest)
python scripts/verify_all.py
```
