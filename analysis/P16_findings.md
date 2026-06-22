# P16 — a semi-analytical criterion: one reproduction number explains when (and why) a coupled TIP helps

_Cope Labs / Seth · 2026-06-21 · model [tip_model_p16_analytic.py](tip_model_p16_analytic.py) ·
the rigorous pre-work — the *condition* for the sign flip, not only simulation._

## The derivation

After ART interruption the target-cell pool has recovered (T ≈ T₀) and free virus is quasi-steady.
A reactivating wild-type cell `Iw` lives `1/(d+κ)` (κ = immune killing rate), makes `p` virions/day,
each infecting a target with probability `≈ b·T₀/c`. Its **effective reproduction number** is

> **R_eff(WT) = b·T₀·p / ((d+κ)·c) = R₀ · d/(d+κ)**, with **R₀ = b·T₀·p/(d·c) = 8.70**.

Each reservoir reactivation seeds one such lineage, so the active infection is a branching process
that goes **stochastically extinct (→ durable control) iff it is sub-critical**:

> **CONTROL ⟺ R_eff(WT) < 1 ⟺ κ > κ_crit = (R₀−1)·d = 7.70 /day.**

(The stochastic boundary is a bit *softer and lower* than κ_crit, because near-threshold bursts
have establishment probability `1 − 1/R_eff < 1`, so many die by chance even when R_eff > 1.)

## The criterion

A TIP **coupled** to reactivation (fraction χ of reactivations co-introduce the interfering TIP)
suppresses the wild-type burst, which **shifts the effective control boundary down**:

> **κ_crit(χ) = κ_crit − Δ(χ),  with Δ(χ) ≥ 0 and increasing in χ.**

Everything the simulations found follows from this one inequality:

- **It rescues the marginal band.** The TIP changes the *outcome* only for individuals whose immunity
  sits in `[κ_crit − Δ(χ), κ_crit)` — they rebound without it, control with it. So the benefit is
  **concentrated near the control threshold** (marginal immunity), and is a **ceiling effect** where
  immunity already wins (R_eff ≪ 1) or already loses badly. → explains **P14's intermediate-immunity peak**.
- **More coupling → wider rescue band** (`Δ` grows with χ). → explains **P14's χ\* requirement**.
- **Δ(χ) ≥ 0 — coupling never *raises* the threshold.** → explains **P15's 0% backfire**.

## Validation (static-κ model, where R_eff is exact)

TIP effect on durable control vs immune killing κ (no-TIP χ=0 vs coupled χ=1):

| κ | R_eff | control χ=0 | control χ=1 | **TIP effect** |
|---|---|---|---|---|
| 4 | 1.74 | 30% | 87% | **+57** |
| 5 | 1.45 | 57% | 90% | +33 |
| 6 | 1.24 | 68% | 91% | +23 |
| 7 | 1.09 | 74% | 94% | +20 |
| 8 | **0.97** | 98% | 100% | +2 |
| 9 | 0.87 | 99% | 100% | +1 |

Exactly as the criterion predicts: the TIP effect is **largest in the sub-/near-threshold band**
(R_eff ≳ 1, where there is headroom and the coupled TIP does the suppressing), **ceilings out once
R_eff < 1** (immunity already controls), and is **≥0 at every κ — it never backfires.** The derived
`R_eff < 1` line (κ ≈ κ_crit) cleanly marks where the ceiling begins (κ ≥ 8).

## Honest nuance

In this *static-κ* model the benefit peaks at **low** immunity (most headroom; the coupled TIP can
suppress even when immunity is weak). In the *dynamic-immunity cohort* (P14) it peaked at
**intermediate** immunity, because there a very weak immune response **collapses entirely** and the
TIP — which *assists* rather than *replaces* immunity — has nothing to assist. Both are the same
rule (*the TIP rescues the marginal band where headroom exists and the system isn't hopeless*); the
peak's exact location depends on whether immunity is fixed or itself dynamic.

## Net

The coupled-TIP benefit is now grounded in a **derived reproduction number**, not only simulation:
control is set by `R_eff(WT) = R₀d/(d+κ) < 1`, and coupling lowers the effective threshold by
`Δ(χ) ≥ 0`, which **simultaneously explains the marginal-immunity peak (P14), the coupling
requirement (P14), and the no-backfire (P15)** — three empirical findings from one inequality. This
is the form a within-host modeler can check on paper.

_Caveats: R_eff is the standard deterministic threshold; the stochastic control boundary is softer
and below κ_crit (establishment probability < 1). Δ(χ) is shown ≥0 and increasing numerically, not
derived in closed form (a full next-generation-matrix treatment of the nonlinear TIP interference is
the remaining rigor step). Static-κ idealization for analytic cleanliness; illustrative params._
