# P3b/P3c — the P(cure) phase map, and validation (dt→0 + exact SSA)

_Cope Labs / Seth · 2026-06-20 · models [tip_model_p3b_phasemap.py](tip_model_p3b_phasemap.py),
[tip_model_p3c_ssa.py](tip_model_p3c_ssa.py) · figure [p3b_phasemap.png](p3b_phasemap.png)._

## P3b — the escape as a probability surface

P3 profiled one stochastic point and found the immune-compatible TIP is a low-probability
outcome. P3b turns that into the surface the design objective actually lives on: **P(cure)**
and **P(TIP-loss)** over (reservoir size R, dual-cell visibility ν), 80 replicates/cell,
single vectorized tau-leaping pass (ψ=22, slow reservoir).

**P(cure) = P(WT-producing cells go extinct), rows ν, cols R/A0:**

| ν \ R | 0× | 0.5× | 1.0× | 1.4× | 1.9× |
|---|---|---|---|---|---|
| 0.10 | 0% | 5% | 22% | 28% | 29% |
| 0.25 | 0% | 1% | 9% | **31%** | 29% |
| 0.40 | 0% | 0% | 2% | 0% | 5% |
| 0.55 | 0% | 0% | 1% | 0% | 0% |
| ≥0.70 | 0% | 0% | 0% | 0% | 0% |

- **A cure requires BOTH a strong immune-evasion regime (ν ≲ 0.25) AND a substantial reservoir
  (R ≳ 1×A0).** Either alone gives ≈0%. The reservoir is the antigen floor that keeps CD8
  armed while the TIP suppresses active WT; evasion is what lets the TIP's carrier cells
  survive that CD8 — both are necessary, and the cure corner is their intersection.
- **Max P(cure) = 31%** (ν=0.25, R≈1.4×A0). Even the best corner loses ~2/3 of the time.
- **Zero stable-coexistence cells anywhere on the grid** — confirms P3's no-coexistence
  finding across the whole (ν, R) plane, not just the one point. The outcome is always the
  lottery; only the odds move.

This is the honest analogue of the deterministic phase diagrams (P1/P1.5/P2 plotted a
log-reduction or a yes/no escape): the stochastic phase diagram is a **probability surface**,
and the design problem is to climb it — maximize P(cure)/P(TIP-loss) — not to find a fixed
point that does not exist.

## P3c — validation: the result is not a tau-leaping artifact

Two independent checks of the P3 headline at (ν=0.2, R=5000):

**(1) dt→0 convergence (full scale, 160 reps).** Tau-leaping → exact SSA as dt→0, so this
directly bounds the discretization error:

| dt | cure | TIP-lost | coexist | CD8 in cures |
|---|---|---|---|---|
| 0.020 | 15% | 85% | **0%** | **82%** |
| 0.010 | 14% | 86% | **0%** | **82%** |
| 0.005 | 10% | 90% | **0%** | **82%** |

**(2) Exact Gillespie SSA, reduced scale** (ω=0.03 system-size factor, preserving the
concentration dynamics; 30 reps; counts = concentration·ω). An event-by-event algorithm with
no time discretization at all: **cure 20%, TIP-lost 80%, coexist 0%, CD8 82%** — and the
lottery resolves fast (absorbing within ~4–35 days at reduced scale).

**What is validated:** the **structural** claims are method-invariant across tau-leaping
(dt 0.005–0.03) and exact SSA — **no stable coexistence (0% everywhere)**, and **CD8 retained
at ~82% in the controlled outcomes** (exact to the point). The cure **magnitude** is
**order-15% (≈10–20% across methods)**, with a mild dt-trend (15→10% as dt→0) that puts the
honest point estimate slightly below P3's first-pass 16% but well within the same regime: a
**low-probability** outcome, not a design point. So P3's headline stands, with the precise
percentage softened to a range and the qualitative core independently confirmed.

## Net

The escape is a **probability surface** with a cure corner at low ν / high R, peaking ~31%
and zero elsewhere; coexistence never occurs; and the reservoir-keeps-CD8-primed mechanism
(~82%) is exact-algorithm-validated. The next grind is the obvious one: a **latent
replication-competent reservoir** so "WT controlled" means reservoir eradication (closing the
model↔clinic loose analogy that has shadowed every phase), and then optimal-control of the
TIP dose/timing to climb the P(cure) surface.

_Caveats: ω-reduced SSA validates structure, not the full-size percentage (smaller systems
absorb faster and noisier — hence 20% vs ~13%); illustrative params; "cure" = active-infection
extinction, still NOT reservoir eradication; the surface is one ψ slice._
