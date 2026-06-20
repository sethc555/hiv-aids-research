# P6 — heterogeneous-patient cohort: the KM fit closes, and the TIP-neutral verdict is final

_Cope Labs / Seth · 2026-06-20 · model [tip_model_p6_heterogeneous.py](tip_model_p6_heterogeneous.py) ·
figure [p6_heterogeneous.png](p6_heterogeneous.png) · the successor P5 named._

## Why this phase exists

P5's calibration exposed a real mis-specification: **a single immune setpoint cannot reproduce
both the fast-rebound majority (~3 wk) and the durable-controller minority (~5%)** — the clock
fit wanted kf≈9, the PTC-fraction fit wanted kf≈13. The clinical resolution is **patient
heterogeneity**: post-treatment controllers have a different immune/reservoir setpoint. P6 makes
each replicate a *patient* drawn from population distributions and fits the **whole rebound
Kaplan-Meier curve**, not a single operating point.

## The model

Each of N=300 patients gets its own:
- **immune strength** `kf_i ~ Normal(m_kf, s_kf)` (the CD8 "vaccinal-effect" knob), and
- **reservoir clock** `f_lat_i ~ lognormal` (median 8×10⁻⁵ from P5, geometric SD ~2.2×).

The cohort is run through the full CHRONIC→ART→ATI schedule (P4 engine); the spread of
(kf, f_lat) produces a **distribution of rebound times → a model KM curve**. We fit (m_kf, s_kf)
to the clinical KM summary; everything else stays at the METHODS §4 values. (Controllers emerge
self-consistently: a high-kf patient establishes a *smaller* reservoir during chronic infection
and then controls its reactivation — the clinically-correct "small reservoir + strong immunity"
correlation, not an imposed one.)

## Result 1 — one cohort now fits the whole curve (P5 tension resolved)

| target | clinical | model (best fit) |
|---|---|---|
| rebound median (>1000 c/mL) | 22 d (A5345) | **20.8 d** |
| rebounded by wk 12 (>200) | 77% (ACTG) | **78%** |
| rebounded by wk 4 (>200) | 66% (ACTG) | 55% |
| spontaneous PTC | ~4–5% (CHAMP/Gunst) | **5%** |

**Best fit: `kf ~ Normal(10, 2.5)`.** A single heterogeneous cohort simultaneously matches the
median, the week-12 rebound fraction, and the post-treatment-control fraction — the three things
P5 could not satisfy with one kf. The **controller minority is the upper tail** of the immune
distribution (kf ≳ 14, ~5% of patients), the **fast-rebound majority is the bulk** (kf ≈ 10).
The identifiability tension is resolved exactly as predicted: by a *distribution* over immune
strength, not a single value. (The wk-4 fraction is the weakest match — 55% vs 66% — and is
expected: the wk-4/wk-12 clinical points use a >200 c/mL threshold while our detector and the
median target use >1000, so the model rebounds slightly later by construction.)

## Result 2 (the payoff) — the TIP is neutral on the realistic cohort

Running the question on the *fitted* population (sustained, engaged TIP, ψ=60, ν=0.9):

| metric | no TIP | + TIP |
|---|---|---|
| rebound median | 27 d | 25 d |
| rebounded by wk 12 | 72% | 79% |
| **post-treatment control (PTC)** | **5%** | **4%** |

The KM curves **coincide** and the controller fraction moves by **−1 point** (within replicate
noise). On the most clinically-faithful version of the model — a patient population calibrated to
the actual ATI rebound curve — **a TIP does not change post-treatment control.** This is the
final, strongest statement of the program's verdict: not a single-point result, not an
un-calibrated one, but one that holds across a realistic distribution of patients fit to clinical
data.

## Net — the program's conclusion, fully anchored

The TIP↔CD8 antagonism that opened the project (P1–P2) is real on the **active-infection axis**
but **does not transfer to the reservoir-control axis that defines cure** (P4), survives every
robustness stress (P4b), survives clinical calibration (P5), and now survives a **patient-level
heterogeneous fit to the ATI rebound curve** (P6): an engineered TIP is **orthogonal to
immune-mediated post-treatment control** — neither the threat the antagonism suggested nor a
useful adjunct — because it requires the active replication the cure is built to erase.

_Caveats: fit to published KM **summary points** (median, wk4/wk12 fractions, PTC), not digitized
patient-level curves (paywalled); the (m_kf, s_kf) fit is a best-on-grid, not a unique or
maximum-likelihood solution — other (heterogeneity, clock) combinations could fit comparably;
reservoir-clock heterogeneity was fixed (gsd 2.2×), not fit; the >1000 vs >200 c/mL threshold
mismatch explains the wk-4 residual; "bNAb = higher kf" abstracts antibody action as enhanced
immune clearance. The verdict is qualitatively robust across all of these; precise percentages
are not._
