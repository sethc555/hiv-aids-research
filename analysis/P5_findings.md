# P5 — clinical calibration: the TIP verdict survives anchoring, and the fit exposes a model limit

_Cope Labs / Seth · 2026-06-20 · model [tip_model_p5_calibration.py](tip_model_p5_calibration.py) ·
figure [p5_calibration.png](p5_calibration.png) · methodology in [METHODS.md](METHODS.md) §7._

## Goal

Move the reservoir rebound clock and the CD8 cure-threshold from "illustrative" to "anchored
to clinical ATI data," then re-derive the TIP verdict at the calibrated operating point — and
report the fit's uncertainty honestly.

## Clinical calibration targets (sourced)

| quantity | value | definition | n | source |
|---|---|---|---|---|
| placebo/untreated rebound median | **16 d** (>50), **21 d** (>400), **32 d** (>10⁴ c/mL) | confirmed VL crossing | 382 | Gunst, *Nat Commun* 2025;16:906 |
| placebo rebound median (modern ART) | **22 d** (range 13–230) | first VL ≥1000 c/mL | 45 | Li, *Clin Infect Dis* 2022 (A5345) |
| rebounded by wk4 / wk12 | **~66% / ~77%** | ≥200 c/mL confirmed | 235 | ACTG pooled (PMC4911279) |
| spontaneous post-treatment control | **4%** chronic / **13%** early | VL≤400 ≥24 wk off-ART | 460/148 | Namazi, *JID* 2018 (CHAMP) |
| **bNAb durable control (RIO)** | **24%** (7/29) vs **6%** placebo | VL<1000 beyond 96 wk | 34/arm | Lee/Fidler, *Lancet HIV* 2026 |
| reservoir reactivation (modeling) | **~0.125/d** (Davenport) … **~4/d** (Hill) | successful reactivation | — | Pinkevych 2015 / Hill 2014 PNAS |
| reservoir t½ | **44 mo** (δ=5.2e-4/d) | LR decay on ART | — | Siliciano, *Nat Med* 2003 |

Our detection analog is `Vw>1000 c/mL`, so the **22 d** A5345 median is the matched rebound-clock
target; the **4–6% / 24%** PTC fractions anchor the CD8 strength `kf`.

## Method (see METHODS §7)

A derived fact: at ATI the reservoir reactivation **flux** = `f_lat × (chronic infection flux)`
(reservoir size ~`f_lat·flux/react`, and reactivation `= react·size` cancels `react`). So the
**clock is set by `f_lat`**, fit to the placebo rebound median; **`kf`** (the CD8 "vaccinal
effect") is anchored to the durable-control fractions. Control is measured by **one consistent
metric — "never rebounded"** (`Vw` stays <detection through ATI), matching how trials measure
time-to-rebound. (An earlier run mixed two control definitions and produced a spuriously large
TIP effect; unifying the metric removed it — a methods lesson logged here.)

## Results

**1. The rebound clock calibrates well.** `f_lat = 8×10⁻⁵` gives a placebo rebound median of
**~19 d** (at placebo-level CD8), squarely in the clinical band (16–22 d at the >1000 threshold).

**2. The control curve brackets the clinical fractions.** Never-rebound control vs CD8 strength:
`kf`=11→0%, 12→1%, 13→3%, 14→13%, 15→29%, 16→41%, 18→75%. So **placebo PTC ~5% ≈ kf 13** and
**bNAb ~24% ≈ kf 15** — both clinical anchors fall on the curve.

**3. The TIP verdict survives calibration.** With the consistent metric, a sustained, engaged
TIP (ψ=60, ν=0.9) shifts control by only **+2–5 percentage points** everywhere — at the
calibrated placebo point 3%→6%, at the bNAb point 29%→31%. **P4b's neutrality holds at the
clinically-anchored operating points.** (Any larger apparent "help" at very low `kf` is
TIP-*dependent* suppression requiring perpetual dosing — the P4b "second ART" regime, not a
durable immune cure.)

**4. The fit exposes a genuine model limitation (the most useful result).** The two placebo
anchors are **mutually inconsistent under a single `kf`**: the rebound-median target (22 d)
wants `kf≈9`, but the PTC-fraction target (5%) wants `kf≈13`, where the median is ~80 d. A
homogeneous-parameter model **cannot** produce the real clinical pattern — a **fast-rebound
majority (~95%, median ~3 wk) plus a durable-controller minority (~5%)** — at one immune
setpoint. Clinically that bimodality reflects **patient heterogeneity** (controllers have a
distinct immune/reservoir setpoint). So the calibration's clean message is: *to fit ATI cohorts
quantitatively the model needs a distribution over `kf` (and/or reservoir size), not a single
value.* This is a mis-specification the calibration **surfaces**, not a bug.

## Net

Calibration does three things. (i) It **anchors the rebound clock** to clinical reality
(`f_lat=8e-5` → ~19 d median, matching A5345's 22 d). (ii) It **confirms the TIP verdict is not
an artifact of un-calibrated parameters** — at the operating points that reproduce the clinical
PTC fractions, a TIP moves durable control by only a few points (P4b neutrality holds). (iii) It
**identifies the next modeling requirement**: population heterogeneity in immune strength, since
no single `kf` matches both the fast-rebound majority and the controller minority. The
qualitative conclusion of the whole program — *a TIP is essentially neutral to immune-mediated
post-treatment control* — is now anchored, not just illustrated.

_Caveats: fit to summary statistics (medians, KM endpoints, PTC fractions), not patient-level
rebound curves; coarse grids (f_lat 5 pts, kf 7 pts) with ~±a-few-% replicate noise at 150 reps;
"bNAb = higher kf" abstracts antibody action as enhanced immune clearance; the
clock/PTC-fraction inconsistency means the single-`kf` calibration is a compromise, not a true
joint fit — the heterogeneous-`kf` successor is the right tool. Reactivation-rate literature
spans ~30× (Davenport vs Hill); our `f_lat`-based clock fit sidesteps that dispute but the
implied per-cell rate inherits its uncertainty._
