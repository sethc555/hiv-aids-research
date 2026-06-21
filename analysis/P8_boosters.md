# P8 — booster vaccines / immunotherapy: the positive mirror, a titration rule, and where the cures live

_Cope Labs / Seth · 2026-06-20 · model [tip_model_p8_boosters.py](tip_model_p8_boosters.py) ·
figures [p8_titration.png](p8_titration.png), [p8_synergy.png](p8_synergy.png) · run on the
P6-calibrated heterogeneous cohort._

## Why this phase exists

The whole program found what does **not** help: a TIP is neutral to cure because it is
*parasitic on active virus* — it needs the very replication the cure removes. A booster vaccine
or therapeutic antibody is the **opposite**: it raises immune control (`kf`) **directly**, on
the axis the cure actually turns on. So it should help where the TIP could not. P8 quantifies
how — and asks "is there more room?"

## 1. Boosting: a single shot only delays; sustained dosing banks cures

On the calibrated cohort (`kf ~ Normal(10, 2.5)`; "cure" = durable treatment-free remission,
never rebounded through 1200 d):

| strategy | durable cure | rebound delay |
|---|---|---|
| no booster | 3% | 26 d |
| **single** boost | 1% | 81 d |
| repeat q90d × 2 yr, then **stop** | 9% | 75 d |
| repeat **q45d** × 2 yr, then stop | **14%** | 342 d |
| repeat q90d **forever** | 16% | 65 d |

A single boost only **delays** rebound (26→81 d) — the RIO/bNAb pattern, not a cure. Repeated
boosting **banks durable cures** (3%→14%), and q45d×2yr-then-stop (14%) nearly matches dosing
forever (16%): most of the benefit is a *treatment-free* cure won during the held window, not
lifelong suppression. **Frequency beats amplitude** (q45d ≫ q90d) — the first hint of the
titration rule.

## 2. The titration rule — there is a level to maintain, and it's the *trough*

| (1) maintenance threshold (hold kf=L forever) | kf 12–14 | **15** | 16 | 18 |
|---|---|---|---|---|
| durable cure | 3–5% | **14%** | 21% | 43% |

| (2) cure = level × duration (hold kf=15, then stop) | 1 yr | 2 yr | 3 yr | 4 yr |
|---|---|---|---|---|
| durable cure | 1% | 5% | 10% | 14% |

| (3) trough rule — **peak fixed at 18**, vary decay → trough | trough 10 | 12 | **14** |
|---|---|---|---|
| durable cure | 4% | 9% | **49%** |

Three rules fall out:
- **Maintenance threshold ≈ kf 15** (~1.5× the cohort-median immune strength): below it almost
  nothing, above it control climbs steeply.
- **Level sets the ceiling, duration fills it:** holding at 15 caps achievable cure ~14%; you
  reach it over ~3–4 years. Higher level → higher ceiling (18 → 43%).
- **Control tracks the TROUGH, not the peak (the key result):** every row of (3) hits the same
  peak (kf 18), yet cure swings 4%→49% purely on how far immunity decays between doses. Hitting
  a high peak is useless if the valley drops below threshold. **Titration target = keep the
  trough above ~kf 15** — exactly the bNAb-trough logic ("rebound doesn't occur while both
  antibodies stay above ~10 µg/mL"), with dosing interval set by the antibody/immunity half-life.

## 3. "More room": reservoir reduction × boosting is where the big cures live

Boosting alone caps out (~14–43%) because the maintained level only sets a ceiling against a
**fixed-size reservoir**. The untested lever is shrinking the reservoir first (early ART,
shock-and-kill, gene therapy) — the CML "deep response × duration" idea. Holding immunity at
each level against a reduced reservoir:

| reservoir ↓ \ maintained kf | 11 | 13 | 15 | 17 |
|---|---|---|---|---|
| **1× (no reduction)** | 2% | 4% | 14% | 31% |
| 3× smaller | 10% | 18% | 39% | 65% |
| 10× smaller | 31% | 48% | 69% | 84% |
| **33× smaller** | **73%** | 78% | 87% | **96%** |

This is a strong **synergy**, and it's where the room is:
- **Neither lever alone gets you there:** full reservoir + strong immunity (kf 17) = only 31%;
  big reservoir reduction + weak immunity (kf 11, 33×) = 73%; *together* → 96%.
- **Shrinking the reservoir does two things at once:** it **lowers the immune level needed**
  (at 33× smaller, even kf 11 cures 73%) **and raises the ceiling**. The two multiply rather than
  add — exactly the Hill 2014 reservoir-reduction threshold meeting the immune-control axis.

## Net — the project's positive program

Read against the negative result that opened the work: a **TIP is neutral** (it depends on the
virus the cure erases), but **immune boosting is not** — it acts directly on the control axis,
with a clear, transferable dosing logic: *keep the trough above the maintenance threshold (≈1.5×
baseline immunity), at an interval set by immunity's half-life, sustained for ~3–4 years to bank
treatment-free cures — and combine it with reservoir reduction, because the reservoir×immunity
synergy, not either lever alone, is where durable cure rates go from ~15% to ~90%.* That is the
same prescription the cured-disease parallels give (CML TFR: depth **and** duration of response;
bNAb trials: trough levels and combination) — now on one quantitative axis.

_Caveats: `kf` and the reservoir scale are dimensionless model knobs (not a clinical vaccine dose
or a measured log-reduction), so the specific thresholds ("kf 15", "33×") and percentages are
illustrative; the robust, transferable findings are the STRUCTURE — single-shot delays vs
sustained cures, trough-not-peak, frequency-over-amplitude, level-sets-ceiling/duration-fills-it,
and the reservoir×immunity synergy. Same heterogeneous cohort and illustrative reservoir
dynamics as P4–P6; "cure" = active-infection extinction with no rebound, the functional (not
provably sterilizing) endpoint discussed in [P4_findings.md](P4_findings.md). NB (AUDIT2 #6): the
reservoir×immunity "synergy" is, in a stochastic-extinction model, partly **structural** —
fewer rebound seeds × faster killing is super-additive in P(extinction) by construction — so the
robust claim is the **direction** (combination ≫ either lever), not synergy-as-discovery._

**Cross-disease references** (AUDIT2 #4 — these underpin the "cured-disease prescription" framing):
the CML treatment-free-remission "depth **and** duration of response" rule — Saussele et al.,
EURO-SKI, *Lancet Oncol* 2018 (~50% TFR); the CML cure-by-stochastic-extinction result —
Lenaerts, Pacheco, Traulsen, Dingli, *Haematologica* 2010 (P(extinction)→1−1/N); extinction
through an ODE-invisible trough — Baar & Bovier, *Sci Rep* 2016; "lessons from cancer immunology"
for HIV T-cell reservoir control — Mylvaganam, Maus, Walker, *Front Immunol* 2019; bNAb
trough/combination logic — RIO trial, Lee/Fidler et al., *Lancet HIV* 2026.
