# METHODS — within-host TIP-vs-immunity modeling (P1 → P5)

_Cope Labs / Seth · consolidated methods for the TIP↔CD8 modeling program. This documents the
equations, parameters, numerical methods, and audit discipline shared across the phase scripts.
Per-phase scientific results live in the `P*_findings.md` files; this file is the how._

> Scope/status: illustrative within-host models for hypothesis generation and methods, not
> calibrated clinical predictions (P5 adds a ballpark calibration to ATI rebound timing). Read
> [AUDIT.md](AUDIT.md) before quoting any number.

---

## 1. The biological object

A **therapeutic interfering particle (TIP)** is an engineered, conditionally-replicating
defective HIV genome. It carries no functional structural genes, so it produces particles
**only in a cell co-infected by wild-type (WT)** HIV, where it parasitizes WT proteins and
diverts packaging to itself (diversion fraction `ρ`), amplified by a mobilization advantage
`ψ`. The therapeutic hypothesis (Weinberger/Pitchai; Rast/de Boer) is that a TIP lowers WT
viral load. The question this program probes is whether that **helps or fights** a CD8 /
"vaccinal-effect" immune-control cure, since the TIP suppresses WT by removing the very
infected-cell antigen pool that primes CD8.

---

## 2. Core within-host model (P1)

State (concentrations per mL): `T` uninfected target CD4⁺; `Iw` WT-only productive; `It`
TIP-only carrier (no output); `Id` dually-infected (makes WT and TIP); `Vw`, `Vt` free WT /
TIP virions.

```
T'  = λ − d_T·T − b·T·Vw − b_t·T·Vt
Iw' = b·T·Vw            − (d+κ)·Iw − b_t·Iw·Vt
It' = b_t·T·Vt          − d·It      − b·It·Vw
Id' = b_t·Iw·Vt + b·It·Vw − (d+κ)·Id
Vw' = p·Iw + (1−ρ)·p·Id − c·Vw
Vt' = ψ·ρ·p·Id          − c·Vt
```

- `κ` = immune (CTL) killing pressure on **productive** cells — the de Boer knob.
- Production split: a dual cell makes WT at the reduced rate `(1−ρ)p` and TIP at `ψρp`.
- Mass-action infection; standard Perelson-class structure.

**Invariants checked:** WT basic reproductive ratio `R₀ = b·T₀·p/(d·c) = 8.70`; TIP invasion
threshold `ψ ≈ 7.5` (verified by next-generation matrix **and** direct ODE). `R₀_TIP` scales
**sub-linearly (~√ψ)**, not `0.134·ψ` — the TIP transmission cycle has two stages
(Vt→infect→Id→produce Vt), so the spectral radius is a geometric mean (audit correction,
2026-06-20).

---

## 3. Model lineage (what each phase changes, and why)

| phase | change | reason |
|---|---|---|
| **P1** | static immunity `κ` | baseline; test "immunity shrinks TIP efficacy" |
| **P1.1** | dynamic effector ODE `E` | feedback rescue? → instead oscillation + early-TIP immune-blunting |
| **P1.2** | memory + exhaustion compartments | realistic immunity; antagonism persists, statistics fragile |
| **P1.3** | **QSS immunity** `E*(A)=E_max·A/(A+K_E)` | remove effector-lag oscillation; + spatial test |
| **P1.4** | graded visibility `ν`, defective-clone antigen `A_def` | Simonetti-grounded escape test |
| **P1.5** | **QSS free virus** too → 4-cell model | remove free-virus lag; test if escape is stable |
| **P2** | replace static `A_def` with a **dynamic, self-maintaining, TIP-proof reservoir** | decouple the antigen pool — is the antagonism by-construction? |
| **P2b** | **Holling-II** saturating killing `k·E/(1+E/E_sat)` | last mean-field variation; bounded predation |
| **P3** | **stochastic tau-leaping** of the cell model | the escape regime oscillated → discreteness resolves it |
| **P3b** | stochastic phase map over (reservoir, ν) | P(cure)/P(TIP-loss) surface |
| **P3c** | dt→0 convergence + exact reduced-scale SSA | validate the stochastic result |
| **P4** | **latent replication-competent reservoir** + ART→ATI | close the model↔clinic gap: "control" = reservoir-rebound suppression |
| **P4b** | sustained dosing, high ψ, 3× reactivation | robustness of the TIP-neutral verdict |
| **P4c** | record ATI rebound time | check rebound clock is clinically realistic |
| **P5** | fit reactivation rate + CD8 strength to clinical ATI curves | quantitative (not just qualitative) verdict |

### Quasi-steady-state (QSS) reductions
- **Immunity (P1.3+):** effectors equilibrate fast relative to cell turnover →
  `E*(A) = E_max·A/(A+K_E)`, with per-cell killing `κ = k·E*`. Antigen pool `A` (see §4).
- **Free virus (P1.5+):** virions clear at `c=23/day ≫ d=1/day` → `Vw* = (p/c)(Iw+(1−ρ)Id)`,
  `Vt* = (p/c)(ψρ Id)`. Collapses to a 4-variable cell-only competition (`T,Iw,It,Id`).

### Two-pool antigen (P2) — the central methodological move
The CD8-priming antigen is split into the **active** pool (suppressible by the TIP) and a
**reservoir** pool (not): `A_eff = Iw + ν·Id + L`, where `L` is a dynamic, self-maintaining,
TIP-proof clone (Simonetti defective clones), `dL = g·L(1−L/R) + s·(Iw+Id)`. This tests whether
the TIP↔CD8 antagonism is an artifact of a single shared antigen pool (it is not — see
[P2_findings.md](P2_findings.md)).

### Stochastic engine (P3+)
Demographic **tau-leaping**: integer populations; each reaction channel fires
`Poisson(propensity·dt)` per step; superinfection consumption capped at available cells;
populations clipped at 0. Free virus kept at QSS (computed from cell counts each step). The
deterministic limit-cycle troughs become real **extinction** events — the qualitative reason
the stochastic model answers what the ODEs could not.

### Reservoir + ATI (P4) — 7 compartments
Adds `L_def` (defective antigen clone, primes CD8, TIP-proof), `L_lat` (latent
replication-competent provirus — **transcriptionally silent, absent from `A_eff`** — reactivates
to `Iw`), and `L_latd` (engineered latent dual → reactivates to `Id`). Schedule:
**CHRONIC** (`art=1`) → **ART** (`b·art`, `art=10⁻³`) → **ATI** (`art=1`). A fraction `F_lat` of
new infections enter latency. "Post-treatment control" = productive cells `(Iw+Id) < 50` at the
end of a 500-day ATI; the latent reservoir persisting in controlled replicates makes this a
**functional**, not sterilizing, cure.

---

## 4. Parameters

Standard Perelson-class within-host HIV values; **illustrative, not patient-fitted** (except the
P5 calibration targets). Defined in `tip_model.py` (`P`, `QS`) and `tip_model_p4_reservoir.py`.

| symbol | value | meaning | basis |
|---|---|---|---|
| λ | 1e4 /mL/day | target cell production | Perelson-class |
| d_T | 0.01 /day | target death (→ T₀=1e6) | standard |
| d | 1.0 /day | productive infected-cell death | standard |
| c | 23 /day | virion clearance | standard |
| p | 2000 /cell/day | virion production | standard |
| b (=b_t) | 1e-7 mL/day | infection rate (→ R₀=8.70) | tuned to R₀ |
| ρ | 0.90 | TIP packaging-diversion fraction | TIP design |
| ψ | 22 (P1.5+) | TIP mobilization advantage | above invasion threshold ~7.5 |
| E_max | 3e4 | max effector level | QSS immunity |
| K_E | 3000 | antigen at half-max effector | QSS immunity |
| k | 3e-5 /day | per-effector killing constant | QSS immunity |
| ν | 0–1 | dual-cell immune visibility (1 = full) | escape lever |
| g, R | 0.02 /day, 5000 | defective-clone timescale, size | P2/P3 reservoir |
| F_lat | 5e-4 | fraction of infections entering latency | P4 (illustrative) |
| A_react | 1e-3 /cell/day | latent reactivation rate | P4 (→ realistic rebound clock, P4c) |
| D_L, P_L | 1e-3, 1e-3 /day | latent death, homeostatic proliferation | P4 (long-lived reservoir) |
| kf | sweep | CD8 "vaccinal-effect" strength (× baseline killing) | P4 cure knob |

---

## 5. Numerical methods & reproducibility

- **ODE integration:** `scipy.integrate.solve_ivp`, `LSODA`. P1 phase sweeps use `rtol=1e-6,
  atol=1e-2, max_step=2.0` (a known checkerboard at thresholds — qualitative structure robust,
  cell-level values noisy). The QSS builds (P1.3/P1.5/P2) tighten to `rtol=1e-8, atol=1e-3,
  max_step=5.0`.
- **Tau-leaping:** fixed `dt` (0.03–0.04 production; convergence checked at 0.02/0.01/0.005).
  Replicates vectorized as array rows; per-replicate parameters (ν, kf, ψ, R) supported so a
  whole grid runs in one step-loop.
- **Exact SSA (P3c):** event-by-event Gillespie with `searchsorted` channel selection, run at a
  reduced system size (volume factor `ω`: counts = concentration·ω, rates rescaled to preserve
  the deterministic dynamics) so it is tractable; validates *structure*, not the full-size %.
- **Determinism:** `numpy` RNG seeded per run; figures/`.npz`/`raw_cache.json` committed.
  Scripts write outputs next to themselves (`HERE = os.path.dirname(...)`), so they are
  path-portable. Dependencies: `numpy`, `scipy`, `matplotlib` (see `requirements.txt`).

---

## 6. Audit discipline (the non-negotiables)

1. **Convergence gate:** never call a steady state / "escape" real unless converged — report
   tail max/min over the last window; `>1.05` = oscillating, excluded.
2. **Stable-only / absorbing-only reporting:** quantitative claims read only off converged ODE
   cells or absorbing stochastic outcomes; oscillating-regime magnitudes are not quoted.
3. **Distributions, not single runs:** stochastic claims are probabilities over ≥80 replicates;
   robustness via dt-halving and an independent exact algorithm.
4. **Engagement check:** "neutral" is only reported after confirming the agent (e.g. the TIP)
   was actually active (Vt>0 / engaged in a stated fraction of replicates), not merely absent.
5. **Adversarial review:** every claim went through self-audit + multi-agent audit; retractions
   and corrections are logged in [AUDIT.md](AUDIT.md), not silently edited.

---

## 7. Calibration methodology (P5)

**Goal:** move the rebound clock and CD8 cure-threshold from "illustrative" to "anchored," then
re-derive the TIP verdict at the anchored operating point.

- **Targets (clinical, sourced in [P5_findings.md](P5_findings.md)):** (i) untreated/placebo ATI
  time-to-rebound (median and, where available, Kaplan-Meier fractions); (ii) bNAb-arm
  delayed-rebound timing (RIO-class); (iii) the spontaneous post-treatment-control fraction.
- **What is fit:** the latent **reactivation rate `A_react`** sets the rebound clock → fit to the
  untreated/placebo median. The **CD8 strength `kf`** sets the rebound delay / control fraction
  → fit to the bNAb arm and the PTC fraction. Everything else stays at the §4 values.
- **Objective:** minimize squared error between model median time-to-rebound (and KM fractions)
  and the clinical targets, over a grid/1-D search in (`A_react`, `kf`); rebound time = first
  crossing of a detection-analog `Vw` threshold, measured exactly as in P4c.
- **Identifiability caveats (documented, not hidden):** rebound timing constrains the *product*
  of reservoir size and reactivation rate more than each alone; the detection threshold choice
  shifts absolute times; fitting to summary statistics ≠ fitting patient-level curves;
  parameters are correlated. P5 reports the fit **and** its uncertainty, and treats the
  recovered TIP verdict as conditional on these caveats.

---

_Last updated 2026-06-20. Equations cross-checked against the phase scripts; any divergence
between this file and a script is a bug — the script is ground truth._
