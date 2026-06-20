# P1 — TIP vs immunity: first-pass within-host result

_Cope Labs / Seth · 2026-06-19 · model: [tip_model.py](tip_model.py) ·
figures: [phase_diagram.png](phase_diagram.png), [timecourse.png](timecourse.png)_

## Question

Does an immune (CTL) response against infected cells shrink the regime where a
**therapeutic interfering particle (TIP)** suppresses wild-type (WT) HIV? This is
the live disagreement in the TIP modeling community:

- **Dodd & de Boer 2026 (J Theor Biol, vol 619; online 2025), pessimistic:** *"Even a moderate immune
  response against virally infected cells drastically decreases the range of
  parameter values for which therapy is effective."*
- **Weinberger / Pitchai 2024 (Science), optimistic:** TIPs reduced viral load and
  disease durably in macaques (which have intact immune responses).

## Model

Standard Perelson-class within-host HIV ODE, extended so the TIP is produced
**only in dually-infected cells `I_d`** (where it hijacks WT proteins and diverts
packaging, fraction `rho=0.9`, amplified by mobilization advantage `psi`).
Immunity enters as a killing pressure `kap` (/day) added to the death rate of
**productive** infected cells `I_w, I_d` — the de Boer mechanism in reduced
(quasi-steady-CTL) form. WT R0 = 8.7. Full equations in the script header.

**Experiment:** marginal TIP benefit = log10(WT setpoint without TIP) −
log10(WT setpoint with TIP), swept over immune pressure `kap ∈ [0, 1.2]` ×
TIP advantage `psi ∈ [1, 20]`. The "without TIP" baseline already includes
immunity at each `kap`, so the metric isolates the TIP's **marginal** value on
top of whatever immunity alone achieves.

## Result — immunity is first-order, not a detail

| immune pressure `kap` (/day) | infected-cell lifespan | `psi` needed for ≥1-log TIP benefit |
|---|---|---|
| 0.0 (none) | 1.0 d | **~12** |
| 0.4 (moderate) | 0.7 d | **~17.6** |
| 0.6 | 0.63 d | **none ≤20** |
| 0.8 | 0.55 d | **none ≤20** (max benefit 0.02 log) |

- Without immunity, a strong TIP (`psi=12`) drives a **1.44-log** WT reduction and
  persists — TIPs work, reproducing the Weinberger result.
- The fraction of the `psi`-range giving ≥1-log benefit **collapses from 28%
  (kap=0) to 0% by kap=0.6**.
- The persistence edge (white contour) and the `psi`-required line (red) both rise
  steeply with `kap`: every increment of immune killing demands a markedly
  better-engineered TIP, and past a *moderate* level **no achievable `psi`
  rescues it**.

**Verdict: this independent model lands on de Boer's side.** A moderate CTL
response — exactly the kind present in the Pitchai macaques and in any treated
human — sharply shrinks (and then closes) the TIP's efficacy window. Immunity
belongs in the *center* of TIP design, not the appendix.

**Constructive output:** the red curve is a **design spec** — `psi_required(kap)` —
the mobilization advantage a TIP must hit to overcome a given immune pressure.
It says where the engineering target is, and where it becomes unreachable.

## The crux this exposes (→ the real open question)

The reduced model treats immunity as a **fixed** killing pressure `kap`. But a
working TIP *lowers viral load*, which lowers antigen, which should *lower* immune
activation — a negative feedback that could **partially rescue** the TIP (fewer
dual cells killed once WT is suppressed). My static-`kap` model omits exactly this
loop — and it is precisely where the Weinberger camp's optimism would live.

So the genuine, unresolved question is **not** "does immunity hurt TIPs" (yes, both
de Boer and this model agree) but: **does the TIP→antigen→immunity feedback create
a rescued regime that a static-immunity model misses?** Neither a constant-`kap`
model (here) nor a worst-case immune model settles it. That is the contribution
opening.

## Honest caveats

- **Reduced immunity.** Constant `kap`, not a dynamic effector `E(t)` with
  antigen-driven proliferation. Captures the mechanism, omits the feedback (above).
- **Thresholds are parameterization-dependent.** The absolute `psi≈12` invasion
  threshold shifts with `rho`, `b`, burst size; the *qualitative* collapse with
  immunity is robust, the exact numbers are not calibrated to data.
- **Conservative toward the TIP.** TIP entry `bt=b` (no superinfection advantage);
  carriers `I_t` spared by immunity. Both make the TIP look *better* than worst
  case, so the pessimistic result is not an artifact of stacking the deck.
- Numerical checkerboard near thresholds (LSODA, marginal-metric sensitivity);
  structure is clear but the boundary is ±1 grid cell noisy.

## Next — P1.1 (the actual novel step)

Add a **dynamic immune compartment** `E(t)` with antigen-driven proliferation
(`dE/dt = a·E·(I_w+I_d)/(I_w+I_d+K) − d_E·E`) and killing `k·E·(I_w+I_d)`. Then
test the open question: **sweep immune responsiveness `a` × TIP advantage `psi`
and ask whether the TIP→antigen→immunity feedback opens a rescued efficacy regime
that the static-`kap` map (this result) declares dead.**

If a rescued regime exists, that reconciles de Boer (static) with Pitchai
(macaques, live immunity) and is a publishable, no-lab result. If it doesn't, that
strengthens de Boer and reframes TIP design around immune evasion. Either outcome
is a real answer to a question two labs currently disagree on.
