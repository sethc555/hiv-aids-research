# P15 — global sensitivity: the coupled-TIP benefit is robust, not a tuned corner

_Cope Labs / Seth · 2026-06-21 · model [tip_model_p15_sensitivity.py](tip_model_p15_sensitivity.py) ·
pre-work for contacting TIP modelers — answers "did you cherry-pick parameters?"_

## Why this phase exists

P11–P14 establish the coupled-TIP benefit at chosen operating points. The obvious reviewer
objection: *is the result a tuned corner?* P15 varies **all** key parameters **jointly** at random
over plausible ranges and asks whether the benefit (and the no-backfire) survive, and which
parameters drive it.

## Method

30 Monte-Carlo draws, all six knobs varied together: immune strength `kf∈[9,14]`, dual-cell
visibility `ν∈[0.1,1]`, TIP mobilization `ψ∈[10,100]`, latency fraction `f_lat∈[2e-5,2e-4]`
(reservoir clock), antigen floor `RDEF∈[2000,10000]`, reactivation rate `A_react∈[5e-4,2e-3]`.
For each draw: `effect = P(control | χ=1) − P(control | χ=0)` (coupled TIP vs the de Boer limit).

## Result — robust across the joint space; never backfires

| outcome over 30 joint draws | fraction |
|---|---|
| **HELPS** (≥ +5 pts) | **63%** (mean effect **+17.5 pts**) |
| ~neutral (−5…+5) | 37% |
| **HARMS** (≤ −5 pts) | **0%** (worst case −1 pt = noise) |
| **no backfire** | **100% of draws** |

The coupled-TIP benefit is **not a tuned corner**: across randomly combined parameters it helps or
is neutral every time, and **never harms**. This is the strongest robustness statement in the
project — it varies the load-bearing knobs *together*, not one-at-a-time (the exact gap
[AUDIT2.md](AUDIT2.md) flagged as the prior audits' blind spot).

## What drives the benefit (parameter importance)

Spearman rank correlation of each parameter with the TIP effect:

| parameter | ρ | reading |
|---|---|---|
| **immune strength `kf`** | **−0.56** | dominant: the TIP helps **most where immunity is marginal** (room to push); little where immunity already wins (ceiling) |
| antigen floor `RDEF` | −0.32 | smaller floor / lower baseline control → more room for the TIP to help |
| visibility `ν` | −0.31 | more-evasive TIP cells help slightly more (they persist and keep suppressing WT) |
| mobilization `ψ` | +0.13 | stronger TIP → marginally more help |
| `f_lat`, `A_react` | ≈0 | reservoir clock / reactivation rate barely matter to the *effect* |

So the benefit is governed by **how much headroom immunity leaves** (marginal controllers gain
most), consistent with P14's intermediate-immunity peak — and it is **insensitive to the
reservoir's exact timing**, which is reassuring given those parameters are the least certain.

## Net

Combined with P14 (the phase boundary + falsifiable prediction), P15 gives the package a within-host
modeler expects: a continuous phase diagram that reduces to the prior model, a stated testable
prediction, and a **global** robustness check showing the coupled-TIP benefit holds across joint
parameter variation and never backfires. It remains a **modeling hypothesis** under the project's
standing caveats — illustrative params, one model family, functional-not-sterilizing control — but
the "cherry-picked" and "fragile" objections are now answered with evidence.

_Caveats: 30 Monte-Carlo draws (illustrative — not a full Sobol decomposition); coupling fixed at
χ=1 (a well-coupled TIP — the achievable-coupling question is separate, see P14's χ*); ranges are
plausible-but-illustrative; rank correlations indicate direction/importance, not a variance
decomposition._
