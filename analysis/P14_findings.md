# P14 — the coupling phase diagram: where a TIP flips from neutral (de Boer) to helpful, and a testable prediction

_Cope Labs / Seth · 2026-06-21 · model [tip_model_p14_coupling_phase.py](tip_model_p14_coupling_phase.py) ·
figure [p14_coupling_phase.png](p14_coupling_phase.png) · pre-work for contacting TIP modelers._

## Why this phase exists

P11–P13 established that a *coupled* TIP helps and never backfires. P14 makes that a **phase
boundary** a modeler can interrogate: the TIP's effect on durable post-treatment control as a
**continuous** function of coupling `χ` (fraction of reservoir reactivations carrying the TIP)
and immune strength `kf`. This does three things reviewers ask for: shows the model **reduces to
de Boer** in the right limit, locates **where** the effect lives, and yields a **falsifiable
prediction**.

## Result 1 — the model reduces to de Boer at χ=0 (reduction-to-prior-model check)

P(durable control), rows = coupling χ, cols = immune strength kf:

| χ \ kf | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|
| **0.0 (de Boer limit)** | 0% | 4% | 24% | 57% | 87% |
| 0.4 | 0% | 4% | 30% | 62% | 87% |
| 0.8 | 0% | 9% | 32% | 79% | 92% |
| **1.0** | 2% | 10% | 54% | 87% | 96% |

At **χ=0 the TIP is inert** — exactly Dodd & de Boer's "no TIP–reservoir interaction, immunity
decides" outcome. Our model departs from theirs **continuously** as χ rises. So the two are not in
conflict: theirs is the χ→0 face of this surface.

## Result 2 — the help is largest at INTERMEDIATE immune strength, and needs real coupling

TIP effect (points above the χ=0 / de Boer baseline):

| χ \ kf | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|
| 0.4 | 0 | −1 | +6 | +5 | +1 |
| 0.6 | +1 | −1 | +4 | **+13** | +2 |
| 0.8 | 0 | +4 | +7 | **+22** | +6 |
| 1.0 | +2 | +6 | **+30** | **+30** | +9 |

Two clean, non-obvious structural facts:
- **The benefit peaks in the "marginal controller" band (kf≈11–12).** Where immunity is too weak
  (kf 9–10) the TIP can't rescue a losing fight; where immunity already wins (kf 13, 87% control)
  there's little left to add (ceiling). The TIP matters **most where control is on a knife-edge** —
  a small, well-timed push decides it. (At kf=12 the effect is a clean monotone dose-response in
  coupling: 0 → +5 → +13 → +22 → +30.)
- **Substantial coupling is required:** the coupling threshold `χ*` for a meaningful (≥+8 pt)
  benefit is ≈0.6–1.0. A loosely-coupled TIP (χ≲0.4) does little. This is a real, falsifiable
  constraint, not a free win.

## Result 3 — the falsifiable prediction (what a lab could test)

> **A TIP engineered to co-reside in / co-reactivate with the latent reservoir (high χ) is
> predicted to improve durable post-treatment control after ART interruption, while a
> non-reservoir-coupled TIP (χ≈0) is predicted to be neutral.** The benefit is predicted to be
> (i) **largest in intermediate / "marginal" controllers**, (ii) **contingent on substantial
> coupling** (χ≳0.6), and (iii) **monotonically increasing** with coupling.

**Direct test:** in a humanized-mouse or NHP ATI model, compare three arms — (a) no TIP, (b) a
standard (non-reservoir-targeting) TIP, (c) a reservoir-co-residing TIP — measuring time-to-rebound
and the fraction maintaining control. The model predicts **(c) > (b) ≈ (a)**, with the gap widest
in animals with intermediate baseline immune control. A null result in arm (c), or (b)≈(c),
would falsify the coupling mechanism.

## What this is, and what's still needed

This is a **modeling hypothesis with a stated boundary and a falsifiable prediction** — the form a
within-host modeler can engage with — not a demonstrated effect. Before claiming priority, the
honest remaining pre-work: (1) a **global sensitivity analysis** (does the χ-help survive joint
variation of reactivation rate, reservoir size, ψ, visibility, exhaustion?); (2) a **semi-analytical
criterion** — the effective reproduction number of a reactivating lineage with vs without the
coupled TIP, to give the *condition* for the sign flip, not just simulation; (3) a **mechanistic
coupling model** (co-packaging / co-integration) replacing the coarse χ knob; and ultimately the
wet-lab test above. All of it inherits the project's standing caveats (illustrative params,
functional-not-sterilizing control, one model family — see [AUDIT2.md](AUDIT2.md)).

_Caveats: stochastic, N=180; χ is a coarse single coupling knob; "control" = active-infection
extinction with the reservoir persisting; the intermediate-immunity peak and χ* threshold are
empirical over this grid, not proven._
