# P11 — coupling the TIP to the rebound flips the verdict: a TIP can HELP

_Cope Labs / Seth · 2026-06-21 · model [tip_model_p11_coupled.py](tip_model_p11_coupled.py) ·
the [AUDIT2.md](AUDIT2.md) #1/#2 follow-up. **This phase materially revises the project's headline.**_

## Why this phase exists

AUDIT2 (the multi-agent re-audit) found that the whole "TIP is neutral / orthogonal to the cure"
verdict rested on **one structural choice**: the model made the TIP and the latent reservoir
**non-interacting compartments** — a reactivating provirus (`L_lat`→`Iw`) never passed through a
TIP-accessible co-infected state. But a TIP is *designed* to mobilize on wild-type co-infection,
so a reactivating reservoir virus *should* be able to carry/acquire it. P11 adds exactly one knob:

> **χ = the fraction of reservoir reactivations that emerge as TIP-carrying dual cells** instead
> of pure wild-type cells. χ=0 reproduces the old (decoupled) model; χ=1 means the TIP rides every
> rebound — the engineered TIP co-resides with the reservoir, as intended.

## Result — the verdict flips, monotonically (3-seed confirmed)

At the marginal immune level (kf=11, baseline antigen floor), post-treatment control vs coupling:

| coupling χ | 0.0 | 0.5 | 0.75 | 0.9 | 1.0 |
|---|---|---|---|---|---|
| **P(durable control)** | **22%** | 32% | 40% | 45% | **58%** |

A clean, monotonic **dose-response**: the more the TIP is coupled to the rebound, the more it
**helps** — roughly tripling control (22%→58%) at the immune level where it matters most. At a
stronger immune level (kf=13) the effect is smaller (88%→96%) because immunity already controls.
**Crucially, CD8 retention stayed ~100% throughout** — the coupled TIP suppressed the rebound
*without* starving the immune response. So in this model the TIP acts as a genuine **adjunct**:
it intercepts reactivating virus, immunity stays primed, and durable control improves.

**The earlier "TIP is neutral" was an artifact of the decoupled construction (χ=0).** Once the
TIP can do the one thing it is engineered to do — engage co-infecting wild-type — it stops being
inert and starts to help.

## Two honest limits this same run shows

1. **The TIP assists, it does not replace immunity.** With a weak antigen floor (RDEF=500, i.e.
   immunity not maintained), control was **0% at every χ and every kf** — a coupled TIP cannot
   rescue a collapsed immune response, only amplify a working one.
2. **The help is graded, and partial coupling is weak.** χ=0.5 only lifts control modestly
   (22%→32%); the large gains need high coupling (χ≥0.9). A real engineered TIP would have to
   co-reside in a large fraction of reactivating reservoir cells to deliver the benefit.

## What is still NOT tested (the remaining open channel)

No **antagonism / over-suppression harm** appeared here — but only because this model's immunity
is quasi-steady with a maintained antigen floor, so CD8 *cannot* be starved or exhausted. The P1
antagonism (TIP suppresses WT → less antigen → CD8 wanes) and the cross-disease "over-suppression
hurts" effect require a **dynamic, exhaustible immune memory** model, which is still unbuilt (the
abandoned P9). So P11 establishes the **help** side of a coupled TIP; whether strong coupling also
reopens the **harm** side under exhaustible immunity is the next, genuinely open question.

## Revised project verdict

- **Decoupled (χ=0):** a TIP is neutral to post-treatment control. _(old headline)_
- **Coupled (χ>0, the realistic engineered design):** a TIP **helps**, monotonically in coupling
  strength, *provided immunity is maintained* — up to ~tripling control at the marginal immune
  level. _(P11)_

So the real-world status is **not** "orthogonal." It is: **a TIP's value depends entirely on how
well it couples to the reactivating reservoir** — kept apart it does nothing; engineered to ride
the rebound (and paired with a maintained immune response) it is a useful adjunct. The audit was
right that the neutral verdict was one structural assumption deep — and removing that assumption
changes the answer.

_Caveats: illustrative stochastic model (same Perelson-class kernel, N=200 cohort); "control" =
active-infection extinction with reservoir persisting (functional, not sterilizing); the HELP is
shown, the potential exhaustion-driven HARM is untested (needs the dynamic-immunity successor);
χ is a coarse single coupling knob, not a mechanistic co-packaging model._
