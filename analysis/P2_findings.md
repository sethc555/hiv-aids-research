# P2 — the two-pool (decoupled-antigen) test: the antagonism is *not* merely by-construction

_Cope Labs / Seth · 2026-06-20 · model [tip_model_p2_twopool.py](tip_model_p2_twopool.py) ·
figure [p2_escape.png](p2_escape.png) · the multi-agent audit's identified next step._

## Why this phase exists

The multi-agent audit's load-bearing finding ([AUDIT.md](AUDIT.md) point #1) was that the
TIP↔CD8 antagonism is **partly by-construction**: P1–P1.5 prime CD8 from a **single antigen
pool** proportional to active infection (`A = Iw + ν·Id [+ static A_def]`), so "TIP lowers
active infection → starves antigen-driven CD8" is semi-definitional. The audit downgraded
the claim from "discovered" to "encoded, then probed," and concluded the real-world status
is **genuinely open**. P1.5 carried a *static* decoupled offset `A_def` but could not
separate the two candidate explanations for the absence of any stable immune-compatible TIP:

- **(i) construction** — the antagonism is the single-pool assumption; decouple the pool and
  it should disappear; vs
- **(ii) immune-fragility** — the TIP lives at R0 ≈ 1, so any CD8 that kills its
  dual-infected factories (`kd = ν·k·E`) sinks it; this would **survive** decoupling.

## The build

P2 replaces P1.5's static `A_def` with a **dynamic, self-maintaining latent/defective
reservoir** `L` — the strongest honest decoupling, grounded in Simonetti (JCI 2023 / Nat
Commun 2026): ~10⁶–10⁷-cell defective clones that persist by **antigen-independent clonal
proliferation** and present HLA antigen **without producing infectious virus**. In the model
`L`:

- **primes CD8 independently of active WT** — `Aeff = Iw + ν·Id + L` (two pools), and
- is **TIP-independent** (a TIP interferes with WT *packaging* in dual cells; it does nothing
  to a defective clone that only transcribes antigen), and is **CD8-resistant**
  (proliferation-maintained, never killed) — i.e. the *best case for the TIP*: CD8 gets a
  persistent priming source the TIP can neither touch nor starve.

`dL = g·L·(1 − L/R) + s·(Iw+Id)`. The reservoir size `R` and the dual-cell visibility `ν` are
swept at ψ=22, on the same 4-cell QSS competition (WT R0 = 8.70 preserved) and the **same
audit discipline as P1.5** (tail max/min reported; no "escape" counted unless converged).

The decisive extra lever P2 isolates is reservoir **timescale**: a large, slow antigen floor
should pin CD8 ≈ constant and **damp** the P1.4/P1.5 immune-feedback oscillation. So I ran
the grid twice — `g_slow = 0.02/day` (half-life ~months, genuinely decoupled) vs
`g_fast = 2/day` (tracks instantly ≈ P1.5's static `A_def`, the control).

## Result — decoupling does **not** rescue the TIP (explanation (ii) wins)

**Slow and fast reservoirs give the same answer: 0/169 stable immune-compatible escape**
(WT ↓ ≥0.5 log AND CD8 ≥70% AND converged); max WT reduction over **stable** cells = **0.21
log** in both. The slow/fast equivalence is itself a finding — the reservoir timescale does
**not** damp the instability, refuting the sub-hypothesis that a slow antigen floor would
stabilize the escape regime.

The grid slices show the mechanism cleanly (slow reservoir, ψ=22):

| regime | WT ↓ | CD8 kept | stability |
|---|---|---|---|
| ν=1, R=0 (no evasion, no reservoir) | 0.17 log | 86% | **stable, TIP inert** |
| ν=1, R≥0.5×A0 (reservoir floors CD8) | ≈0 / slightly **negative** | 100%+ | oscillating, TIP inert |
| ν→0, R=0 (dual cells evade CD8) | 0.55–0.74 log | **38–52%** | oscillating |

- Where the reservoir **keeps CD8 primed** (large `R`), that primed CD8 **kills the TIP's
  dual-infected factories** → the TIP can't establish (inert; WT ↓ ≈ 0, sometimes <0 as the
  TIP slightly *helps* WT via carrier competition). The decoupled pool doesn't free the TIP —
  it **arms the immune response against it**.
- The only TIP-effective regime is **dual-cell evasion** (ν<1, R=0) — and there the TIP's
  suppression of active WT *does* starve the active pool, **crashing CD8 to 38–52%** (the
  original antagonism, now isolated to the pool the TIP can actually cut) — **and it
  oscillates**.

So the antagonism is **not** merely the single-pool construction. Giving CD8 a fully
decoupled, persistent, TIP-proof priming source — the most generous decoupling available —
still admits no stable immune-compatible TIP, because the obstruction is **mechanical**: the
TIP's R0 ≈ 1 fragility means the very CD8 the cure strategy needs kills the TIP's production
cells. **Two pools, same wall.**

## What this does to the audit's open question

The audit left the real-world status "genuinely open" *on the grounds that the antagonism
was encoded by the single pool*. P2 removes that ground: **decoupling the encoding does not
open an escape.** The open question therefore narrows from the broad "is it just
construction?" to a single sharp, mechanistic condition —

> **does a stable regime exist in which the TIP's carrier/dual cells evade CD8 (ν<1) while
> the rest of the immune response still controls WT?**

Mean-field ODEs answer *no* (P1.4, P1.5, and now P2 with decoupled antigen all oscillate
exactly there). That is precisely the regime where discreteness matters — the Simonetti
clones are individuals, not a mean field — so it **sharpens the motivation for the
stochastic/agent-based successor** rather than another ODE reduction.

## Net

P2 converts the audit's biggest caveat from a standing doubt into a tested statement: the
TIP↔CD8 antagonism **survives full antigen decoupling**, so it is driven by the TIP's
intrinsic immune-fragility (R0 ≈ 1), not by the single-pool construction. The escape
question is now pinned to one falsifiable condition (stable ν<1 carrier evasion) that
mean-field cannot deliver — handing a **better-posed** problem to the stochastic/ABM tool.

_Caveats: mean-field ODE, illustrative params; `L` is modeled CD8-resistant and TIP-proof
(the generous case — a CD8-killable reservoir would only strengthen the antagonism); the
0/169 result is empirical over this grid at ψ=22, not a proven theorem; reservoir antigen
enters CD8 priming linearly (`Aeff = Iw + ν·Id + L`) — a saturating cross-pool presentation
is untried and is the one remaining mean-field variation before the successor model._
