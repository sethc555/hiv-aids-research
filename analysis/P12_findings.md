# P12 — coupled TIP under a wear-down-able immune memory: it still helps (no backfire)

_Cope Labs / Seth · 2026-06-21 · model [tip_model_p12_exhaustible.py](tip_model_p12_exhaustible.py) ·
the last [AUDIT2.md](AUDIT2.md) open item. Closes the harm-channel question P11 left open._

## Why this phase exists

P11 showed a **coupled** TIP (one that rides the reactivating reservoir) *helps* post-treatment
control — but used quasi-steady immunity with a constant antigen **floor**, so CD8 could never be
starved. AUDIT2's remaining worry: with a **wear-down-able immune memory and no floor**, a coupled
TIP that mops up the rebound might *starve* the antigen CD8 needs (the P1 antagonism) and
**backfire**. P12 tests exactly that.

## The build (and why it had to be stochastic)

The QSS+floor immunity is replaced by a **dynamic effector memory `E`** earned from *active*
antigen only and decaying without it (no floor): `dE = ρE·A/(A+KA)·(1−E) − δE·E`, `A = Iw + ν·Id`,
killing `= kclear·E` (and `ν·kclear·E` on the TIP's dual cells). A first **deterministic** attempt
was degenerate — a smooth ODE has one blurred setpoint, so everything relaxes back and the TIP
looks neutral by construction. Control-vs-rebound is genuinely a **stochastic** outcome (the
reservoir dies out or re-establishes), so this is built in the tau-leaping engine. Slow immune
recovery (`ρE` small) puts the cohort in the **marginal regime** (partial baseline control ~72–82%)
where a TIP effect is actually visible; with fast recovery immunity self-rescues to ~98% and
nothing is distinguishable.

## Result — the coupled TIP still HELPS; the backfire does not appear

Control with no TIP (χ=0) vs fully coupled TIP (χ=1), across post-ART memory `e_boost` and TIP
dual-cell visibility `ν` (N=200, 700-day ATI):

| memory \ visibility ν | 1.0 | 0.5 | 0.2 |
|---|---|---|---|
| **0.6** | 76→**89%** | 76→**93%** | 82→86% |
| **0.5** | 72→**91%** | 78→**90%** | 80→86% |
| **0.3** | 72→**89%** | 72→**86%** | 76→**84%** |

**Helps in 7/9 cells, neutral in 2 (both at the most-evasive ν=0.2), harms in none.** Visibility
sets the *magnitude* of the help, not its sign: a visible coupled TIP (ν≥0.5) clearly helps
(+12–18 pts); a strongly-evading one (ν=0.2) gives mild help to neutral (+4–8 pts). The
hypothesized antagonism — TIP suppresses the rebound, starves the memory, control collapses —
**did not materialize** in any tested regime, because even evasive dual cells present enough
antigen (`ν·Id`) and the TIP's suppression of the rebound outweighs the starvation cost.

## What this settles, and what it doesn't

- **Settles:** the P11 "coupled TIP helps" finding is **robust across two structurally different
  immune models** — P11's maintained-floor QSS immunity *and* P12's wear-down-able dynamic memory.
  The structural reversal AUDIT2 feared (on the immunity axis) did **not** happen here.
- **Doesn't settle:** only specific regimes were swept. A more aggressive exhaustion mechanism
  (effectors actively *degraded* by sustained antigen, not merely unmaintained), or visibility
  below 0.2, could in principle still tip it to harm — untested. And, as always, these are
  illustrative parameters; `ρE` was chosen to expose the marginal regime.

## Net (the project's current bottom line)

Putting P11 + P12 together: **a TIP is not orthogonal to the cure, and — when coupled to the
reactivating reservoir — it helps post-treatment control, a result that survives both a
maintained and a wear-down-able immune model.** The earlier "neutral" verdict was an artifact of
modeling the TIP and the reservoir as non-interacting. The remaining caveats are quantitative
(how much coupling is achievable; extreme-exhaustion regimes), not a lurking sign flip.

_Caveats: illustrative stochastic model, N=200; "control" = active-infection extinction with the
reservoir persisting (functional, not sterilizing); coupling χ and visibility ν are coarse single
knobs; only the marginal regime and ν≥0.2 were swept; a degradation-type exhaustion term is
untried._
