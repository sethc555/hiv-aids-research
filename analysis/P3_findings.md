# P3 — the stochastic successor: the escape is a *probability*, not a fixed point

_Cope Labs / Seth · 2026-06-20 · model [tip_model_p3_stochastic.py](tip_model_p3_stochastic.py) ·
figure [p3_outcomes.png](p3_outcomes.png) · the tool the whole ODE program pointed to._

## Why this phase exists

Every mean-field result — P1–P1.5, P2 (antigen decoupling), P2b (Holling-II killing) —
hit the **same wall**: the only regime where the TIP suppresses WT is dual-cell immune
evasion (ν<1), and there the ODEs **oscillate**, so the question "does an immune-compatible
TIP exist?" was **undecidable in the continuum** (no fixed point to read off). P1.5 and P2
both named the successor explicitly: a **stochastic** model, where (a) the limit-cycle
troughs the ODE rides through become real **extinction** events, and (b) the Simonetti
defective clones are **discrete individuals**, not a mean field.

## The build

A demographic **tau-leaping** version of the P2 cell model: the same 5 integer populations
(T, Iw, It, Id, L), QSS free virus, two-pool antigen (`Aeff = Iw + ν·Id + L`), slow
self-maintaining reservoir. 200 replicates per condition, run from the WT+reservoir setpoint
after a TIP bolus, classified by **absorbing state**. Audit discipline: report the full
outcome **distribution** (never a single run); the headline was confirmed **dt-robust**
(16% vs 17% at dt = 0.03 vs 0.015).

## Result — the oscillation was a continuum artifact; the answer is a lottery

The mean-field limit cycle **does not survive discreteness**. Stochastically there is **no
stable coexistence (0% at every ν)** — every replicate falls into one of two absorbing
outcomes, a **lottery** between them (reservoir R = 5000 ≈ 1×A0, ψ = 22):

| ν (dual-cell visibility) | WT controlled ("cure") | TIP lost, WT rebounds | coexist | CD8 kept (in controlled) |
|---|---|---|---|---|
| 0.2 (strong evasion) | **16%** | 84% | 0% | **82%** |
| 0.4 | 2% | 98% | 0% | — |
| 0.6–1.0 | 0% | 100% | 0% | — |

Three things the ODEs could not have told us:

1. **The immune-compatible TIP outcome genuinely EXISTS** — in the WT-controlled replicates,
   CD8 is retained at **82%**. When the TIP works, it *is* immune-compatible. The continuum
   could never show this because it never reached an absorbing state.
2. **But it is a low-probability event, not a design point** — ≤16%, and only under **strong
   immune evasion** (ν ≤ 0.2); by ν = 0.4 it is ~2%, and it vanishes without evasion. You do
   not *design* an immune-compatible TIP; you **win** one, against odds, in a regime
   (immune-evading carrier cells) that is itself the risky property.
3. **The reservoir is causal.** Re-running ν = 0.2 with **R = 0 gives 0% control**; R = 5000
   gives 16% control with CD8 at 82%. The persistent, TIP-proof reservoir antigen (P2's
   decoupled pool) is what keeps CD8 **armed independently of the TIP's suppression**, long
   enough to drive active WT extinct in a stochastic trough. **The reservoir converts the
   antagonism into a low-probability *synergy*** — the most hopeful, and most honest, result
   in the project.

## The reframing this buys

The program's central question — *can a TIP suppress the virus without starving the immune
response that controls it?* — is now **answered**, and the answer changes shape: it is not a
parameter regime you tune into (the ODE framing), it is a **probability you have to win**.
That probability is (i) small, (ii) **increasing with immune evasion** — anti-correlated
with TIP safety — and (iii) **conditional on a persistent reservoir antigen** that keeps CD8
primed. The right design objective is therefore not "find the stable escape" but "**raise
P(cure) / P(TIP-loss)**" — a fundamentally stochastic, optimal-control problem.

## Net

The stochastic successor did exactly what P1.5/P2 said it would: it dissolved the
mean-field oscillation (a continuum artifact) into a quantified outcome distribution, and
turned an *undecidable* question into a *measured probability*. The immune-compatible TIP is
real but rare, evasion-dependent, and reservoir-enabled.

_Update: the phase map + validation ([P3b/P3c](P3b_findings.md)) firm the cure estimate to
**order-15% (≈10–20% across dt and an exact reduced-scale SSA)** and confirm the invariants
(0% coexistence, CD8 ~82%) are method-independent; the qualitative headline below stands._

_Caveats: tau-leaping (dt-robust at two step sizes, but not an exact SSA); illustrative
Perelson-class params; **"WT controlled" = extinction of active-infection compartments**, NOT
eradication of a replication-competent latent reservoir (the model still has none — the same
model↔clinic loose-analogy caveat as [AUDIT.md](AUDIT.md)); a single (R, ψ) point was
profiled in depth — the full stochastic phase map over (R, ν, ψ) and an exact-SSA
spot-check are the next grind._
