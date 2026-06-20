# P4 — the latent reservoir + ATI: closing the model↔clinic gap; the TIP is neutral to cure

_Cope Labs / Seth · 2026-06-20 · model [tip_model_p4_reservoir.py](tip_model_p4_reservoir.py) ·
figure [p4_ati.png](p4_ati.png) · the caveat that shadowed every prior phase._

## Why this phase exists

The audit's standing HIGH caveat ([AUDIT.md](AUDIT.md)): every model "controlled" *active*
infection, but the clinical cure target is the **latent replication-competent reservoir**,
which no prior model contained — so "TIPs antagonize the 2024–26 immune-control breakthroughs"
generalized an active-cell result to a reservoir mechanism that wasn't represented. P4 builds
the reservoir and runs the actual clinical experiment.

## The build

A 7-compartment stochastic (tau-leaping) model with **two biologically distinct reservoirs**:
`L_def` (the Simonetti defective antigen clone — primes CD8, makes no virus, TIP-proof) and
`L_lat` (a **latent replication-competent** provirus pool — transcriptionally **silent**, hence
*absent from the CD8 antigen pool*, long-lived, reactivating stochastically to a productive WT
cell). An engineered TIP can be delivered two ways: as a persistent latent dual pool
(`L_latd`, rides through ART) or as a bolus at treatment interruption. The schedule is the
clinical one: **CHRONIC → ART (block new infection) → ATI (release, watch for rebound)**.
"Post-treatment control" = productive cells (Iw+Id) < 50 at the end of a 500-day ATI; the CD8
"vaccinal-effect" strength is the knob `kf`.

## Result 1 — post-treatment control is a sharp CD8 threshold; the reservoir always persists

Control is **all-or-nothing in CD8 strength**: P(control) jumps 1% → 59% → 100% across
kf = 10 → 12 → 14 (sharp transition at kf ≈ 12, where CD8 killing pushes the rebound's
effective R₀ below 1). **In every controlled replicate the latent reservoir persists**
(~200–470 latent cells retained) — i.e. **functional, not sterilizing, cure**. This is the
model↔clinic gap closed concretely: "control" now means *durable suppression of a persisting
reservoir's rebound*, exactly the post-treatment-control / RIO / Geneva setting — not
active-infection extinction.

## Result 2 (headline) — a TIP is neutral to post-treatment control

The control curve is **unchanged across all TIP arms** (80 reps/cell):

| kf | no TIP | TIP latent/engineered | TIP bolus@ATI (ν=0.9) | TIP bolus@ATI (ν=0.3) |
|---|---|---|---|---|
| 10 | 1% | 1% | 2% | 1% |
| 12 | 59% | 57% | 56% | 57% |
| 14 | 100% | 100% | 100% | 99% |
| ≥16 | 100% | 100% | 100% | 100% |

The threshold's position is set by CD8 strength alone; the TIP — latent or bolus, visible
(ν=0.9) or evasive (ν=0.3) — **does not move it**. A stress test with an unrealistically large
*co-reactivating* engineered TIP reservoir nudges the margin **down** (kf=12: 61% → 51%, within
sampling noise), never up: at most a hint of harm, **never a benefit**.

## The mechanism — TIP and immune-control cure are misaligned by construction

The reason is structural and is the project's deepest point: **the TIP is parasitic on active
WT replication, but the immune-control cure works by eliminating active WT.** So:

- **In the control regime** the TIP has *no substrate* — WT is suppressed, dual cells never
  form (Vt ≈ 0, TIP extinct in nearly all controlled reps). The TIP cannot *help* control
  because control is precisely the absence of the WT it needs.
- **In the rebound regime** the TIP is *mistimed*: bolus carrier cells die in ~1 day, before
  the slow reservoir-driven rebound is large enough to amplify them; the latent TIP reactivates
  too rarely to catch the burst. The TIP cannot *stop* a rebound CD8 couldn't.

So the TIP is orthogonal to reservoir control — it lives on the active-infection axis the cure
is trying to zero out.

## What this does to the project's central claim

It **resolves the audit's HIGH caveat in the TIP's favor, and against the strong framing
both ways.** With the reservoir properly modeled: the active-infection "antagonism" of P1–P2
does **not** translate into harm to the immune-control cure strategy — a TIP is **neutral** to
post-treatment control (at most marginally harmful under a forced large dose, never helpful).
The headline "TIPs antagonize the 2024–26 breakthroughs" is **not supported** once the cure
target is in the model; but neither is "TIP as a useful adjunct to immune-control cure" — the
TIP is simply **inert to reservoir control**, because it needs the very thing the cure removes.

## Net

Building the latent reservoir closes the loose analogy that shadowed P1–P3 and reframes the
verdict: on the active-infection axis a TIP and CD8 trade off (the antagonism, P1–P2,
probabilistic in P3), but on the **reservoir-control axis that defines cure**, the TIP is
**neutral** — it has no substrate where the cure succeeds and can't catch the rebound where it
fails. A TIP is neither a threat to, nor a tool for, the durable-remission cure strategy in
this model.

_Caveats: illustrative Perelson-class + illustrative latency/reactivation rates (F_lat=5e-4,
reactivation 1e-3/day); control is modeled as CD8 killing driving rebound R₀<1 (one of several
real control mechanisms); "control" = productive-cell suppression with the reservoir persisting
(functional cure); the slight-harm signal at the margin is sub-significant; stochastic
tau-leaping (the P3c-validated engine), single (reservoir-size, ψ) point._
