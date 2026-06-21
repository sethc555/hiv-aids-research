# hiv-aids-research

A research log on **HIV cure strategy**, combining (1) a citation-grounded literature
corpus of the 2024–2026 cure frontier and (2) a within-host **mathematical-modeling**
investigation of **therapeutic interfering particles (TIPs)** versus CD8 immune control —
with **multi-agent adversarial audits** of every claim.

> ⚠️ **Status: research log, not findings for clinical or wet-lab use.**
> The models are illustrative mean-field ODEs (plus a tau-leaping stochastic successor) with
> un-calibrated parameters. The central modeling result (a TIP↔CD8 "antagonism") was probed
> for the project's own "by-construction" caveat ([P2](analysis/P2_findings.md): it survives
> antigen-pool decoupling, so it is driven by the TIP's R0≈1 **immune-fragility**, not the
> single-pool construction), then re-posed stochastically ([P3](analysis/P3_findings.md): the
> immune-compatible TIP is a **low-probability outcome**, not a stable design point), then
> tested against the actual cure target ([P4](analysis/P4_findings.md): with a
> replication-competent latent reservoir + ART/ATI, the TIP is **neutral to post-treatment
> control**), and finally **calibrated to clinical ATI data** ([P5](analysis/P5_findings.md):
> the rebound clock matches the ~22-day clinical median and the neutrality survives anchoring to
> the RIO/CHAMP control fractions), and finally fit at the **patient-population level**
> ([P6](analysis/P6_findings.md): a heterogeneous cohort reproduces the whole ATI rebound
> Kaplan-Meier curve, and the TIP still moves post-treatment control by ~1 point). Methods are
> consolidated in [METHODS.md](analysis/METHODS.md). The work is a calibrated research log —
> treat it as hypothesis-generation and methods, audited for honesty — nothing more.

---

## The arc (how this came to exist)

It started from a different question — *could a federated-database / drug-repurposing
engine help cure HIV?* The answer was **no** (the cure-relevant data isn't in federatable
databases), which led to mapping where the field actually is: a decisive shift toward
**durable ART-free remission via CD8/"vaccinal-effect" immune control** (the Geneva & 2nd-
Berlin transplants, the RIO bNAb trial, the UCSF combination-immunotherapy result). The one
**no-wet-lab modeling lane** that was genuinely open was **TIPs** — engineered conditionally-
replicating defective HIV genomes (Weinberger/Pitchai; de Boer). So the project became a
multi-phase model of *whether a TIP helps or fights* an immune-control cure, stress-tested by
audits at every step.

## What's here

```
analysis/          all working files (flat, runnable in place)
  scan.py, analyze.py            Semantic Scholar corpus builder + citation-graph expander
  scan_results.md, bibliography.md, citation_graph.md, MEMO.md   the literature corpus
  tip_model.py                   P1   — static immunity vs TIP (the base model)
  tip_model_dynamic.py           P1.1 — dynamic immune feedback
  tip_model_p12.py, p12_pareto.py  P1.2 — memory + exhaustion immunity
  robustness.py, check_*.py      audit verifications (robustness sweep, controls)
  tip_model_p13_wm.py / _spatial.py / p13_check.py   P1.3 — stabilized + spatial test
  tip_model_p14.py, p14_check.py P1.4 — Simonetti-grounded escape test
  tip_model_p15_stable.py        P1.5 — stable-by-construction build
  tip_model_p2_twopool.py        P2   — two-pool decoupled-antigen test (the audit's #1 caveat)
  tip_model_p2b_satkill.py       P2b  — Holling-II saturating-predation control (closes mean-field)
  tip_model_p3_stochastic.py     P3   — stochastic (tau-leaping) successor: the escape as a probability
  tip_model_p3b_phasemap.py      P3b  — stochastic P(cure)/P(TIP-loss) phase map over (R, nu)
  tip_model_p3c_ssa.py           P3c  — validation: dt->0 convergence + exact reduced-scale SSA
  tip_model_p4_reservoir.py      P4   — latent replication-competent reservoir + ART/ATI (model<->clinic)
  tip_model_p4b_robustness.py    P4b  — robustness: sustained/high-psi/engaged TIP still neutral
  tip_model_p4c_calibration.py   P4c  — calibration: ATI rebound timing vs clinical anchors
  tip_model_p5_calibration.py    P5   — fit to clinical ATI data; TIP verdict survives anchoring
  tip_model_p6_heterogeneous.py  P6   — patient population (kf+reservoir distribution); KM-curve fit
  tip_model_p8_boosters.py       P8   — booster vaccines: titration rule + reservoir×immunity synergy
  METHODS.md                     consolidated methods (equations, params, numerics, audit)
  P1_findings.md … P6_findings.md, P8_boosters.md   per-phase write-ups
  AUDIT.md                       the self-audit + two multi-agent audit overlays
  *.png, *.npz, raw_cache.json   figures, sweep outputs, cached S2 responses
conversation/
  session-transcript.jsonl       the full build-session thread (raw)
NARRATIVE.md                     the same thread as a readable story
```

## Key results (audit-calibrated)

1. **Literature map** ([MEMO.md](analysis/MEMO.md), [citation_graph.md](analysis/citation_graph.md)):
   the cure field has shifted toward immune-mediated remission; the active TIP-modeling
   community is small and theory-first (Weinberger/Singh/Dar/de Boer).
2. **The TIP↔CD8 antagonism** ([P1–P1.3](analysis/P1.3_findings.md)): in the simplest
   well-mixed model a TIP can only suppress wild-type by removing the infected-cell antigen
   pool that primes CD8 — so suppression and CD8 maintenance trade off. **Audit caveat:**
   this is *partly by-construction* (single shared antigen pool), not a discovered law.
3. **It does not survive realistic refinements cleanly** ([P1.4](analysis/P1.4_findings.md),
   [P1.5](analysis/P1.5_findings.md)): dynamic feedback, spatial structure, and immune-
   evasion all **re-introduce oscillation**; the apparent "escapes" were transient artifacts
   (caught by convergence checks). The stable-by-construction build (P1.5) shows
   *immune-pressure + TIP-evasion + a stable fixed point are mutually incompatible* in these
   mean-field ODEs — so the escape question needs a **stochastic/agent-based** successor.
4. **The antagonism is *not* merely by-construction** ([P2](analysis/P2_findings.md)): the
   audit's central caveat was that a single antigen pool makes the antagonism semi-definitional.
   P2 decouples it — a dynamic, self-maintaining, TIP-proof reservoir antigen pool (Simonetti
   clones) primes CD8 independently of active WT. Even so, **0/169 stable immune-compatible
   TIP**: a primed CD8 pool kills the TIP's dual-infected factories regardless of where the
   antigen came from. The obstruction is the TIP's R0 ≈ 1 **immune-fragility**, not the
   construction — which *narrows* the open question to one falsifiable condition (stable ν<1
   carrier evasion) for the stochastic/ABM successor.
5. **The escape is a *probability*, not a fixed point** ([P3](analysis/P3_findings.md)): the
   stochastic (tau-leaping) successor dissolves the mean-field oscillation — it was a continuum
   artifact (**0% stable coexistence**). Outcomes are a **lottery** between WT control ("cure")
   and TIP loss. An immune-compatible TIP **does exist** (CD8 kept at 82% when it works) but is
   a **low-probability event (~16%)**, only under strong immune evasion, and **only with the
   reservoir present** (R=0 → 0%) — the reservoir converts the antagonism into a rare *synergy*.
   The design objective shifts from "find the stable escape" to "raise P(cure)/P(TIP-loss)."
   The [P3b](analysis/P3b_findings.md) phase map makes this a surface — cure needs **both**
   strong evasion (ν≲0.25) **and** a reservoir (R≳1×A0), peaking ~31%, zero coexistence
   anywhere — and [P3c](analysis/P3b_findings.md) validates it (dt→0 + an exact reduced-scale
   SSA: invariants method-independent; cure magnitude order-15%).
6. **With the latent reservoir modeled, a TIP is *neutral* to the actual cure**
   ([P4](analysis/P4_findings.md)): the standing model↔clinic caveat is closed by adding a
   replication-competent latent reservoir + an ART→ATI (treatment-interruption) schedule.
   Post-treatment control is a **sharp CD8-strength threshold**, and the reservoir **persists
   in every controlled case** (functional, not sterilizing, cure). A TIP — latent or bolus,
   visible or evasive — **does not move that threshold** (at most a sub-significant nudge
   *down*, never up). Mechanism: the TIP is parasitic on active WT, but the cure works by
   *removing* active WT, so the TIP has no substrate where the cure succeeds and can't catch
   the rebound where it fails. So the active-infection "antagonism" does **not** translate
   into harm to the immune-control cure — the TIP is **inert to reservoir control**.
   *Robustness ([P4b](analysis/P4b_findings.md)):* the verdict holds even when the TIP is
   continuously dosed, high-ψ, and genuinely engaged (threshold shift within noise); the only
   regime where a TIP suppresses rebound is strong-evasive-continuously-dosed — TIP-dependent
   chronic suppression (a "second ART"), not a durable cure. *Calibration ([P4c](analysis/P4b_findings.md)):*
   the model's ATI rebound clock is clinically realistic **unfitted** — median 7.4 d at the
   macaque-control setting, 10–21 d at the human-ATI setting, abolished above threshold (RIO/bNAb).
7. **Calibrated to clinical ATI data, the verdict holds** ([P5](analysis/P5_findings.md)): the
   rebound clock fits A5345's ~22-day median, and at the operating points reproducing clinical
   post-treatment-control fractions (placebo ~5%, bNAb/RIO ~24%) a sustained, engaged TIP shifts
   durable control by only **+2–5 points** — P4b neutrality is anchored, not just illustrated.
   The fit also **exposes a model limit**: no single immune setpoint reproduces both the
   fast-rebound majority and the controller minority, so quantitative cohort fitting needs
   **population heterogeneity** — the cleanly-identified next step. Methods consolidated in
   [METHODS.md](analysis/METHODS.md).
   *Heterogeneous fit ([P6](analysis/P6_findings.md)):* a patient population `kf ~ Normal(10, 2.5)`
   + a reservoir-clock distribution closes the whole rebound **Kaplan-Meier curve** at once
   (median 20.8 d, wk12 78%, **PTC 5%** — matching A5345/ACTG/CHAMP), resolving the P5 tension.
   On this clinically-fit cohort the TIP moves post-treatment control by **−1 point (5%→4%, KM
   curves coincide)** — the **final, patient-level confirmation** that a TIP is orthogonal to the
   immune-control cure.
8. **The one genuinely emergent (non-tautological) result**: P1.1's **early-TIP immune-
   blunting** — a dynamically-driven failure mode, not pinned by the single-pool construction.
9. **The positive mirror — what *would* help** ([P8](analysis/P8_boosters.md)): unlike the TIP,
   a booster vaccine/antibody acts **directly** on immune control, so it helps. The model gives a
   titration rule: a **single boost only delays rebound**; sustained dosing **banks treatment-free
   cures**; control tracks the **trough not the peak** (maintain ≈1.5× baseline immunity, interval
   set by immunity's half-life); and the real room is the **reservoir×immunity synergy** — neither
   lever alone exceeds ~30–40%, but reservoir reduction + boosting together reach **~70–96%** cure.
10. **Real-world anchor**: Simonetti's 5′-leader-defective non-suppressible viremia
   (JCI 2023 / Nat Commun 2026) — nature makes the TIP's packaged-but-Env-deficient payload
   in huge clones; see the [P1.3 addendum](analysis/P1.3_findings.md).

## Auditing

Every claim was put through adversarial review — see [AUDIT.md](analysis/AUDIT.md):
a same-day self-audit, then a multi-agent audit (math/literature/logic/adversary, each
independently verified), then a sanitized re-run of the dimensions an automated content
filter had blocked. The audits **retracted real overclaims** (a fabricated quote, a
mis-derived R0 scaling, a rehabilitated statistic) and **corrected the project's own
audit** (a wrong self-critique). The literature facts were independently verified; the math
was fully re-derived. This honesty trail is the most reusable part of the repo.

## Running it

```bash
pip install -r requirements.txt
cd analysis
python3 tip_model.py                 # P1 base model + phase diagram
python3 tip_model_p15_stable.py      # P1.5 stable build + escape test
python3 tip_model_p2_twopool.py      # P2 two-pool decoupled-antigen test
python3 tip_model_p3_stochastic.py   # P3 stochastic successor (outcome distribution)
# corpus rebuild needs a Semantic Scholar API key in $S2_API_KEY:
S2_API_KEY=... python3 scan.py
```
Models are plain numpy/scipy/matplotlib. `raw_cache.json` lets the corpus analysis run
without re-querying Semantic Scholar.

## Provenance & caveats

Built in a single interactive session (2026-06-19/20); the full thread is in
`conversation/`. Parameters are standard-but-illustrative (Perelson-class within-host HIV);
no calibration, no clinical validation. The work is **defensive/basic research** on a
published cure concept — read the audit before citing any number.
