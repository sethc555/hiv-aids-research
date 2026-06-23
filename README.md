# hiv-aids-research

> 📄 **Read the manuscript:** [**rendered web version**](https://sethc555.github.io/hiv-aids-research/) · [PDF](docs/manuscript.pdf)
> · **preprint DOI** [10.5281/zenodo.20801699](https://doi.org/10.5281/zenodo.20801699) · code archive DOI [10.5281/zenodo.20799761](https://doi.org/10.5281/zenodo.20799761)
>
> **Final result (post-audit):** a therapeutic interfering particle (TIP) *coupled* to latent-reservoir
> reactivation is predicted to **assist** immune-mediated HIV-1 post-treatment control — it **extends, does
> not contradict,** Dodd & de Boer (2025) (a non-coupled TIP is neutral and recovers their limit).
> *Illustrative modeling hypothesis — not validated findings, not a cure.* The P2–P6 narrative just below
> is the **earlier, superseded "neutral" stage**; see [P11](analysis/P11_findings.md)–[P16](analysis/P16_findings.md)
> and [paper/MANUSCRIPT.md](paper/MANUSCRIPT.md) for the final verdict. [AUDIT3](analysis/AUDIT3.md) adds the
> contested-biology review: the immune axis is **effector-agnostic** (CD8/NK/antibody), not CD8-specific, and
> "no systematic backfire" replaces the earlier absolute "never backfires."

A research log on **HIV cure strategy**, combining (1) a citation-grounded literature
corpus of the 2024–2026 cure frontier and (2) a within-host **mathematical-modeling**
investigation of **therapeutic interfering particles (TIPs)** versus cytotoxic (CD8/NK) immune control —
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
>
> **VERDICT REVISED ([P11](analysis/P11_findings.md) + [AUDIT2.md](analysis/AUDIT2.md), 2026-06-21):**
> the "neutral / orthogonal" headline was an **artifact of a decoupled model**. The re-audit found
> the verdict hinged on the TIP and the reservoir being modeled as **non-interacting** (a
> reactivating provirus never met the TIP — yet a TIP is *designed* to mobilize on WT co-infection).
> P11 added one knob — χ = the fraction of reservoir reactivations that carry the TIP — and the
> answer **flips**: as χ goes 0→1, post-treatment control rises monotonically (~22%→58% at the
> marginal immune level), **CD8 stays primed**, i.e. a coupled TIP **HELPS** as an adjunct.
> The feared *harm* channel — a coupled TIP starving an exhaustible immune memory (the P1
> antagonism) — was then tested directly ([P12](analysis/P12_findings.md)): with a wear-down-able
> dynamic immune memory (no floor), the coupled TIP **still helps and did not backfire** (helps in
> 7/9 regimes, mild/neutral only at extreme evasion). So "coupled TIP helps" is robust across two
> structurally different immune models. The last corner — **active immune exhaustion** (burden
> *degrades* CD8) plus a **stealthy** TIP — was then tested too ([P13](analysis/P13_findings.md)):
> still **no backfire anywhere**; the TIP in fact helps *most* there (+~50 pts), because by holding
> the viral burden down it **protects the immune response from exhaustion**. So "coupled TIP helps,
> no systematic backfire" now holds across **three** structurally different immune models. Caveats: it
> needs substantial coupling (χ≈0.5 is weak) and **assists but cannot replace** immunity (0%
> control if immunity isn't maintained). **Bottom line: a TIP is not orthogonal to the cure;
> coupled to the reactivating reservoir it is a useful adjunct — most valuable, not least, when the
> immune system is fragile.**

---

## The arc (how this came to exist)

It started from a different question — *could a federated-database / drug-repurposing
engine help cure HIV?* The answer was **no** (the cure-relevant data isn't in federatable
databases), which led to mapping where the field actually is: a decisive shift toward
**durable ART-free remission via immune control of the reservoir** — a multifactorial mechanism (CD8,
NK, and antibody/"vaccinal-effect" effectors; the dominant one is debated and varies by setting — see
[AUDIT3](analysis/AUDIT3.md)) — evidenced by the Geneva & 2nd-Berlin transplants, the RIO bNAb trial,
and the UCSF combination-immunotherapy result. The one
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
  tip_model_p11_coupled.py       P11  — couple TIP to the rebound (audit follow-up): verdict FLIPS to "helps"
  tip_model_p12_exhaustible.py   P12  — coupled TIP under wear-down-able immunity: still helps, no backfire
  tip_model_p13_exhaustion_damage.py  P13 — harshest (active exhaustion + stealth): no backfire, TIP helps more
  tip_model_p14_coupling_phase.py P14  — coupling phase diagram + reduction to de Boer + falsifiable prediction
  tip_model_p15_sensitivity.py   P15  — global sensitivity: benefit robust across joint params, 0% backfire
  tip_model_p16_analytic.py      P16  — semi-analytical criterion: R_eff=R0*d/(d+k); coupling shifts threshold
  tip_model_deboer_bridge.py     recreate Dodd & de Boer 2025 + bridge: why ours extends (not flips) theirs
  METHODS.md                     consolidated methods (equations, params, numerics, audit)
  NOVELTY.md                     prior-art / novelty assessment + the de Boer reconciliation
  P1_findings.md … P6_findings.md, P8_boosters.md   per-phase write-ups
  AUDIT.md                       the self-audit + two multi-agent audit overlays
  *.png, *.npz, raw_cache.json   figures, sweep outputs, cached S2 responses
NARRATIVE.md                     a readable narrative of the build session
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
   curves coincide)** — the **patient-level confirmation** (conditional on the structural caveats
   in [AUDIT2.md](analysis/AUDIT2.md)) that a TIP is neutral-to-the-cure / orthogonal to the
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

## Verifying the claim chain

Every load-bearing number in this README is **machine-checkable**. Run:

```bash
cd analysis && python3 verify_claims.py     # 14/14 checks, exit 0
```

It independently re-derives each headline claim from the phase modules and asserts it within a
tolerance — R0=8.70, the TIP invasion threshold (~7.5) and immunity-collapse (P1), P3's 0%
stable coexistence / ~order-15% cure / ~82% CD8, P4's sharp control threshold and TIP-neutrality,
P6's spontaneous-control fraction (~5%) and clinical rebound median (~20 d), and P8's trough rule
and reservoir×immunity synergy. A regression in any script trips a FAIL. Stochastic checks use
reduced replicate counts (direction/order, not last digit); the precise figures live in the
findings docs. Claims the audit *retracted* (oscillation-fragile P1 magnitudes) are deliberately
**not** asserted — the harness checks the robust versions the audit upholds.

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

Built in a single interactive session (2026-06-19/20). Parameters are standard-but-illustrative (Perelson-class within-host HIV);
no calibration, no clinical validation. The work is **defensive/basic research** on a
published cure concept — read the audit before citing any number.

## License

Dual-licensed, intentionally:
- **Code** (the `analysis/` scripts, `verify_claims.py`, build tooling) — **MIT** (see [LICENSE](LICENSE)).
- **Text, figures, and prose** (`paper/`, `docs/`, the manuscript, `ABSTRACT.md`, `*_findings.md`,
  this README) — **CC BY 4.0** (https://creativecommons.org/licenses/by/4.0/).

This is the standard split for a repository that is both software and a manuscript (MIT is a
software license and is a poor fit for prose; CC BY is the right fit for the text). The Zenodo
deposit (DOI `10.5281/zenodo.20799761`) is labeled MIT as the archive-level license; this section
is the authoritative statement of the per-component licensing.
