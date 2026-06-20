# hiv-aids-research

A research log on **HIV cure strategy**, combining (1) a citation-grounded literature
corpus of the 2024–2026 cure frontier and (2) a within-host **mathematical-modeling**
investigation of **therapeutic interfering particles (TIPs)** versus CD8 immune control —
with **multi-agent adversarial audits** of every claim.

> ⚠️ **Status: research log, not findings for clinical or wet-lab use.**
> The models are illustrative mean-field ODEs with un-calibrated parameters. The central
> modeling result (a TIP↔CD8 "antagonism") is, by the project's own audit, **partly a
> property of the model's single-antigen-pool construction**, and the model's "immune
> clearance" term is only a **loose analogy** to the clinical "vaccinal effect." Treat
> everything here as hypothesis-generation and methods, audited for honesty — nothing more.

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
  P1_findings.md … P1.5_findings.md   per-phase write-ups
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
4. **The one genuinely emergent (non-tautological) result**: P1.1's **early-TIP immune-
   blunting** — a dynamically-driven failure mode, not pinned by the single-pool construction.
5. **Real-world anchor**: Simonetti's 5′-leader-defective non-suppressible viremia
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
