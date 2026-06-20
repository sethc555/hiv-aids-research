# HIV cure — state of the field & where a (no-lab) modeler fits

_Synthesis of two research passes (web deep-research + Semantic Scholar citation
scan, 186→275 papers) and a citation-graph expansion of 5 anchor papers.
Cope Labs / Seth · 2026-06-19._

Companion data: [scan_results.md](scan_results.md) (raw 18-topic scan),
[bibliography.md](bibliography.md) (filtered foundation doc),
[citation_graph.md](citation_graph.md) (anchor references/citations + active labs).

---

## Bottom line

There is **no scalable cure**; the latent reservoir remains the barrier and the
field has decisively reframed the goal from *eradication* (sterilizing cure) to
*durable ART-free remission via immune control* (functional cure). Every 2024–26
breakthrough — the wild-type/heterozygous-CCR5 transplant remissions, the RIO
bNAb trial, the UCSF combination immunotherapy — works through the **same
mechanism: CD8⁺ T-cell-mediated control ("vaccinal effect")**, not the
intervention everyone designed for.

For an outside modeler with no wet lab, the citation graph is unambiguous: the
**predict-control** problem is owned by large, sample-rich clinical consortia and
is not enterable as primary work. The **one open, modeling-first niche** is the
small **stochastic-latency / therapeutic-interfering-particle (TIP) school**
(Weinberger–Singh–Dar–Razooky), which works in theory/simulation, is actively
publishing unresolved questions, and just crossed the NHP threshold.

---

## What's settled (don't re-litigate)

1. ART suppresses but never clears; rebound ~2–4 weeks off ART from the latent
   reservoir (resting memory CD4⁺ T cells + tissue compartments).
2. >90% of proviruses are defective; the rare **intact** ones drive rebound, and
   the reservoir is maintained by **clonal proliferation**, not ongoing replication.
3. The ~6–7 transplant *cures* / ~10 *remissions* (Martínez-Picado 2025, "power of
   ten") are non-scalable existence proofs done for blood cancers.

## What genuinely changed in 2024–26

- **CCR5Δ32 dogma broken.** Geneva patient (wild-type donor, Nat Med 2024) + 2nd
  Berlin/B2 (heterozygous donor, Nature 2025) → remission via **graft-versus-
  reservoir immunity + ADCC**, not coreceptor loss.
- **Immunotherapy delivered the first non-transplant durable remissions.** RIO
  (Lancet HIV 2026): LS-bNAbs, one participant 4 yr off ART. UCSF (Nature 2025):
  vaccine+bNAb+TLR9 → 7/10 post-intervention control. Both = vaccinal effect;
  **reservoir size did not predict control** (RIO).
- **Shock-and-kill is effectively dead** (RIVER null; vesatolimod 0.2-wk rebound
  delay; doesn't even induce expression in macaques).
- **EBT-101** CRISPR excision: Phase 1/2a complete, n=6, all 3 ATI → rebound.
- **Measurement is contested**: IPDA fails ~28% (within-subtype-B polymorphism / probe mismatch), underestimates;
  Q4ddPCR / 5-region Rainbow ddPCR emerging.
- **Block-and-lock (dCA)** still preclinical 8 yr on ("Tat inhibitors not clinically
  available"). **Anti-HIV CAR-T** preclinical, exhaustion-limited.

---

## The structural finding: the field splits into two communities

The citation graph around the anchors resolves cleanly into two non-overlapping
crowds. This is the map that decides where an outsider can stand.

### Community A — clinical viral-dynamics & correlates-of-control (sample-rich, crowded)
Perelson/Ke (Los Alamos) · Schiffer/Reeves/Cardozo-Ojeda (Fred Hutch) ·
Deeks/Peluso (UCSF) · Siliciano (JHU) · Ndung'u (AHRI) · Jonathan Li (Harvard).
Hub paper: **Gunst 2025 Nat Commun** — IPD meta-analysis of 24 ATI studies, **48
citations within a year**, all the predict-control work orbits it. Owns: reservoir
decay models, ATI biomarkers, post-treatment-control correlates. **Entry barrier:
needs primary patient samples.** An outsider can only contribute meta-analysis or
assay-harmonization here, not new correlates.

### Community B — stochastic latency & TIPs (theory-first, small, OPEN)
Weinberger (Gladstone) · Abhyudai Singh (Delaware, *control theory*) · Dar (UIC) ·
Razooky · Simpson · Soltani/Hansen/Schaffer/Sauro (systems biology). Roots in the
bistability canon (Arkin 1998 phage-λ, Ozbudak 2004 lac, Alon design principles).
This community **does not need a wet lab** — it lives in single-cell models, control
theory of the Tat positive-feedback circuit, and population-dynamics design of TIPs.
It is small (the same ~10 names recur) and is **actively publishing open problems**.

> The author-frequency ranking makes it concrete: across all papers building on the
> anchors, the top names are **Weinberger (31), Dar (22), Singh (21), Simpson (13),
> Razooky (11)** — a dynamical-systems school, not a clinical one.

---

## Where a no-lab modeler fits — ranked, with named open questions

### 1. TIP within-host dynamics & resistance-proofing  ★ best fit
Therapeutic interfering particles (engineered defective HIV that parasitizes WT
virus, R0>1, single-dose, self-spreading) crossed NHP proof-of-concept in
**Pitchai 2024 Science**. The follow-on is *almost entirely modeling*, and the
central question is **unresolved and contested**:

- **Does adaptive immunity neutralize the TIP advantage in vivo?** — **Dodd & de
  Boer 2026 (J Theor Biol, vol 619; online 2025)** argue it may. Rob de Boer is a top theoretical
  immunologist; an open within-host-dynamics debate with a major lab on record is
  exactly an enterable problem.
- **Single-cell-resolution TIP inhibition efficacy** — Sazonov/Grebennikov 2025
  (Viruses), a computational-virology group; room for better multiscale models.
- **Trial-design simulation** for TIP/cure interventions — Schiffer/Reeves 2026.

Why it fits: pure theory/simulation, small community, live disagreement, and it
sits on the Perelson viral-dynamics foundation — tractable for a dynamical-systems
person. Realistic output: a within-host model resolving the TIP-vs-immunity
question, or a design-space map for resistance-proof TIP parameters.

### 2. Control-theoretic / personalized modeling of the Tat latency switch
The Tat positive-feedback circuit is a noisy bistable switch (Weinberger 2006/07,
Burnett 2009). The field cooled (only 3 recent citations of the 2007 anchor) — but
one is **Rasi, Emili, Conway 2025 (npj Systems Biology & Applications),
"Mathematical modeling … for personalized anti-latency therapies" (8 cites)** —
reopening it as a *design* problem: use the circuit model to choose per-patient
latency interventions. Abhyudai Singh's presence (control-systems engineer, 21
appearances) confirms the door is open to control-theory framing — "collapse the
bistability so the population falls into one basin." Lower momentum than TIPs but
the most intellectually distinctive lane.

### 3. Cross-study reservoir-assay harmonization (unglamorous, real)
IPDA/Q4ddPCR/Rainbow/NFL-seq give non-comparable numbers (White 2022, Reeves 2023
"misclassified defective proviruses"). Nobody owns the normalization layer that
would let the scattered cohorts be pooled. No primary samples needed; complements
Community A without competing for samples.

### NOT recommended
- **Predict-control as primary work** — owned by Community A (Gunst hub, Ma 2026
  Immunity, Fisher 2026 Nat Immunol). Enter only via #3.
- **Any wet-lab-gated modality** — gene therapy, CAR-T, transplant, bNAb
  engineering. No compute leverage without a bench.
- **A dbvision-style federation engine** — the cure-relevant data (reservoir
  dynamics, integration-immune-escape) lives in papers/raw repos, not federatable
  DBs; the gap is data that doesn't exist cleanly, not siloed data.

---

## Concrete first project options

- **P1 (TIP-vs-immunity).** Reproduce the Pitchai/Chaturvedi TIP model, add an
  adaptive-immune compartment, test whether Dodd & de Boer's pessimism holds across
  parameter space. Deliverable: a within-host model + the conditions under which
  TIPs retain R0>1 against immunity. Smallest path to a real result.
- **P2 (latency-switch control).** Take the Tat-circuit model (Weinberger/Razooky),
  pose it as a control problem, search for noise-modulating inputs that collapse
  bistability toward a single fate. Deliverable: a control-theoretic latency map.
- **P3 (assay harmonization).** Build the cross-assay normalization model over the
  reservoir-measurement literature already in [bibliography.md](bibliography.md).

Recommendation: **P1.** It is the live, contested, modeling-first question in the
one community an outsider can join, and it needs nothing but a laptop and the
viral-dynamics canon (Perelson 1996, Finzi 1997, Chomont 2009 — all in
[citation_graph.md](citation_graph.md)).
