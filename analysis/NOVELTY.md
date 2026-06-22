# Novelty & prior art — how this work relates to the TIP / DIP literature

_2026-06-21. A literature/novelty assessment (Semantic Scholar + web), the honest reconciliation
with the closest published model (Dodd & de Boer 2025), and a bridge experiment that recreates
their result and pinpoints **why** ours differs. Read with the project's status caveat: these are
**illustrative model results**, so "novel" means *novel as a computational hypothesis*, not an
established finding — and [AUDIT2.md](AUDIT2.md) shows how assumption-sensitive the verdict is._

## Novelty verdicts (per the indexed literature, through June 2026)

| # | Finding | Verdict |
|---|---|---|
| 1 | A TIP modeled against a **latent replication-competent reservoir** in an **ART→ATI / post-treatment-control** cure setting | **Apparently novel** — TIP models and reservoir/ATI-rebound models exist as *separate* literatures; no paper found joining them |
| 2 | A **reservoir-coupled** TIP (co-reactivates with rebounding WT) that **helps** CD8 / immune-mediated post-treatment control | **Apparently novel as a TIP claim** — but see the de Boer reconciliation below; this is an *extension*, not a refutation |
| 3 | A TIP **sparing the immune response from exhaustion** by lowering antigen burden (benefit largest when immunity is exhaustible) | **Principle is known, application is novel** — antigen-burden→exhaustion is textbook (Mueller & Ahmed 2009); defective-virus "spares immunopathology" is known in flu/SARS-CoV-2 (Dimmock; DVG work); applying it to a TIP, to HIV, and showing exhaustion-*dependence* is not in the indexed literature |

## The closest paper, read in full, and the honest reconciliation

**Dodd & de Boer, "Immune responses may make HIV-1 therapeutic interfering particles less
effective," _J Theor Biol_ 2025** (DOI 10.1016/j.jtbi.2025.112317; PubMed 41319726). Read via the
[Utrecht portal](https://research-portal.uu.nl/en/publications/immune-responses-may-make-hiv-1-therapeutic-interfering-particles/)
+ Semantic Scholar abstract (no open-access full text). Their model is an ODE of **active
infection only** — **no latent reservoir, no ATI, no rebound/post-treatment control**. Their
headline: *"even a moderate immune response against virally infected cells drastically decreases
the range of parameter values for which therapy is effective"* — i.e. **immunity hurts the TIP's
suppression of active virus.**

**We do NOT flip this.** (An earlier auto-summary called it an "inversion"; that overclaims.) It
is a **different question in a different setting**, and where the models overlap we *agree*:

| | Dodd & de Boer 2025 | This project |
|---|---|---|
| Question | does immunity help/hurt **the TIP**? | does a TIP help/hurt **the immune cure**? |
| Setting | active infection, no reservoir | latent reservoir + ART→ATI (durable remission) |
| Answer | immunity *hurts* the TIP | a *coupled* TIP *helps* the cure |

## The bridge experiment — recreate their result, then show *why* ours differs

[tip_model_deboer_bridge.py](tip_model_deboer_bridge.py) runs **both** metrics on the **same**
within-host TIP machinery (our [tip_model.py](tip_model.py) is a faithful structural recreation of
their active-infection ODE):

- **Metric A — de Boer's question** (marginal TIP suppression of active viral load vs immune
  strength `kap`): **1.44 → 0.10 → 0.00 → 0.00 log** as immunity rises. **We reproduce de Boer:
  immunity collapses the TIP's active-suppression benefit.**
- **Metric B — our question** (coupled TIP's effect on *durable post-treatment control* vs immune
  strength, in the reservoir/ATI setting): TIP effect **+3 / +33 / +10 points** across rising
  immunity. **A coupled TIP raises durable control: it helps the cure.**

**Same TIP, same "immunity kills infected cells" mechanism — opposite sign.** Why:
1. **The objective flips:** "suppress active viral load *now*" (they compete with immunity over the
   same active WT → immunity shrinks the TIP's niche) vs "prevent the reservoir rebounding *later*"
   (the TIP and immunity *cooperate* against the reservoir).
2. **The ingredients that enable cooperation are absent from their model:** the latent reservoir,
   the ATI, and the TIP coupling to reservoir reactivation. Add them and competition becomes
   cooperation; remove them and we recover de Boer.
3. Their own conclusion is "TIP success depends *even more strongly on the immune response*… hard
   to predict." Our contribution **refines that "it depends"** into *how*: **decoupled → neutral/hurt
   (their regime); coupled to the reservoir → helps (the cure regime).**

## Bottom line

The genuinely new pieces are the **reservoir/ATI cure setting (#1)** and the **coupled-TIP-helps /
exhaustion-sparing mechanism (#2, #3)** — framed honestly as **extending and refining** Dodd & de
Boer (consistent with them where the models overlap), not overturning them, and applying a known
exhaustion principle (Mueller & Ahmed) to a new context. All of it is a **modeling hypothesis**
that the audit shows is assumption-sensitive — the next step to claim priority for real is a
full-text read of de Boer's parameters, then wet-lab/clinical testing, not more simulation.

## Key references
- Dodd & de Boer 2025, _J Theor Biol_ — [PubMed 41319726](https://pubmed.ncbi.nlm.nih.gov/41319726/)
- Pitchai/Tanner et al. 2024, _Science_ 385:eadn5866 (TIP reduces SHIV viremia in NHPs) — [PMC11545966](https://pmc.ncbi.nlm.nih.gov/articles/PMC11545966/)
- Sazonov/Grebennikov/Bocharov 2025, _Viruses_ (single-cell TIP efficacy) — [PMC12567885](https://pmc.ncbi.nlm.nih.gov/articles/PMC12567885/)
- Karki, Bull & Krone 2022, _Virus Evolution_ (DIP + immunity modeling) — DOI 10.1093/ve/veac033
- Conway & Perelson 2015, _PNAS_ (post-treatment control model, no TIP) — [PMC4418889](https://pmc.ncbi.nlm.nih.gov/articles/PMC4418889/)
- Mueller & Ahmed 2009, _PNAS_ (antigen burden → CD8 exhaustion) — DOI 10.1073/pnas.0809818106
- Dimmock et al., DI-influenza immune modulation — [PMC3163266](https://pmc.ncbi.nlm.nih.gov/articles/PMC3163266/)
- SARS-CoV-2 DVG protective immunity / reduced immunopathology, _Cell_ 2021 — DOI 10.1016/j.cell.2021.09.027
