# Summary, citation, and disclosures

_The honest, citable one-page summary of this project — the version to put at the top of a
preprint or link in an email. Read the [README](README.md) caveat and [AUDIT2.md](analysis/AUDIT2.md)
before quoting any number._

## Title (working)

**A reservoir-coupled therapeutic interfering particle can assist immune-mediated HIV-1
post-treatment control: a within-host modeling study**

## Abstract

**Background.** Therapeutic interfering particles (TIPs) are engineered, conditionally-replicating
defective HIV-1 genomes proposed as single-administration antivirals. Existing within-host models
(notably Dodd & de Boer, *J Theor Biol* 2025) treat *active* infection and find that an immune
response *reduces* the parameter range over which a TIP is effective. Whether a TIP helps or hurts
the **immune-mediated cure** — durable ART-free remission, which is governed by the *latent
reservoir* and post-treatment control — has not, to our knowledge, been modeled.

**Methods.** We built a within-host model that adds a latent replication-competent reservoir and an
ART → analytical-treatment-interruption (ATI) schedule, with stochastic (tau-leaping) dynamics, and
calibrated its rebound timing and post-treatment-control fractions to clinical data (A5345, ACTG,
CHAMP, RIO). We introduce a single coupling parameter χ — the fraction of reservoir reactivations
that co-introduce the TIP — and test across three structurally distinct immune models
(quasi-steady with an antigen floor; waning memory; actively-exhaustible).

**Results.** A *non-coupled* TIP (χ=0) is neutral to post-treatment control and **recovers the
de Boer active-infection limit**. A TIP **coupled** to reservoir reactivation raises durable control
**monotonically in χ**, helps most for "marginal controllers," and **never backfires** — robust
across global (joint) parameter sampling and all three immune models, and protective even when
immunity is actively exhausted. A derived effective reproduction number, **R_eff = R₀·d/(d+κ)**,
explains the effect: coupling lowers the immune threshold for control by Δ(χ) ≥ 0.

**Conclusions.** We predict that a reservoir-co-residing TIP could **assist** (not antagonize)
immune-mediated post-treatment control — testable in an ATI animal model comparing a
reservoir-targeting TIP vs a standard TIP vs none. This is an **illustrative modeling hypothesis**:
parameters are illustrative beyond calibrated rebound timing, the result is conditional on the
coupling assumption, "control" is functional (not sterilizing), and a closed-form Δ(χ) remains open.

## The one falsifiable prediction

A TIP engineered to co-reside in / co-reactivate with the latent reservoir (high χ) should improve
durable post-treatment control after ATI, while a non-reservoir-coupled TIP (χ≈0) should be neutral;
the benefit should be largest in intermediate/marginal controllers and require substantial coupling
(χ≳0.6). A null result in the coupled arm would falsify the mechanism.

## Reproducibility

Everything is open and machine-checkable. `cd analysis && python3 verify_claims.py` re-derives every
headline number from the code (22/22), including a reproduction of the Dodd & de Boer result. The
project includes a same-day self-audit and a multi-agent adversarial re-audit ([AUDIT.md](analysis/AUDIT.md),
[AUDIT2.md](analysis/AUDIT2.md)) that retracted several of the author's own overclaims, plus a
prior-art / novelty assessment ([NOVELTY.md](analysis/NOVELTY.md)).

## How to cite

> Seth C. (2026). *A reservoir-coupled therapeutic interfering particle can assist immune-mediated
> HIV-1 post-treatment control: a within-host modeling study.* GitHub: sethc555/hiv-aids-research.

For a citable DOI, archive a release of the repository via [Zenodo](https://zenodo.org) (GitHub →
Zenodo integration), or post as a preprint (see below).

## Disclosures

- **Status:** illustrative within-host modeling / hypothesis generation. **Not** validated
  experimental or clinical findings, and **not** medical advice or evidence of a cure.
- **AI assistance:** this project was developed with substantial help from an AI coding/analysis
  assistant (Anthropic Claude) for implementation, derivation, auditing, and drafting, under the
  author's direction. AI tools are not authors; their use is disclosed here per current norms.
- **Competing interests:** none.
- **Data/code availability:** all code, figures, and the literature corpus are in this repository.
