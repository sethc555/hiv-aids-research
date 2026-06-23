# A reservoir-coupled therapeutic interfering particle can assist immune-mediated HIV-1 post-treatment control: a within-host modeling study

**Author:** Seth Cope¹  _(¹ Independent researcher. ORCID: [0009-0000-5520-915X](https://orcid.org/0009-0000-5520-915X))_

_Illustrative within-host modeling study / hypothesis generation — not validated experimental or
clinical findings, and not medical advice. Preprint: Zenodo DOI
[10.5281/zenodo.20801699](https://doi.org/10.5281/zenodo.20801699)._

---

## Abstract

**Background.** Therapeutic interfering particles (TIPs) are engineered, conditionally-replicating
defective HIV-1 genomes proposed as single-administration antivirals. Existing within-host models
treat active infection; Dodd & de Boer (2025) found that an immune response *reduces* the parameter
range over which a TIP is effective. Whether a TIP helps or hurts the immune-mediated cure — durable
ART-free remission, governed by the latent reservoir and post-treatment control (PTC) — has not been
modeled. **Methods.** We add a latent replication-competent reservoir and an ART → analytical
treatment interruption (ATI) schedule to a within-host TIP model, with stochastic (tau-leaping)
dynamics, calibrated to clinical rebound timing and PTC fractions (A5345, ACTG pooled, CHAMP, RIO).
We introduce one coupling parameter, χ, the fraction of reservoir reactivations that co-introduce the
TIP, and test across three structurally distinct immune models. **Results.** A non-coupled TIP (χ=0)
is neutral to PTC and recovers the de Boer active-infection limit. A TIP coupled to reservoir
reactivation raises durable control monotonically in χ, helps most for marginal controllers, and
shows no systematic backfire — robust across global parameter sampling and all three immune models, and
protective even under active immune exhaustion. A derived effective reproduction number,
R_eff = R₀·d/(d+κ), explains the effect: coupling lowers the immune threshold for control.
**Conclusions.** We predict that a reservoir-co-residing TIP could assist immune-mediated PTC,
testable in an ATI animal model. This is an illustrative hypothesis conditional on the coupling
assumption.

---

## 1. Introduction

Antiretroviral therapy (ART) suppresses HIV-1 but does not eliminate the latent reservoir, so
treatment is lifelong and the virus rebounds on interruption. A leading cure goal is durable
ART-free remission — *post-treatment control* (PTC), in which the reservoir persists but rebound is
contained after ART interruption. Its mechanism is multifactorial and not settled: PTC has been associated
with a smaller reservoir, cytotoxic CD8 responses, NK-cell activity, and humoral immunity, with the dominant
correlate differing between individuals and settings (Mesquita & Li 2024; Blazkova et al. 2021). CD8 control
is one contributor — clearest under early ART and antibody/combination immunotherapy (Passaes et al. 2024) —
but the canonical spontaneous-PTC cohort (VISCONTI) controls with *weak* CD8 responses, and
broadly-neutralizing-antibody trials such as RIO act substantially through antibody and reservoir effects
(the "vaccinal" effect; Tipoe & Fidler 2022). We therefore model the immune axis generically, as
cytotoxic-effector killing of antigen-expressing cells (CD8 and/or NK/ADCC), not a CD8-specific pathway.

A therapeutic interfering particle (TIP) is an engineered defective HIV-1 genome that replicates only
in cells co-infected by wild-type (WT) virus, diverting WT packaging to itself; a single dose reduced
SHIV viremia by >3 log₁₀ (>1000-fold) in non-human primates, the TIP conditionally replicated for
>6 months, and "TIP-treated animals exhibited significantly improved immune responses, with no evidence
of increased inflammation" (Pitchai et al., *Science* 2024) — an empirical hint of TIP–immune interplay
that has not been modeled. *In vitro*, the same study found that upon ART cessation the TIP reactivated
together with HIV and interfered with its outgrowth — a non-clinical hint of the reservoir co-reactivation
this model formalizes as the coupling χ (the in vivo NHP work included no ART-interruption arm).

Within-host TIP theory has, to date, treated active infection only. Dodd & de Boer (2025, *J Theor
Biol*) derive the TIP basic reproduction number and show that an immune response against infected
cells "drastically decreases the range of parameter values for which therapy is effective" — i.e.
immunity hurts the TIP's suppression of active virus. Their model contains no latent reservoir and no
treatment interruption. The cure-relevant question — does a TIP help or hurt *immune-mediated
post-treatment control*? — therefore remains open. We address it.

## 2. Model and methods

We use a Perelson-class within-host model (target cells T; WT-only, TIP-only, and dually-infected
productive cells; quasi-steady free virus), with WT basic reproduction number R₀ = b·T₀·p/(d·c) =
8.70. To this we add (i) a latent replication-competent reservoir L_lat that reactivates to a
productive cell, (ii) a defective antigen-presenting clone (Simonetti-type) that primes CD8
independently of active WT, and (iii) an ART → ATI schedule (chronic infection → ART suppression →
interruption). Dynamics are stochastic (tau-leaping); durable post-treatment control is defined as
active-infection extinction with the reservoir persisting (functional, not sterilizing).

The central new parameter is the **coupling** χ ∈ [0,1]: the fraction of reservoir reactivations that
emerge as TIP-carrying (dual) cells rather than pure WT cells — i.e. how well an engineered TIP
co-resides in / co-reactivates with the reservoir. χ=0 reproduces a non-interacting (de Boer-limit)
TIP; χ=1 is a fully reservoir-coupled TIP.

We test robustness across three structurally distinct immune models: quasi-steady killing with a
maintained antigen floor; a dynamic waning effector memory (no floor); and an actively-exhaustible
memory (burden-driven degradation). In all three, the killed compartment is antigen-expressing productive
cells, so the killing rate κ denotes cytotoxic-effector pressure generically — CD8 and/or NK/ADCC — which
keeps the result independent of the unresolved question of which effector dominates post-treatment control
(a deliberate generalization; κ is antigen-driven, fitting adaptive CD8 best). The rebound clock and PTC fractions are calibrated to clinical
data (untreated/placebo ATI median ~16–22 d; ACTG rebound-by-week-4/12; CHAMP spontaneous PTC ~4–13%;
RIO bNAb durable control ~24% — the 7/29 ATI sub-analysis, vs ~21% (7/34) by the trial's primary
endpoint, HR 0.09). Parameters and equations are documented in `analysis/METHODS.md`; all
code is open and a single script (`analysis/verify_claims.py`) re-derives every headline number,
including a reproduction of the Dodd & de Boer result.

## 3. Results

**3.1 A non-coupled TIP is neutral; coupling flips it to helpful.** With the reservoir modeled but
the TIP decoupled (χ=0), a TIP does not change post-treatment control — recovering the de Boer
"immunity decides" outcome. Coupling the TIP to reservoir reactivation makes control rise
**monotonically** with χ: at a marginal immune level, durable control climbed 22% → 32% → 40% → 45% →
58% as χ went 0 → 0.25 → 0.5 → 0.75 → 1, with CD8 maintained throughout (the coupled TIP intercepts
the rebound without starving immunity).

**3.2 The benefit is robust across immune structure, with no systematic backfire.** The result holds
across all three immune models. Under a wear-down-able memory, a coupled TIP still helped (7/9 regimes)
and did not harm control. Under active, burden-driven immune exhaustion combined with a stealthy
(low-visibility) TIP — the configuration most likely to backfire — there was still no backfire: the TIP
helped *most* there (e.g. 17% → ~65% control), because by suppressing the rebound it keeps the antigen
burden low and **spares the immune response from exhaustion**. (Across the broader phase scans, occasional
small negative nudges of ≤1 point appear within stochastic noise — see `analysis/AUDIT2.md` — but no
configuration showed a systematic backfire.)

**3.3 Phase boundary and reduction to prior work (Fig. 1).** As a continuous function of χ and immune
strength, the TIP effect reduces to the de Boer limit at χ=0 and grows with χ; the benefit is
concentrated at *intermediate / marginal* immune strength (where control is on a knife-edge) and
requires substantial coupling (χ ≳ 0.6).

**3.4 Global sensitivity.** Across 30 Monte-Carlo draws varying all six key parameters jointly, the
coupled-TIP effect helped (63%) or was neutral (37%) and **harmed in 0%** — no backfire across the
sampled space. The dominant driver is immune strength (the TIP helps where immunity is marginal); the
reservoir's reactivation timing barely matters.

**3.5 A derived criterion (Fig. 2).** Post-ART, a reactivating WT lineage has effective reproduction
number **R_eff(WT) = R₀·d/(d+κ)** (κ = immune killing rate); control corresponds to sub-criticality,
R_eff < 1, i.e. κ > κ_crit = (R₀−1)d = 7.70/day. A coupled TIP shifts the effective control threshold
*down* by Δ(χ) ≥ 0. This single inequality explains all three empirical findings: the benefit is
concentrated just below threshold (marginal immunity), widens with coupling, and — since Δ ≥ 0 —
never raises the threshold (no backfire). In a static-κ model where R_eff is exact, the TIP effect is
large in the sub-/near-threshold band and ceilings out once R_eff < 1, exactly as predicted.

## 4. Discussion

A TIP modeled against the latent reservoir in a post-treatment-control setting is, to our knowledge,
new; the existing TIP literature treats active infection, and reservoir/ATI models contain no TIP.
Our result does **not** contradict Dodd & de Boer: at zero coupling we reproduce their finding, and
the help appears only in a regime (reservoir + ATI) and via a mechanism (TIP–reactivation coupling)
their model does not contain. The mechanism is intuitive: on the active-infection axis a TIP and
immunity *compete* over the same wild-type, but on the reservoir-control axis a coupled TIP and
immunity *cooperate* — the TIP caps each reactivation burst while immunity clears it, and by keeping
the burden low it protects immunity from exhaustion.

**Falsifiable prediction.** A TIP engineered to co-reside in / co-reactivate with the latent reservoir
(high χ) should improve durable post-treatment control after ATI, while a non-reservoir-coupled TIP
(χ≈0) should be neutral; the benefit should be largest in intermediate/marginal controllers and
require substantial coupling. Direct test: a humanized-mouse or NHP ATI study comparing (a) no TIP,
(b) a standard TIP, (c) a reservoir-targeting TIP, measuring time-to-rebound and control fraction; the
model predicts (c) > (b) ≈ (a). A null result in arm (c) would falsify the coupling mechanism.

**Limitations.** This is an illustrative within-host model: parameters are illustrative beyond
calibrated rebound timing; "control" is functional (active-infection extinction with the reservoir
persisting), not sterilizing; coupling χ is a coarse single knob, not a mechanistic co-packaging
model; Δ(χ) is shown ≥0 numerically, not in closed form (a full next-generation-matrix treatment of
the nonlinear TIP interference is the natural next step). The clinical target is itself uncertain:
post-treatment control is multifactorial and its mechanism unresolved — reservoir size, CD8, NK, and humoral
immunity are all implicated and vary by individual — so the effector-killing axis here is a deliberate
generalization, not a validated CD8 mechanism; whether immunity stays primed through suppression (the antigen
"floor" that makes control achievable, modeled via defective-clone antigen presentation) is a further
load-bearing assumption; and there is currently no evidence that a therapeutically delivered TIP localizes to
the pre-existing replication-competent reservoir, so χ>0 is an engineering goal, not an established property
(see `analysis/AUDIT3.md`). The conclusions are conditional on the
coupling assumption. A multi-agent adversarial audit (in the repository) retracted several of the
author's own earlier overclaims and is the reason the framing is "conditional/assist," not "cure."

## Figures

- **Fig. 1** — Coupling phase diagram: P(durable control) and the TIP effect vs (coupling χ, immune
  strength). `analysis/p14_coupling_phase.png`.
- **Fig. 2** — Semi-analytical criterion: coupling lowers the control threshold (left), and the TIP
  effect is concentrated in the sub-threshold band where R_eff = R₀d/(d+κ) > 1 (right).
  `analysis/p16_analytic.png`.
- **Supp.** — ATI rebound / post-treatment-control calibration (`analysis/p4_ati.png`); global
  sensitivity (`analysis/p15_sensitivity.npz`).

## Data and code availability

All code, figures, the literature corpus, the verification harness, and the audit trail are openly
available: https://github.com/sethc555/hiv-aids-research (archived at Zenodo, DOI:
[10.5281/zenodo.20799761](https://doi.org/10.5281/zenodo.20799761)). Running
`cd analysis && python3 verify_claims.py` re-derives every headline number (22/22), including a
reproduction of Dodd & de Boer.

## Disclosures

**Status:** illustrative within-host modeling / hypothesis generation; not validated experimental or
clinical findings, not medical advice. **AI assistance:** this work was developed with substantial
help from an AI coding/analysis assistant (Anthropic Claude) for implementation, derivation, audit,
and drafting, under the author's direction; AI tools are not authors. **Competing interests:** none.
**Funding:** none.

## References (key — complete in `analysis/NOVELTY.md`)

1. Dodd GK, de Boer RJ. Immune responses may make HIV-1 therapeutic interfering particles less
   effective. *J Theor Biol* 2025; DOI 10.1016/j.jtbi.2025.112317.
2. Pitchai FNN, Tanner EJ, … Weinberger LS. Engineered deletions of HIV replicate conditionally to
   reduce disease in nonhuman primates. *Science* 2024; 385:eadn5866.
3. Karki B, Bull JJ, Krone SM. Modeling the therapeutic potential of defective interfering particles
   in the presence of immunity. *Virus Evolution* 2022; 8:veac047.
4. Conway JM, Perelson AS. Post-treatment control of HIV infection. *PNAS* 2015; 112:5467–5472.
5. Mueller SN, Ahmed R. High antigen levels are the cause of T cell exhaustion during chronic viral
   infection. *PNAS* 2009; 106:8623–8628.
6. Lee MJ, Fidler S, Frater J, et al. (RIO trial). *Lancet HIV* 2026.
7. Li JZ, et al. (A5345) *Clin Infect Dis* 2022; 74:865. Namazi G, et al. (CHAMP) *J Infect Dis*
   2018; 218:1954.
8. Mesquita FS, Li Y, Li JZ. Viral and immune predictors of HIV post-treatment control. *Curr Opin HIV
   AIDS* 2024; DOI 10.1097/COH.0000000000000898.
9. Passaes C, Desjardins D, et al. Early ART favors post-treatment SIV control associated with expanded
   enhanced-memory CD8+ T cells. *Nat Commun* 2024; DOI 10.1038/s41467-023-44389-3.
10. Tipoe T, Fidler S, et al. How broadly neutralizing antibodies might induce HIV remission: the
    'vaccinal' effect. *Curr Opin HIV AIDS* 2022; DOI 10.1097/COH.0000000000000731.
11. Blazkova J, Gao F, et al. Distinct mechanisms of long-term virologic control in two HIV-infected
    individuals. *Nat Med* 2021; DOI 10.1038/s41591-021-01503-6. (Plus, on CD8 quality/persistence:
    Hersperger & Migueles 2011, *Curr Opin HIV AIDS*; Blankson 2024, *Nat Immunol*; Kiani et al. 2025,
    *Nature*. See `analysis/bibliography.md` and `analysis/AUDIT3.md`.)
