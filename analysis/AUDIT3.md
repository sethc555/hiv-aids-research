# AUDIT3 — contested-biology review (2026-06-22)

_A third audit, with a different lens from [AUDIT.md](AUDIT.md) and [AUDIT2.md](AUDIT2.md). Those asked
"is the math right?" and "is the model **structure** load-bearing?" This one asks the question that forced
a rework of the sibling type-1-diabetes study: **are the underlying biological premises actually settled,
or are we treating contested biology as established fact?** Method: re-read the model's biological
assumptions ([METHODS.md](METHODS.md), [NOVELTY.md](NOVELTY.md), the model code) against the current indexed
literature (Semantic Scholar Graph API + web, June 2026); each candidate weighed for (a) is it contested,
(b) is it load-bearing, (c) is it disclosed. Supporting references are in [bibliography.md](bibliography.md)._

**One-line verdict:** good shape on this axis. The single most speculative premise — the TIP↔reservoir
coupling χ — is already labelled a conditional hypothesis throughout, which is the correct treatment of
unestablished biology. **One genuine issue:** the manuscript frames post-treatment control as *immune /
"notably CD8"*-mediated more confidently than the literature supports — PTC mechanism is multifactorial and
unresolved, and the canonical spontaneous-PTC cohort (VISCONTI) controls with **weak** CD8. This is a
**framing/citation** fix (the model's killing term is effector-agnostic), not a model rebuild — and the
supporting literature was *already* in `bibliography.md` but never propagated into the manuscript.

---

## 🔴 HIGH — F1. "CD8 / immune-mediated PTC" stated as the settled paradigm

**As written** (`paper/MANUSCRIPT.md`, Intro): *"The leading cure paradigm is durable ART-free remission via
immune (notably CD8 / 'vaccinal-effect') control of the reservoir."* The whole contribution is "a TIP
**assists immune-mediated PTC**," so this premise is load-bearing for clinical relevance.

**What the literature actually says — contested, multifactorial:**
- PTC mechanism is multifactorial and differs between individuals — reservoir size, CD8, NK, and humoral
  immunity all implicated, no single settled driver (Mesquita & Li 2024 review; Blazkova & Gao 2021,
  *"Distinct mechanisms ... in two individuals"*).
- The canonical spontaneous-PTC cohort, **VISCONTI, has *weak* HIV-specific CD8** — and NK (HLA-B35/KIR)
  and humoral arms are implicated instead (Mouquet 2024; Webb 2025).
- Reservoir **size** is a co-equal determinant (Pasternak 2021; Charre 2024); Conway & Perelson's own PTC
  model (ref 4) has a CTL-strength band where the *reservoir*, not immunity, sets the outcome.

**Honest counterweight — CD8 is *not* irrelevant:** it contributes clearly under **early ART** (Passaes
2024, SIV, enhanced-memory CD8) and **antibody/combination immunotherapy** (Kiani 2025, CD8 stemness
precedes post-intervention control; Moriarty 2025, CD8α depletion abrogates SIV control). So the defensible
statement is *"CD8 is one of several effectors — dominant in some settings, weak in spontaneous control,"*
not *"the paradigm."*

**Why it does not sink the result:** the model's `κ` is a death rate on **antigen-expressing productive
cells** — nothing in the math is CD8-specific. NK and ADCC also kill antigen-expressing infected cells, so
reading `κ` as **generic cytotoxic-effector pressure (CD8 and/or NK/ADCC)** *generalizes* the result and is
better supported than a CD8-specific claim. (Partial caveat: `κ` is antigen-*driven*, which fits adaptive
CD8 better than innate NK; antibodies act on free virus/entry — the generalization is good, not perfect.)

**Fix [applied to manuscript]:** soften "the leading paradigm … notably CD8" → multifactorial/unresolved
with the reviews cited; state `κ` as effector-agnostic in §2; add a PTC-mechanism limitation.

## 🟡 MEDIUM — F2. RIO / bNAb enlisted for "CD8 / immune control"

**As written:** RIO and combination immunotherapy cited as support for *immune (CD8)* control.
**Reality:** bNAb-induced remission is substantially **humoral / "vaccinal."** RIO's own correlates are
**autologous neutralizing antibodies + reservoir clearance** (Fumagalli 2025, already in `bibliography.md`);
the "vaccinal effect" is antibody-*induced* and may boost CD8 only downstream (so mentioning CD8 is not
wrong, just imprecise) (Tipoe & Fidler 2022; Naranjo-Gomez 2023). **Fix [applied]:** cite the vaccinal-effect
work; do not equate RIO with direct CD8 killing.

## 🟢 LOW / SELF-CORRECTION — F3. "antigen floor keeps CD8 primed" vs CD8 waning

The first draft of this audit flagged the maintained antigen floor (load-bearing per AUDIT2 #3) as
contradicting the observed decline of HIV-specific CD8 on ART. **The literature corrected me:** CD8
frequency declines after ART start but then **stabilizes**, and long-term ART can *rejuvenate* HIV-specific
CD8 (Blankson 2024); and it is CD8 **quality**, not frequency, that tracks control (Hersperger & Migueles
2011). So "CD8 stays primed" is biologically defensible — **downgraded.** The residual, genuine caveat is
narrower: the floor's *source in the model* (defective-clone antigen presentation, Simonetti-type) sustaining
priming through suppression is the unverified part — not CD8 persistence per se. Logged as a limitation, not
a retraction of the result. _(This self-correction is kept on the record deliberately — the trail is the
credibility.)_

## ✅ CONFIRMED WELL-HANDLED — F4. the χ coupling is the biggest leap, and it is labelled as such

The load-bearing novel premise — a TIP co-residing in / co-reactivating with the latent reservoir (χ>0) —
has **zero support in the indexed literature**: the TIP field is entirely active-infection (Tanner 2019
engineering; Sazonov/Grebennikov 2025 single-cell efficacy; Dodd & de Boer 2025) and nascent (single-digit
citation counts). Confirming the *absence* is the point — and the manuscript already labels χ "conditional,"
"a coarse single knob, not a mechanistic co-packaging model," with a falsifiable arm-(c) test. **This finding
is the template the other three were brought up to.** One addition for candor (manuscript Limitations): state
plainly that there is currently *no* evidence a therapeutically delivered TIP localizes to the pre-existing
replication-competent reservoir — achieving χ>0 is an unsolved engineering problem, not a tuning choice.

---

## What this does NOT change
- **No numbers move.** This audit concerns biological framing, not computation; `verify_claims.py` (22/22)
  and `attestation.json` are unaffected, and the math/numerics were already verified in AUDIT/AUDIT2.
- The **core result stands** as a conditional modeling hypothesis: *given reservoir coupling and a maintained
  effector/antigen floor, a reservoir-coupled TIP can assist effector-mediated post-treatment control.*
  Reading "immune" as effector-agnostic makes it **more** robust, not less.

## The meta-point (why this lens was worth running)
AUDIT/AUDIT2 are deep on model structure and numerics but took the **biological backdrop** (CD8-mediated PTC)
as given. The supporting-yet-complicating literature was already in `bibliography.md`; the failure was not
collecting it but **not letting it touch the framing.** Same lesson as the T1D rework — gather the
contested-biology evidence, then make sure the manuscript reflects it.

## Corrections checklist
- [x] Intro: "leading paradigm … notably CD8" → multifactorial/unresolved, cited (F1)
- [x] §2: state `κ` as effector-agnostic — CD8 and/or NK/ADCC (F1)
- [x] Limitations: PTC-mechanism uncertainty + antigen-floor-source caveat + "no evidence for χ>0 yet" (F1, F4)
- [x] References: add Mesquita 2024, Passaes 2024, Tipoe & Fidler 2022
- [x] `bibliography.md`: add the 6 refs below (dated appendix), preserving the generated bib's provenance
- [ ] _(optional)_ verbatim re-check the two Pitchai/*Science* 2024 quotes used in the Intro ("suppressed HIV
      rebound after ATI"; "improved immune responses") — flagged, not yet done
- [ ] _(public propagation, after sign-off)_ rebuild Pages (`pubkit pages`) + refresh the Zenodo preprint;
      the title/abstract use "immune-mediated" (already effector-agnostic) and need no change

## References added to bibliography.md by this review
Passaes 2024 (*Nat Commun*); Kiani 2025 (*Nature*); Blazkova & Gao 2021 (*Nat Med*); Blankson 2024 (*Nat
Immunol*); Hersperger & Migueles 2011 (*Curr Opin HIV AIDS*); Tipoe & Fidler 2022 (*Curr Opin HIV AIDS*).
Already present and reused: Mesquita 2024, Charre 2024, Mouquet 2024, Pasternak 2021, Webb 2025, Moriarty
2025, Naranjo-Gomez 2023, Fumagalli 2025 (RIO correlates).

---
_This audit changed framing, not findings._
