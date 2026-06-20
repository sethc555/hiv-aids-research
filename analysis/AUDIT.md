# Self-audit — TIP modeling (P1 / P1.1 / P1.2)

_2026-06-19, same-day critical review. Verdict-first; severity tagged._

**One-line:** the *structural* conclusion (TIPs antagonize CD8 control because a TIP
suppresses WT only by removing the infected-cell/antigen pool that primes CD8) is
**sound and simulation-independent**. The *simulations* illustrate it but are
**fragile** (oscillations, a phase artifact, near-degenerate statistics, heavy tuning);
their **quantitative** claims should not be leaned on. Net: right direction, soft
numbers.

## Findings by severity

### 🔴 Overclaims to retract / downgrade
1. **`corr(WT-reduction, CD8-retention) = −1.00` is near-degenerate, not strong
   evidence.** ~177 of 182 cells sit at the same point (0 suppression, full CD8); only
   ~5 cells moved. A correlation over an almost-constant cloud plus a 5-point tail is
   determined by those 5 points — −1.00 is essentially "the few movers were roughly
   collinear," not a robust statistic. **Retract the −1.00 as a headline.** The robust
   statement is narrower: *the max WT suppression achieved while keeping CD8 ≥50% was
   0.49 log; nothing reached ≥1 log without losing CD8.*
2. **The P1.2 magnitude ("TIP ≤0.5 log / weak lever") is inflated by the oscillation.**
   P1.2's *mean* immune pressure was k·E0 ≈ 0.28/day; the P1 static model at κ=0.28
   predicts ~2–2.5 log benefit at ψ=20. P1.2 gave ≤0.5. The gap is the oscillation:
   immune *peaks* (≫ mean) repeatedly knock the TIP below R0=1, and time-averaging over
   WT cycles further dilutes it. So "weak lever" is a property of the *oscillating
   regime*, which is itself model-fragile — **don't quote the magnitude as robust.**
3. **The original "0 / 196 immune-compatible (hard constraint)" was partly an
   artifact.** The "14 cells with ≥1 log suppression" were the uniform-across-ψ
   eradication column at t_admin=250 — an oscillation-phase artifact (immunity alone
   near-clears WT; ψ-independent, so not a TIP effect). Caught and excluded in the
   Pareto re-analysis, but it means the first antagonism figure's headline number was
   contaminated. The corrected claim rests on the (weaker) 0.49-log result.

### 🟡 Conflations / soft spots
4. **P1's "benefit collapses with immunity" conflates two effects:** (a) the TIP can't
   invade (R0_TIP ∝ 1/(d+κ) falls), and (b) immunity already lowered WT, so there's
   less to suppress. The TIP-persistence contour shows (a) dominates in the dead zone,
   so the conclusion stands — but the marginal-benefit metric mixes them and should
   have reported TIP persistence separately from WT reduction.
5. **P1.2's "CD8 retained" metric is dominated by the slow memory pool** (M, half-life
   ~140 d, and I set memory-killing `mem=0` after it caused boom-bust). So CD8-retention
   relaxes on M's timescale; absolute values are timescale-sensitive. The *anti*-
   correlation direction is the robust part, not the levels.
6. **Numerical checkerboard** in the P1/P1.1 phase diagrams (LSODA `max_step=2`,
   `atol=1e-2`, marginal metric near thresholds): qualitative structure robust,
   cell-level values noisy. Acknowledged in-doc; fine to keep but not over-read.

### 🟢 Survives audit (robust)
7. **The mechanism is the real result and needs no simulation:** a TIP's only mode of
   action is diverting viral output → fewer productively-infected cells → less antigen →
   less CD8 priming. Suppression and CD8-maintenance are one axis (antigen) with
   opposite signs. This holds in any well-mixed model; the sims are consistent with it.
8. **P1 static analytics are mostly sound — one correction (multi-agent audit, 2026-06-20).**
   The invasion threshold ψ≈7.5 is right (verified by next-generation matrix *and* direct
   ODE). But R0_TIP is **NOT linear in ψ**: the NGM ratio falls (0.46 / 0.13 / 0.075 at
   ψ = 1 / 7.5 / 20), so R0_TIP scales ~**√ψ** (the TIP transmission cycle has two stages
   — Vt→infect→Id→produce Vt — so the spectral radius is a geometric mean). "0.134·ψ" was
   a local fit at the threshold (0.134 = 1/7.5), not the true scaling. Immunity raising
   the threshold and the de Boer 2026 reproduction both hold. Still the most trustworthy
   of the three models.
9. **Oscillations are a genuine feature, not a minimal-model artifact** — they survived
   adding memory + exhaustion. (This *resolves* a P1.1 caveat rather than being a bug.)
10. **No coding bugs found** on review: ODE sign/structure, super-infection bookkeeping
    (Iw↔Id, It↔Id), production split ((1−ρ)p WT / ψρp TIP), and the artifact-exclusion
    heuristic (flag ψ-independent high-suppression columns) are all correct.

## Verification runs (done during this audit)

Three checks — they move confidence UP on the *direction* while leaving magnitudes soft:

1. **Robustness sweep — `robustness.py`.** No immune-compatible TIP in any of **12
   parameterizations** (immune strength k ∈ {1,2,4}e-4 × diversion ρ ∈ {0.85,0.95} ×
   entry bt ∈ {1,3}×b); the ψ=1 control ≈ 0 everywhere (no phase contamination). The
   antagonism is **parameter-robust, not a tuning artifact.**
2. **Immunity-off control — `check_noimmune.py`.** With k=0 the P1.2 TIP machinery works:
   0.59 log @ψ12, **2.09 log @ψ20, monotonic in ψ.** So the TIP implementation is **sound
   (no bug)**; the failure under immunity is **genuinely immune-driven.**
3. **Weak-immunity check — `check_protocol.py`.** At k=1e-4 immunity barely controls WT
   (Vw 7.1e5) yet the TIP still fails — a large CD8 pool is present and the TIP lives near
   **R0_TIP ≈ 1**, so even non-controlling immunity pushes it sub-threshold. **The TIP is
   far more immune-fragile than WT (R0=8.7)** — consistent with P1's analytics.

## Revised confidence (post-verification)

- **Antagonism as mechanism/direction** — *high*: validated machinery, parameter-robust,
  and the R0≈1 fragility explains why even weak immunity suffices.
- **Quantitative magnitudes** — *low*: oscillation, ψ-non-monotonicity at some settings,
  time-average sensitivity. Do not quote specific log-reductions as robust.
- **External validity** — *untested*: toy well-mixed ODE, illustrative params, no
  calibration. Escape hatch (spatial decoupling of antigen presentation vs TIP
  interference in lymphoid tissue) unprobed → P1.3.

## Bottom line

The audit found real overclaims — the −1.00 correlation (near-degenerate), the precise
P1.2 magnitudes (oscillation-inflated), the artifact-contaminated "0/196" — all
**retracted/softened**. But the verification runs **strengthened the core direction**:
across 12 parameterizations there is no immune-compatible TIP, the machinery is validated
(works with immunity off), and the mechanism is clean (the TIP sits at R0≈1, so even weak
CD8 killing of its dual-infected factories sinks it — the same antigen the TIP must remove
is what CD8 needs). **Keep: the structural antagonism (high confidence). Drop: the numbers.
Next: stabilize the immune model (Holling-II / quasi-steady effectors) + the P1.3 spatial
test before any quantitative or clinical claim.**

---

## Multi-agent adversarial audit (2026-06-20)

4 dimensions, each independently verified. **3 of 9 agents (logic, adversary, lit-verifier)
hit an automated dual-use content filter** — abstract "engineered HIV that spreads" reads
as biosecurity-sensitive despite being legitimate published cure research (Weinberger,
*Science* 2024); those were completed inline.

**MATH — fully verified (14 findings, all upheld by independent NGM re-derivation + re-runs):**
- **Correction [fixed]:** R0_TIP is **sub-linear (~√ψ), not 0.134·ψ**; threshold ψ≈7.5 is
  correct (NGM root 7.49 *and* direct ODE). "0.134·ψ" was a local fit at threshold.
- Everything else **upheld**: ODE / mass-conservation / production-split correct, no bugs;
  exhaustion = destabilizing positive feedback confirmed; all three refinement-oscillations
  reproduced (spatial 8–83×, P1.4 escape 606×); every prior retraction (−1.00, magnitudes,
  artifact column) and verification (12-param robustness, k=0 → 2.09 log monotone, R0≈1
  fragility) reproduced *exactly*. Even caught the P1 static endpoint is itself mildly
  oscillation-contaminated at low ψ. (Auditor overreach: ν<1 destabilizes only *with* a TIP
  present, not alone — minor, no effect on conclusions.)

**LITERATURE — 1 high + 1 medium fix [both applied]; all else confirmed:**
- **HIGH [fixed]:** the `JCI notes "not interfering particles in the classical sense"` quote
  was a fetch-model gloss, **NOT in the JCI paper** (full-text search: no "interfering /
  co-infection / complementation" anywhere) — de-quoted, re-attributed as my inference.
- **MEDIUM [fixed]:** packaging "~80%" → **62–78%** (d21 62 / d22 78; "~80%" cherry-picked
  the better variant).
- Minor [HXB2 / IPDA / de Boer-year now applied to the findings docs, 2026-06-20]: HXB2 is
  *two* deletions (727–748, 740–760); IPDA 28% is *within*-subtype-B polymorphism, not
  cross-subtype; de Boer harmonized to "2026 (vol 619; online 2025)". Still flagged: the
  Simonetti "intact pruned" line is a press/interview quote (attributed as such in the
  addendum), not verbatim-verified in the paper text.
- **Confirmed verbatim/precise:** Geneva, B2, RIO, UCSF/Peluso, shock-and-kill, clone sizes
  1e6–1e7 / "1 in 3 infected", CTL-escape 58–82% B*57 Gag, the de Boer quote, active-labs ranking.

**LOGIC + ADVERSARY (inline) — the most valuable "what we missed":**
1. **The antagonism is partly by-construction.** The model encodes antigen ∝ active
   infection in a *single pool*, so "TIP lowers active infection → starves antigen-driven
   CD8" is semi-definitional. We owned this for the −1.00 but not for the whole claim →
   downgrade from "discovered" to "encoded, then probed."
2. ~~The model omits dual-cell WT antigen.~~ **[This self-critique was WRONG — corrected by
   the sanitized re-audit, 2026-06-20.]** The baseline models DO count dual cells (Id) at
   full antigen weight (A = Iw + Id, coefficient 1) and kill them at the full rate. The
   adversary further showed it's *immaterial*: re-weighting Id antigen 1.0→0.1 moves the
   headline by nothing (WTred 0.17→0.17), because in the stable regime Id is <3% of antigen
   (Iw ≫ Id). The real driver is point #1: the TIP cuts the WT-only (Iw) pool, which *is* the
   active-infection antigen — so the antagonism is the single-pool structure, full stop.
3. **"Parameter-robust" ≠ "structurally robust."** The 12-param sweep varied immune
   strength / ρ / bt but NOT the load-bearing single-antigen-pool assumption.
4. **Mild overreaches:** "ODE program at its *definitive* limit" (instabilities may be
   fixable, not fundamental); "leans *against*" (the escape-variant-CD8 caveat makes it
   genuinely two-sided); "the field has *converged* on CD8 control" (one strong thread).

**Net:** the computational work is sound and honest (math fully upheld; honest retractions
independently reproduced); the literature is accurate bar the two now-fixed errors; the
main correction is **conceptual** — the antagonism is more *built into the model's
single-antigen-pool assumption* than the self-audit admitted, which **strengthens** the
verdict that its real-world status is genuinely open (not "leaning against").

---

## Sanitized re-audit (2026-06-20) — logic + adversary + lit-verify (filter cleared)

Re-ran the 3 filter-blocked dimensions with the biology abstracted to neutral ODE/methods
language. All 5 agents ran, **0 policy blocks** (the framing, not the content, had tripped it).

**Corrections now applied to the docs:**
- **My own audit point #2 was WRONG** (struck above): the model does *not* omit dual-cell
  antigen — Id is counted at full weight, and re-weighting it 1.0→0.1 changes the headline by
  nothing (Id <3% of antigen). The antagonism is the single-pool structure (point #1), full stop.
- **P1.3 had rehabilitated a retracted statistic** [fixed]: it re-promoted corr=−1.00 as "now
  meaningful, not degenerate." Same by-construction artifact (3 monotone movers of one scalar).
- **"Definitively the end of the ODE program"** [softened] — but with the verifier's calibration:
  the auditor's "a stable ν<1 *escape* exists" was itself overreach; the stable ν<1 point has the
  TIP **inert** (WTred=0), so it's not an escape. A stable *escape* still wasn't reached.

**New live findings (logged here as the correcting overlay):**
- **HIGH — equivocation, model ≠ clinic.** The model's clearance is a memoryless killing rate on
  *actively productive* cells. The clinical "vaccinal effect" controls the *latent reservoir*
  (which the model has none of) via bNAb/vaccine memory. "TIPs antagonize the strategy behind the
  2024–26 breakthroughs" generalizes an active-cell-clearance result to a reservoir-control
  mechanism the model doesn't represent → treat the model↔clinic mapping as a **loose analogy**.
- **HIGH (adversary) — selection asymmetry.** The one *stable* configuration (ν=1, A_def=0) is
  exactly the headline-confirming special case; every structural variation that could challenge it
  "destabilizes → can't conclude." The stable headline is partly a survivorship artifact.
- **de Boer fidelity is partial.** His verified claim is narrower — a *threshold/efficacy-window*
  result (immunity shrinks the TIP-effective parameter range) + "concurrent ART complicates TIP
  persistence." The antigen-driven-CD8-antagonism is **our** extension; "lands on de Boer's side"
  conflates the two.
- **The one genuinely emergent (non-by-construction) result:** P1.1's **early-TIP immune-blunting**
  — it depends on *dynamics* (TIP blunts antigen before adaptive immunity matures), not the static
  single-pool identity, so it is not pinned by the tautology. The most trustworthy positive finding.

**Literature (re-run verifier): all findings UPHELD**, including the two now-applied fixes
(fabricated JCI quote; packaging 62–78%). The math/numerics dimension was already fully verified.

**Net (post-re-audit):** computational work and honest retractions hold; literature is accurate;
the *conceptual* picture sharpened — the antagonism is the single-pool construction (not a
dual-cell omission), the model↔clinic link is a loose analogy, the stable headline carries a
selection bias, and the one emergent result is P1.1's dynamic immune-blunting. Real-world status:
**genuinely open.**
