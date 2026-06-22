# Re-audit (multi-agent + self-verified, 2026-06-21)

_A second adversarial audit covering everything added since [AUDIT.md](AUDIT.md): P2–P8, the
clinical calibration, the cross-disease research, and the `verify_claims.py` harness. Method:
five independent read-only auditors (math, literature, logic, code, adversary), each finding
then weighed; the load-bearing numeric claims were re-run by the maintainer under a hard memory
cap. (An earlier attempt to run the auditors as compute-capable parallel agents OOM'd the box —
hence read-only agents + maintainer-run numerics.)_

**One-line verdict:** the computational work is numerically clean and the clinical calibration
is accurate and honestly caveated — but the headline *"a TIP is orthogonal to the cure"* is
**conditional on two structural modeling choices that were under-disclosed relative to how much
they carry**, and the verification harness was testing a **duplicate** of the model, not the
production code. Downgrade "orthogonal / empirically-survived" → **"neutral within one model
family, conditional on a non-coupled reservoir and a maintained antigen floor."**

---

## 🔴 HIGH — confirmed, require correction

1. **`verify_claims.py` tested a re-implementation, not the production engine.** [code; confirmed
   by inspection] Its `_ati_control` hand-copies the ATI loop, and `cohort_burn` calls
   `simulate(..., t_ati=0)` — so `simulate()`'s **entire ATI branch is never exercised** by the
   harness. A regression in the production ATI math (the `react`/`f_lat`/`tip_sustained`
   additions) would not trip a check. The "a regression in any phase script trips a FAIL" claim
   was false for the ATI dynamics. **Fix:** route the P4/P6/P8 checks through the real
   `simulate()` / `rebound_curve()` / `run()`. _[FIXED — see verify_claims.py v2.]_

2. **The TIP↔reservoir compartment decoupling *forces* the neutrality.** [adversary #1] `L_lat`
   reactivates only to `Iw`; the engineered TIP's `L_latd` is a **separate** pool reactivating
   only to `Id` (`tip_model_p4_reservoir.py` lines ~101–103). A reactivating WT provirus therefore
   never transits a co-infected state the TIP could parasitize — which is the *sole* reason "the
   TIP has no substrate where the cure succeeds." But a TIP is *designed* to mobilize on WT
   co-infection; a coupled model (TIP co-packages / co-reactivates with rebounding WT) would let
   it engage exactly the rebound the model says it misses, plausibly reopening the P1 antagonism.
   **The model encodes non-interaction, then reports the non-interaction as a finding.** This is
   the P4 analogue of the P1 single-pool tautology AUDIT.md owned — and was not flagged with equal
   candor. **Fix:** state it as *the* load-bearing assumption; build a coupled-TIP successor.

3. **The `Ldef` antigen floor is load-bearing for control itself — verified.** [adversary #2 +
   maintainer test] `E = Estar(Iw + ν·Id + Ldef)` with `Ldef` self-maintaining at ~`RDEF`=5000,
   so CD8 stays primed **even at zero active infection**. Maintainer re-run varying `RDEF`:
   control = **0% at RDEF=500, 100% at RDEF=50000** (all kf) — i.e. whether post-treatment control
   is even *achievable* is set by the floor, not just the immune strength. The QSS structure also
   removes any antigen-driven CD8 contraction/exhaustion (dropped at P1.3), which is precisely the
   channel by which a TIP could harm the cure. **Fix:** disclose as a load-bearing sensitivity;
   the dynamic-exhaustion successor (the abandoned P9) is the right test.

4. **Cross-disease CML/cancer claims are uncited in-repo.** [literature] P8 / NARRATIVE lean on
   "CML deep-response × duration," Lenaerts 2010, EURO-SKI, Baar–Bovier, Mylvaganam 2019 — none
   appear in `bibliography.md`, `scan_results.md`, or any findings file. The sources are real, but
   as *used here* they are uncredited assertions. **Fix:** add the citations. _[FIXED — refs added.]_

5. **METHODS §7 misstates the fit.** [literature] It says the **reactivation rate `A_react`** sets
   the rebound clock; P5 says (and derives) that **`f_lat`** does and that the fit *sidesteps* the
   reactivation-rate dispute. The consolidated methods doc contradicts the procedure. **Fix.** _[FIXED.]_

## 🟡 MEDIUM — confirmed, qualify the claims

6. **The P8 "synergy" is semi-definitional.** [logic + adversary + code] `resv_scale` multiplies
   the rebound *seed* pool `Llat` while `kf` multiplies the *kill rate*; in a Poisson-extinction
   model fewer-seeds × faster-killing is super-additive **by construction**. "Neither lever alone
   >40%, together ~70–96%" restates `P(control)≈exp(−Llat·scale/f(kf))`. The *direction* is robust;
   "synergy as discovery" is not. (The grid also drops the TIP entirely.) **Fix:** reword to a
   structural property; note `verify_claims` only tests the trivial direction.
7. **P4's headline "neutral" largely means the TIP was *absent*.** [logic] The P4 bolus washes out
   in ~1 day and the latent TIP reactivates too rarely — genuine *engagement* only arrives in P4b.
   The README compresses P4+P4b, importing P4b's credibility into the P4 table. **Fix:** separate them.
8. **`verify_claims` neutrality tolerances are near-vacuous.** [code] "TIP neutral" passes for
   `|ΔControl|` anywhere in `[0,20]` pts (noise at N=120 is several pts), while "helps" is gated
   tightly (`synergy≥10`). Asymmetric — the headline verdict is essentially untestable there. **Fix:** tighten / reframe.
9. **ACTG wk12 rebound fraction is ~83%, cited as 77%.** [literature] PMC4911279 gives 194/235 ≈
   82.6% at confirmed ≥200 c/mL. **Fix the number in P5_findings.** _[FIXED.]_
10. **P4b varies every knob except the structural ones.** [adversary #5] Dosing/ψ/ν/reactivation
    move; the compartment decoupling (#2) and floor (#3) never do — the same "parameter ≠ structural
    robustness" limit AUDIT.md named and never escaped. **Fix:** acknowledge in the robustness framing.

## 🟢 LOW — note

11. **"Orthogonal" rounds a small consistent *negative* TIP nudge to zero.** [logic] Evidence is
    "neutral-to-slightly-harmful" (−1 pt P6; −0.3..−0.65 kf P4b; 61→51% co-reactivating stress).
    **Fix:** soften the README banner word. _[FIXED.]_
12. **RIO "24% (7/29) vs 6% placebo" mixes denominators.** [literature] 7/29 (ATI sub-analysis) vs
    a /34 placebo; primary endpoint is 7/34≈21% vs 2/34≈6%. HR 0.09 is exact. **Fix:** clarify.
13. **Control metric (`Iw+Id<50`) has a built-in sign asymmetry.** [adversary #6] A TIP only *adds*
    `Id`, so in rebound it can only push *over* threshold → "slight harm, never help" is partly the
    metric. **Fix:** note as a metric caveat.
14. **`tip_model.py` hardcodes `/home/seth/...` output paths** (lines ~121/136/149) — the recurring
    portability bug; phase modules use `HERE`. **Fix.** _[FIXED.]_
15. **`tip_model_p8_boosters.py:run()` drops the `ν` weighting** (`Estar(Iw+Id+Ldef)`); fine since
    P8 has no TIP (ν=1), but an undocumented divergence from the shared kernel. **Note.**

## What was checked and UPHELD

- **Numerics:** R0=8.70, TIP invasion threshold ~7.5, QSS reductions, and the stochastic engine
  (tau-leaping propensities, `np.minimum` superinfection caps, the latency binomial split,
  conservation, clipping, per-replicate `kf`/`f_lat` broadcasting, seeded determinism) all read
  **clean**; `verify_claims.py` re-ran **14/14** under cap (caveat #1: against a model copy).
- **P2 "antagonism not by-construction (0/169)"** is genuinely a decoupled-pool test and holds.
- **Clinical calibration numbers are accurate and correctly attributed:** Gunst 2025 (16/21/32 d),
  A5345 (22 d @≥1000), CHAMP (4%/13%), Siliciano (t½ 44 mo), RIO HR 0.09, Hill-vs-Davenport ~30×
  reactivation span — all verified against primary sources. No "Bar 2016" misattribution exists.
- **The honesty discipline is real:** the illustrative-params and functional-vs-sterilizing caveats
  *are* in the README banner, not buried — the audit's job here was to find where *framing* outran
  the experiment, and it did, but the underlying disclosures are largely present.

## The single load-bearing assumption

**That a reactivating reservoir virus never transits a TIP-accessible co-infected state** (#2),
reinforced by **a maintained antigen floor that keeps CD8 primed regardless of the TIP** (#3).
Together they manufacture "no substrate where the cure succeeds; can't catch the rebound where it
fails." The verdict is *one compartment-coupling assumption deep*. The defensible claim is the
**conditional** one: *given a non-coupled reservoir and a maintained antigen floor, a TIP is inert
to post-treatment control* — narrower than the headline.

## The shared blind spot

Every phase, every prior audit, and `verify_claims.py` run the **same** ODE skeleton, R0=8.70,
one-clone exponential-reactivation reservoir, memoryless QSS immunity, and non-interacting
`L_lat`/`L_latd`. **No audit ever instantiated a structurally different model** (coupled-TIP,
dynamic-exhaustion immunity, clonal/spatial reservoir) and asked whether neutrality survives. The
layered honesty trail is real but is **honesty within one model family** — it reads as robustness
but is partly convergence on shared priors.

## Required corrections (checklist)

- [x] `verify_claims.py` → exercise production `simulate()`/`rebound_curve()`/`run()` (#1)
- [x] METHODS §7 `A_react`→`f_lat` clock fix (#5); ACTG wk12 77%→~83% (#9)
- [x] Add cross-disease citations (#4); soften "orthogonal"→"neutral (conditional)" + name the two
      structural assumptions in the README banner (#2,#3,#11); reword P8 "synergy" (#6)
- [x] Fix `tip_model.py` hardcoded paths (#14)
- [x] **Coupled-TIP reservoir built ([P11](P11_findings.md)) — and it FLIPPED the verdict.** With
      χ = fraction of reactivations carrying the TIP, control rises monotonically 22%→58% (χ:0→1)
      at the marginal immune level: a coupled TIP **helps**. The "neutral" headline was a decoupling
      artifact, exactly as finding #2 predicted. README verdict revised.
- [ ] **Dynamic-exhaustion immunity** still unbuilt — P11 used QSS+floor immunity (CD8 can't be
      starved), so it shows the HELP side only; whether strong coupling reopens the HARM side
      (antagonism via exhaustible memory) is the remaining open test (#3, blind-spot).
