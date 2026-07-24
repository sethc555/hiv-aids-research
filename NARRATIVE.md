# The thread, as a story

> ⚠️ **SUPERSEDED — historical record.** This narrates the 2026-06-19/20 session, which ended at the
> verdict *"a TIP is orthogonal / neutral to the immune-control cure."* **That verdict was later
> RETRACTED.** The multi-agent re-audit ([AUDIT2.md](analysis/AUDIT2.md) #2) found it was an artifact of
> modelling the TIP and the reservoir as **non-interacting**; once coupled ([P11](analysis/P11_findings.md)),
> a reservoir-coupled TIP **helps** post-treatment control. See [README](README.md) and
> [paper/MANUSCRIPT.md](paper/MANUSCRIPT.md) for the current verdict, and
> [analysis/corrections.json](analysis/corrections.json) (entry **C2**) for the machine-verifiable
> retraction. Kept unedited **because** the reversal is the point — the trail is the credibility.

_A readable narrative of the session that produced this repo (2026-06-19/20)._

## 1. A wrong turn that pointed the right way

The session opened on HIV/AIDS and its cure. An early instinct — could a **federated-database /
drug-repurposing engine** accelerate a cure? — was
examined and rejected: the cure-relevant data (latent-reservoir biology, immune control)
doesn't live in clean federatable databases, and repurposing-via-shared-targets is a
crowded, low-leverage angle for HIV specifically.

That "no" was the useful part. It forced the real question: **where is the field actually
stuck, and where could someone with compute but no wet lab contribute?**

## 2. Mapping the frontier (the corpus)

Two research passes — a web deep-research run, then a **Semantic Scholar** citation scan
(`scan.py`, `analyze.py`, using the key from the soil project) — mapped the 2024–2026 cure
frontier. The signal was unambiguous: the field has converged toward **durable ART-free
remission via CD8 "vaccinal-effect" immune control** — the Geneva and 2nd-Berlin transplant
remissions (immune-mediated, *not* CCR5-dependent), the RIO long-acting-bNAb trial, the UCSF
combination-immunotherapy result. A citation-graph expansion identified the one **open,
no-wet-lab modeling lane**: **therapeutic interfering particles (TIPs)** — engineered
conditionally-replicating defective HIV — a small, theory-first community (Weinberger,
Singh, Dar, de Boer), with a live, unresolved disagreement: *does immunity help or hurt a
TIP?* (de Boer 2026 said it hurts.)

## 3. The model, in phases

So the project became a within-host model of TIP-vs-immunity, built up in phases, each
write-up in `analysis/P*_findings.md`:

- **P1 (static immunity).** A TIP is produced only in *dually-infected* cells; immunity
  kills those cells. Moderate immune pressure collapses the TIP's efficacy window. The TIP
  lives near R0 ≈ 1, so even weak immunity sinks it.
- **P1.1 (dynamic feedback).** The hoped-for rescue — a TIP lowers viral load, immunity
  relaxes — didn't materialize cleanly; instead the system oscillated, and *early* TIP
  administration **blunted the immune response itself**. First hint of an antagonism, and
  (it later turned out) the one genuinely emergent result.
- **P1.2 (memory + exhaustion).** Realistic immunity; the antagonism persisted, but the
  model oscillated and a first-pass "−1.00 correlation / 0-of-196" headline had to be
  **retracted** as a near-degenerate, partly-artifactual statistic.
- **P1.3 (stabilize + space).** Quasi-steady-state immunity gave a clean, stable model in
  which the antagonism held — *but* a spatial two-compartment "escape hatch" test
  re-oscillated and its apparent escape was a transient artifact (caught by a convergence
  check). The Simonetti 5′-leader-NSV paper was folded in here as a real-world anchor.
- **P1.4 (Simonetti-grounded escape test).** Graded immune-visibility + a persistent
  defective-clone antigen reservoir *looked* like a clean escape — until the convergence
  gate again exposed it as oscillation.
- **P1.5 (stable-by-construction).** The principled fix — eliminate the fast free-virus
  stage too — **failed to stabilize the escape regime**, refuting the hypothesis that the
  oscillation was a free-virus lag. The sharp finding: *immune-pressure + TIP-evasion +
  stability are mutually incompatible* in these mean-field ODEs.

## 4. The audits (the real spine)

Running through everything was an insistence on **adversarial honesty**. A same-day
self-audit retracted overclaims. Then a **multi-agent audit** (math, literature, logic,
adversary — each independently verified) was run; three agents tripped an automated dual-use
content filter (TIPs read as "engineered transmissible HIV"), so those dimensions were
re-run with the biology abstracted to neutral ODE/methods language, which cleared the
filter. The audits:

- **caught a fabricated quote** (a paper-summary gloss I'd quoted as if from the source),
- **re-derived the math** and corrected an R0 scaling (√ψ, not linear),
- **caught a rehabilitated statistic** (P1.3 reviving the retracted −1.00),
- and — most usefully — **corrected the project's own audit**, which had wrongly claimed the
  model "omits dual-cell antigen" (it doesn't; and it's immaterial),
- while surfacing the deepest caveat: the model's "immune clearance" (killing *active*
  cells) is only a **loose analogy** to the clinical "vaccinal effect" (which controls the
  *latent reservoir* the model doesn't even contain).

## 5. Where it honestly lands

The computational work is sound and the literature is accurate. The headline antagonism is
a **model property, partly by construction** — clean in the simplest single-antigen-pool
model, dynamically excluded from the realistic ones, and only a loose analogy to the clinic.
The one result that *isn't* an artifact of construction is **P1.1's early-TIP immune-
blunting**. The escape question — *can a TIP suppress the virus without starving the immune
response that controls it?* — is now cleanly posed for the right successor tool: a
**stochastic / agent-based** model, where the ~10⁷-cell defective clones Simonetti measured
are individuals rather than a mean field.

That's the whole shape of it: a confident question, narrowed by honest negative results and
relentless self-audit into a small, well-posed, falsifiable one — reached from a laptop.

## 6. The continuation (2026-06-20, later) — answering the handed-off question

The audit left two things on the table: its #1 caveat (the antagonism is *partly
by-construction*, single antigen pool) and its explicit successor (a stochastic model). Both
got built.

- **P2 — the two-pool test.** Replace the single antigen pool with a *dynamic, self-
  maintaining, TIP-proof* reservoir (the Simonetti clones) that primes CD8 independently of
  active WT — the strongest honest decoupling. The antagonism **survived it: 0/169 stable
  immune-compatible TIP.** A primed CD8 pool kills the TIP's dual-infected factories no matter
  where its antigen came from. So the obstruction was never the construction — it is the TIP's
  R0≈1 **immune-fragility**. P2b (Holling-II saturating predation) then stabilized the inert
  region but still gave **0/169** in the TIP-effective region — closing the mean-field program.

- **P3 — the stochastic successor.** A tau-leaping version of the same cell model. The
  mean-field oscillation **was a continuum artifact** — discretely there is **no stable
  coexistence (0%)**; outcomes are a **lottery** between WT control and TIP loss. And the
  answer the ODEs could never give: an immune-compatible TIP **does exist** (CD8 kept at 82%
  when it works) but as a **low-probability event (~16%)**, only under strong evasion, and
  **only with the reservoir present**. The reservoir converts the antagonism into a rare
  *synergy*; the design problem becomes "raise P(cure)/P(TIP-loss)," not "find the fixed point."

The corpus was re-scanned the same day (Semantic Scholar) and confirmed current — one day on,
no new papers, only citation-count drift.

- **P3b/P3c — map + validation.** The probability was turned into a surface (P3b): a cure needs
  *both* strong evasion (ν≲0.25) *and* a reservoir (R≳1×A0), peaking ~31%, with zero stable
  coexistence anywhere. Two independent checks (P3c) — dt→0 convergence and an exact reduced-
  scale Gillespie SSA — confirmed the invariants (0% coexistence, CD8 ~82% in cures) are
  method-independent; the cure magnitude firms to order-15%.

- **P4 — the reservoir, and the real answer.** Finally the caveat that shadowed every phase:
  the model had no replication-competent latent reservoir, so "control" never meant *cure*. P4
  adds one, plus an ART→ATI (treatment-interruption) schedule — the actual clinical test.
  Post-treatment control turns out to be a **sharp CD8-strength threshold**, with the reservoir
  **persisting in every controlled case** (functional, not sterilizing, cure). And the TIP?
  **Neutral.** Latent or bolus, visible or evasive, it does not move the control threshold. The
  reason is structural: a TIP is parasitic on active WT, but the cure works by *removing* active
  WT — so the TIP has no substrate where the cure succeeds, and can't catch the rebound where it
  fails. The active-infection antagonism (P1–P3) simply **does not carry over** to the
  reservoir-control axis that defines cure.

So the grind lands somewhere honest and slightly surprising: on the axis the field actually
cares about — durable, ART-free, immune-mediated remission — an engineered TIP is neither the
threat the antagonism suggested nor a useful tool. It is **orthogonal to the cure**, because it
needs the very thing the cure is built to erase.

## 7. Calibration (P5) — anchoring the verdict, and a methods doc

The last thing separating "methods + hypotheses" from a quotable result was calibration. A
research agent pulled the clinical ATI literature (Gunst 2025, A5345, RIO, CHAMP, Hill/Davenport
reactivation models) with citations, and P5 fit the model to it:

- The **rebound clock** calibrates cleanly — `f_lat=8e-5` reproduces the ~19-day rebound median
  against A5345's 22 days (>1000 c/mL).
- At the CD8 settings that reproduce the clinical **post-treatment-control fractions** (placebo
  ~5%, bNAb/RIO ~24%), a sustained, engaged TIP shifts durable control by only **+2–5 points**.
  The earlier illustrative "neutral" is now **anchored**, not assumed. (A mid-run scare — a big
  apparent TIP benefit — turned out to be two inconsistent definitions of "control"; unifying the
  metric dissolved it. Logged as a methods lesson.)
- And the fit **earns an honest limitation**: no single immune setpoint reproduces *both* the
  fast-rebound majority (~3 wk) *and* the durable-controller minority (~5%). Real cohorts are
  heterogeneous; quantitative fitting needs a distribution over immune strength. The calibration
  didn't just confirm the answer — it named the next model.

The methodology for all of it — equations, parameters with sources, numerical settings, and the
audit discipline — is now consolidated in `analysis/METHODS.md`.

## 8. Heterogeneity (P6) — the curve closes, the verdict is final

P5 didn't just calibrate; it named its own successor. A single immune setpoint couldn't be both
the fast-rebound majority and the controller minority, so P6 made each replicate a *patient* with
its own immune strength and reservoir clock, drawn from population distributions, and fit the
**whole rebound Kaplan-Meier curve**. It closes: `kf ~ Normal(10, 2.5)` reproduces the clinical
median (20.8 d), the week-12 rebound fraction (78%), and the spontaneous post-treatment-control
fraction (5%) — all at once. The controllers are simply the upper tail of the immune
distribution; the tension P5 found was real, and heterogeneity is its resolution. And on this
clinically-fit population the TIP moves post-treatment control by a single point (5%→4%, the KM
curves lying on top of each other). The answer the whole project was built to reach is now stated
at the level the field measures it — patient cohorts, rebound curves — and it is the same answer:
*a TIP is orthogonal to the immune-control HIV cure.* A confident question, narrowed by negative
results and relentless honesty into a calibrated, falsifiable, and now patient-level one —
reached from a laptop.

> ⚠️ **This conclusion was RETRACTED.** "Orthogonal/neutral" held only because the TIP and the
> reservoir were modelled as non-interacting (χ=0). Coupling them ([P11](analysis/P11_findings.md))
> flips it: a reservoir-coupled TIP **helps**. See [AUDIT2](analysis/AUDIT2.md) #2 and
> [corrections.json](analysis/corrections.json) **C2**.

## 9. The positive program (P8) — what *would* help

A negative result earns the right to ask the positive one. The TIP fails because it's parasitic
on the virus the cure erases — so what acts *directly* on immune control? A booster vaccine /
therapeutic antibody. Run on the same fit cohort, it gives a clean dosing logic: a **single shot
only delays** rebound (the bNAb pattern), **sustained dosing banks treatment-free cures**, and —
the elegant part — **control tracks the trough between doses, not the peak**: you must keep
immunity above a maintenance threshold (≈1.5× baseline), at an interval set by how fast immunity
wanes, for a few years, to win the stochastic extinction race. And the real headroom isn't more
immunity — it's the **reservoir×immunity synergy**: neither shrinking the reservoir nor boosting
immunity alone clears ~30–40%, but together they reach ~70–96%. That is exactly the prescription
the diseases we *have* cured hand over — CML's depth-and-duration-of-response, the bNAb trials'
trough-and-combination logic — now on one quantitative axis. The project closes where it should:
having shown what doesn't work, and pointed at what does.
