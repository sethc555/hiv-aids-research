#!/usr/bin/env python3
"""DESIGN-TRIAGE ENGINE — the model's verdict map, made QUERYABLE.

Pattern borrowed from the sibling T1D repo's trial_triage.py. Given a TIP design / regime, `triage(...)`
returns a verdict (HELPS / CONDITIONAL / NEUTRAL / NO-BENEFIT) + the mechanistic reason + the fix + a
confidence tag that honours the model's limits + a trace to the result that grounds the rule. The
thresholds are the DERIVED ones imported from analytic.py (kappa_crit, psi*), not hardcoded.

Honest by construction: the model finds NO systematic backfire (P15), so HARMS/CONTRAINDICATED is never
returned; every HELPS verdict is tagged conditional on the coupling premise chi (H3) and the contested
effector-PTC mechanism (H6), and notes that the clinical KM anchors are fitted (loo_validation). Lightweight
(rule engine + the analytic setpoint/NGM; no stochastic runs). Run under the 4 GB cap by habit.
"""
import analytic as A   # KAPPA_CRIT = 7.70, psi_star = 7.49 (derived, self-validated)

KCRIT, PSTAR = A.KAPPA_CRIT, A.psi_star
MARGIN = 3.0   # marginal-immunity band below kappa_crit -- the rescue band a coupled TIP can flip
VERDICTS = ("HELPS", "CONDITIONAL", "NEUTRAL", "NO-BENEFIT")   # note: HARMS is never returned (no backfire)


def _v(verdict, mechanism, fix, confidence, trace):
    return {"verdict": verdict, "mechanism": mechanism, "fix": fix,
            "confidence": confidence, "trace": trace}


def triage(chi, kappa, psi=60.0, nu=0.9, delivery="latent", immune_maintained=True):
    """Query the model's verdict for a TIP design.
      chi      : reservoir-coupling fraction in [0,1] (fraction of reactivations carrying the TIP)
      kappa    : effector killing rate (compare to the derived control threshold kappa_crit=7.70/day)
      psi      : TIP mobilization advantage (must exceed the NGM invasion threshold psi*=7.49)
      nu       : dual-cell immune visibility (low = stealthy / evasive)
      delivery : 'latent' (reservoir-coupled/engineered) | 'sustained' (maintained dosing) | 'bolus'
      immune_maintained : whether effector priming is sustained (the antigen floor, H4)
    """
    # 1. the TIP cannot even invade
    if psi < PSTAR:
        return _v("NO-BENEFIT",
                  f"psi={psi:.1f} < the NGM invasion threshold psi*={PSTAR:.2f}: the TIP cannot persist",
                  "raise the mobilization advantage psi above ~7.5",
                  "HIGH (derived threshold)", "analytic.py H2 / verify_claims p1")
    # 2. no maintained immunity -> no control, TIP or not (the TIP assists, never replaces)
    if not immune_maintained:
        return _v("NO-BENEFIT",
                  "without maintained effector priming there is no post-treatment control, with or without a TIP",
                  "maintain immunity first; the floor is load-bearing AND contested",
                  "MEDIUM (the antigen floor is contested, H4)", "AUDIT2 #3 / assumptions H4")
    # 3. decoupled -> the de Boer limit (inert)
    if chi < 0.1:
        return _v("NEUTRAL",
                  "chi~0: a reactivating provirus never transits a TIP-accessible state -> inert (the de Boer limit)",
                  "engineer reservoir coupling (raise chi) -- BUT this is the UNPROVEN premise (no evidence a "
                  "delivered TIP reaches the pre-existing reservoir)",
                  "result SOLID; the FIX is the open Dodd & de Boer question (H3)", "P11/P14 (chi=0) / assumptions H3")
    # 4. delivery that isn't present when the reservoir rebounds
    if delivery == "bolus":
        return _v("NEUTRAL",
                  "a bolus washes out in ~1 day; the reservoir reactivates too rarely to meet it (P4)",
                  "use sustained dosing or a reservoir-coupled latent TIP so the TIP is present at rebound",
                  "HIGH", "P4 / P4b")
    # 5. coupled, but immunity far below threshold -> nothing to assist
    if kappa < KCRIT - MARGIN:
        return _v("NO-BENEFIT",
                  f"kappa={kappa:.1f} far below threshold (kappa_crit={KCRIT:.1f}): the TIP assists but cannot replace immunity",
                  "raise immune strength; the TIP is an adjunct for MARGINAL controllers, not a substitute",
                  "HIGH (P11: 0% control if immunity not maintained)", "P11 / P16")
    # 6. coupled, but immunity already wins -> no headroom
    if kappa >= KCRIT:
        return _v("NEUTRAL",
                  f"kappa={kappa:.1f} >= threshold ({KCRIT:.1f}): immunity already controls; the TIP ceilings out (no headroom)",
                  "the benefit is for MARGINAL controllers; strong controllers gain little",
                  "HIGH", "P14 / P16 (benefit concentrated at marginal immunity)")
    # 7. marginal immunity but sub-threshold coupling
    if chi < 0.6:
        return _v("CONDITIONAL",
                  f"marginal immunity (kappa={kappa:.1f}, just below {KCRIT:.1f}) but coupling chi={chi:.2f} is sub-threshold (~0.6 needed)",
                  "raise chi toward >=0.6 (substantial reservoir coupling) to clear the rescue band",
                  "MODERATE (needs strong coupling -- the contested H3)", "P14 (chi*) / P16")
    # 8. the sweet spot: strong coupling + marginal immunity
    extra = " even a stealthy/evasive TIP helps here, by holding burden down and sparing CD8 from exhaustion," if nu < 0.3 else ""
    return _v("HELPS",
              f"reservoir-coupled TIP (chi={chi:.2f}) intercepts each rebound burst at marginal immunity (R_eff~1) -- "
              f"the band just below kappa_crit where a small push flips control;{extra} (P13/P14/P16)",
              "(already favourable) -- the headline regime",
              "SOLID within the model; CONDITIONAL on chi (H3) + the contested effector-PTC mechanism (H6); "
              "the clinical KM anchors are fitted (loo)", "P11/P13/P14/P15/P16 / analytic")


# canonical TIP-design regimes (the analog of T1D running its triage on every real trial)
SCENARIOS = [
    ("standard TIP, bolus at ATI",                dict(chi=1.0, kappa=5.0, delivery="bolus")),
    ("standard TIP, sustained, DECOUPLED (chi=0)", dict(chi=0.0, kappa=5.0, delivery="sustained")),
    ("sub-threshold TIP (psi<7.5)",               dict(chi=1.0, kappa=5.0, psi=5.0)),
    ("immunity NOT maintained (no floor)",         dict(chi=1.0, kappa=5.0, immune_maintained=False)),
    ("reservoir-coupled, WEAK immunity",           dict(chi=1.0, kappa=2.0)),
    ("reservoir-coupled, STRONG controller",       dict(chi=1.0, kappa=9.0)),
    ("partial coupling, marginal (chi=0.5)",       dict(chi=0.5, kappa=5.0)),
    ("reservoir-coupled, MARGINAL controller",     dict(chi=1.0, kappa=5.0)),
    ("stealthy coupled TIP, exhaustible immunity",  dict(chi=1.0, kappa=5.0, nu=0.1)),
]


def main():
    print("DESIGN-TRIAGE ENGINE — the model's verdict on each TIP-design regime\n" + "=" * 78)
    order = {v: i for i, v in enumerate(VERDICTS)}
    rows = [(name, triage(**kw)) for name, kw in SCENARIOS]
    for name, r in rows:
        print(f"\n[{r['verdict']:>11}] {name}")
        print(f"   why : {r['mechanism']}")
        print(f"   fix : {r['fix']}")
        print(f"   conf: {r['confidence']}   <- {r['trace']}")
    helps = [n for n, r in rows if r["verdict"] == "HELPS"]
    print("\n" + "=" * 78)
    print(f"FAVOURABLE regimes ({len(helps)}): " + "; ".join(helps))
    print("No regime returns HARMS/CONTRAINDICATED — the model finds no systematic backfire (P15); the")
    print("control metric can only show a TIP pushing OVER threshold, so 'slight harm, never help' would be")
    print("partly the metric (AUDIT2 #13). Every HELPS is conditional on the coupling premise chi (H3).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
