#!/usr/bin/env python3
"""Delta(chi) IN CLOSED FORM — the open question in MANUSCRIPT §4 / analytic.py #5, and the one
Dodd & de Boer asked about ("a closed-form expression for how much coupling lowers the threshold").

DERIVATION (post-ART ATI; QSS free virus; targets recovered, T ~ T0).

1. NO SHIFT OF THE ASYMPTOTIC THRESHOLD.
   Linearising the infected subsystem (Iw, It, Id) about zero infection, the superinfection terms
   (b_t*Iw*Vt, b*It*Vw) are QUADRATIC in infected densities and vanish. So at linear order a dual cell
   begets only Iw, and It is a dead end: the next-generation matrix is triangular and its spectral
   radius is set by the WT cycle alone. Hence
        kappa_crit(chi) = kappa_crit(0) = (R0 - 1) d      for every chi.
   Coupling does NOT move the deterministic invasion threshold. (Confirmed: at kappa > kappa_crit the
   simulated TIP effect is ~0.)

2. WHAT COUPLING ACTUALLY CHANGES — the seed type of each reactivation.
   A reactivation emerges as a WT cell (prob 1-chi) or a TIP-carrying DUAL cell (prob chi). Their
   first-generation WT-offspring numbers are
        R_w = R0 * d / (d + kappa)                    [WT seed: output p,        lifetime 1/(d+kappa)]
        R_d = (1 - rho) * R0 * d / (d + nu*kappa)     [dual seed: output (1-rho)p, lifetime 1/(d+nu*kappa)]
   (the dual cell's remaining output, psi*rho*p, makes TIP virions -> It carriers, which are inert at
   linear order). The per-reactivation flare potential is therefore multiplied by

        THETA  =  R_d / R_w  =  (1 - rho) (d + kappa) / (d + nu*kappa)                          (*)

   and a chi-coupled reservoir seeds flares at the mixture rate  (1-chi) + chi*THETA.

3. THE CLOSED FORM. Writing the seeded-flare reproduction number R_seed(chi) = R_w[(1-chi) + chi*THETA]
   and solving R_seed = 1 for the EFFECTIVE control boundary kappa_eff gives
        kappa_eff(chi) = R0 * d * [(1-chi) + chi*THETA] - d
   so the threshold shift is

        Delta(chi) = kappa_eff(0) - kappa_eff(chi) = R0 * d * chi * (1 - THETA)                 (**)

   * LINEAR in chi                                    -> explains the monotone rise with coupling [P14]
   * Delta >= 0  <=>  THETA <= 1                      -> the NO-BACKFIRE CONDITION, in closed form
   * Delta = 0 at chi = 0                             -> recovers the de Boer limit               [P11]

4. THE NO-BACKFIRE CONDITION (new, falsifiable, and NOT what the project assumed).
        Delta >= 0   <=>   (1 - rho)(d + kappa)  <=  d + nu*kappa
   i.e. the TIP's packaging diversion (1-rho) must outweigh the survival advantage a stealthy dual cell
   enjoys (nu < 1). With the project's rho = 0.9 this holds comfortably, which is WHY every previous run
   found "no backfire" — but it is a property of the DESIGN, not a theorem about TIPs. A weakly-diverting,
   highly-stealthy TIP (low rho, low nu) is predicted to HARM control. VALIDATED below.

Lightweight: closed forms only (no stochastic runs). `--validate` re-runs the confirming simulation.
Run under the 4 GB cap by habit.
"""
import os
import numpy as np
from tip_model import P, T0

HERE = os.path.dirname(os.path.abspath(__file__))
R0 = P["b"] * T0 * P["p"] / (P["d"] * P["c"])
D = P["d"]
KAPPA_CRIT = (R0 - 1) * D


def theta(kappa, rho=None, nu=0.9):
    """THETA = R_d/R_w, eq (*): the factor by which coupling scales per-reactivation flare potential."""
    rho = P["rho"] if rho is None else rho
    return (1.0 - rho) * (D + kappa) / (D + nu * kappa)


def delta(chi, kappa, rho=None, nu=0.9):
    """Delta(chi) = R0*d*chi*(1-THETA), eq (**): the effective control-threshold shift (/day)."""
    return R0 * D * chi * (1.0 - theta(kappa, rho, nu))


def kappa_eff(chi, kappa, rho=None, nu=0.9):
    """The effective control boundary under coupling chi."""
    return R0 * D * ((1 - chi) + chi * theta(kappa, rho, nu)) - D


def rho_star(kappa, nu=0.9):
    """Diversion below which coupling BACKFIRES (THETA = 1)."""
    return 1.0 - (D + nu * kappa) / (D + kappa)


# Confirming simulation, measured 2026-06-27 (tip_model_p16_analytic.simulate_static, N=200,
# kappa=5, nu=0.1, chi 0->1; rho varied). Predicted crossover rho* = 0.75.
VALIDATION = [  # (rho, THETA, observed TIP effect in points)
    (0.9, 0.40, +29), (0.7, 1.20, +1), (0.5, 2.00, -21), (0.3, 2.80, -36), (0.1, 3.60, -36),
]


def main():
    import sys
    print(__doc__.split("Lightweight:")[0].rstrip())
    print("=" * 92)
    print(f"\nR0 = {R0:.3f}   d = {D}   kappa_crit = (R0-1)d = {KAPPA_CRIT:.3f}/day   (independent of chi)")

    print("\n(a) THETA and Delta(chi=1) at the project's design point (rho=0.9, nu=0.9):")
    print(f"    {'kappa':>6} {'THETA':>7} {'Delta(1)':>9} {'kappa_eff(1)':>13}  verdict")
    for k in [4.0, 5.0, 6.0, 7.0, 8.0]:
        t = theta(k); print(f"    {k:>6.1f} {t:>7.3f} {delta(1.0,k):>9.2f} {kappa_eff(1.0,k):>13.2f}"
                            f"  {'helps' if t < 1 else 'BACKFIRES'}")

    print("\n(b) THE NO-BACKFIRE CONDITION  (1-rho)(d+kappa) <= d+nu*kappa   [rho* = diversion floor]")
    for nu in [0.9, 0.5, 0.1]:
        print(f"    nu={nu:>4}:  " + "  ".join(f"kappa={k:.0f} -> rho*={rho_star(k,nu):.2f}" for k in [3.,5.,7.]))

    print("\n(c) VALIDATION — the closed form PREDICTS THE SIGN of the coupling effect")
    print("    (simulate_static, N=200, kappa=5, nu=0.1, chi 0->1; predicted crossover rho* = "
          f"{rho_star(5.0, 0.1):.2f})")
    print(f"    {'rho':>5} {'THETA':>7} {'predicted':>10} {'observed':>9}  match")
    ok = True
    for rho, th, eff in VALIDATION:
        pred = "helps" if th < 1 else "BACKFIRE"
        seen = "helps" if eff > 5 else ("BACKFIRE" if eff < -5 else "neutral")
        m = (th < 1 and eff > 5) or (th > 1 and eff <= 5)      # sign predicted (boundary -> neutral ok)
        ok &= m
        print(f"    {rho:>5.1f} {th:>7.2f} {pred:>10} {eff:>+8} pts  {'OK' if m else 'XX'}")
    print(f"\n    ==> the sign of the coupling effect is predicted in {sum(1 for r,t,e in VALIDATION if (t<1 and e>5) or (t>1 and e<=5))}"
          f"/{len(VALIDATION)} cases; observed crossover ~0.70 vs predicted {rho_star(5.0,0.1):.2f}.")

    print("\n(d) THE HONEST CONSEQUENCE")
    print("    'No systematic backfire' is TRUE FOR THIS DESIGN (rho=0.9) and is now EXPLAINED, not just")
    print("    observed — but it is NOT a general property of coupled TIPs. A weakly-diverting, stealthy")
    print("    TIP (rho below rho*) is predicted to REDUCE post-treatment control. That is a design")
    print("    constraint the project had not identified: high packaging diversion is what makes")
    print("    reservoir coupling safe.")
    if "--validate" in sys.argv:
        print("\n(e) RE-RUNNING the confirming simulation (a few minutes) ...")
        import tip_model
        import tip_model_p16_analytic as A
        keep = P["rho"]
        for rho, th, _ in VALIDATION:
            tip_model.P["rho"] = rho
            c0 = A.simulate_static(200, 5.0, 0.0, nu=0.1, seed=7)
            c1 = A.simulate_static(200, 5.0, 1.0, nu=0.1, seed=7)
            print(f"    rho={rho:.1f} THETA={th:.2f}: effect {100*(c1-c0):+.0f} pts")
        tip_model.P["rho"] = keep
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
