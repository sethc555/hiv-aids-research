#!/usr/bin/env python3
"""verify_delta.py — assert the closed-form Delta(chi) (analytic_delta.py) holds and predicts the
simulated sign of the coupling effect. Lightweight (closed forms + the stored validation table; the
confirming simulation is re-runnable via `analytic_delta.py --validate`). 4 GB cap by habit.
"""
import sys
import analytic_delta as D


def main():
    print("verify_delta — closed-form Delta(chi) re-asserted:")
    k, nu, rho = 5.0, 0.9, 0.9
    checks = [
        # the derivation's own consistency
        ("Delta(0) = 0 (recovers the de Boer limit)", abs(D.delta(0.0, k)) < 1e-12),
        ("Delta is LINEAR in chi", abs(D.delta(1.0, k) - 2 * D.delta(0.5, k)) < 1e-9),
        ("kappa_eff(0) = kappa_crit = (R0-1)d (threshold unshifted at chi=0)",
         abs(D.kappa_eff(0.0, k) - D.KAPPA_CRIT) < 1e-9),
        # THETA and the no-backfire condition
        ("THETA < 1 at the design point (rho=0.9) -> coupling helps", D.theta(k, rho, nu) < 1.0),
        ("Delta > 0 at the design point", D.delta(1.0, k, rho, nu) > 0),
        ("no-backfire condition is exactly THETA<=1: Delta>=0 iff (1-rho)(d+kap)<=d+nu*kap",
         (D.delta(1.0, k, rho, nu) >= 0) == ((1 - rho) * (D.D + k) <= D.D + nu * k)),
        # the backfire regime the closed form predicts
        ("a weak-diversion stealthy TIP is predicted to BACKFIRE (rho=0.3, nu=0.1)",
         D.theta(k, 0.3, 0.1) > 1.0 and D.delta(1.0, k, 0.3, 0.1) < 0),
        ("predicted crossover rho*(kappa=5, nu=0.1) ~ 0.75", abs(D.rho_star(5.0, 0.1) - 0.75) < 0.02),
        # validation against the confirming simulation (stored)
        ("closed form predicts the SIGN of the simulated effect in all validation cases",
         all((t < 1 and e > 5) or (t > 1 and e <= 5) for _, t, e in D.VALIDATION)),
    ]
    npass = 0
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        npass += bool(ok)
    print(f"\n==== Delta(chi) verification: {npass}/{len(checks)} PASS ====")
    if npass < len(checks):
        print("FAILED:", [n for n, ok in checks if not ok])
    return 0 if npass == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
