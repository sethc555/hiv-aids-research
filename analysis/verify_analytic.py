#!/usr/bin/env python3
"""verify_analytic.py — assert the derived closed forms reproduce the model (exit 0 iff all pass).

Companion to analytic.py; mirrors the verify_*.py pattern. Re-asserts each closed form independently of
the report. Lightweight (deterministic setpoint + 3x3 NGM linear algebra + committed-grid load; no
stochastic runs) — run under the 4 GB cap by habit:
    bash -c 'ulimit -v 4194304; timeout 595 python3 verify_analytic.py'
"""
import sys
import analytic as A   # importing computes the values quietly (report only prints when run as a script)


def main():
    print("verify_analytic — derived closed forms re-asserted against the model:")
    checks = [
        ("R0 = b*T0*p/(d*c) = 8.70 (exact)",            abs(A.R0 - 8.70) < 0.05),
        ("kappa_crit = (R0-1)d = 7.70 (exact)",          abs(A.KAPPA_CRIT - 7.70) < 0.05),
        ("R_eff(kappa_crit) = 1 exactly",                abs(A.reff(A.KAPPA_CRIT) - 1.0) < 1e-9),
        ("NGM TIP invasion threshold psi* ~ 7.49",       abs(A.psi_star - 7.5) < 0.5),
        ("R0_TIP sub-linear: ratio(20/7.5) in (1.3,2.0)", 1.3 < A.R0_tip(20.0) / A.R0_tip(7.5) < 2.0),
        ("Delta(chi) >= 0 and increasing (p14 grid)",    A.OUT[4]),
    ]
    npass = 0
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        npass += bool(ok)
    print(f"\n==== analytic verification: {npass}/{len(checks)} PASS ====")
    if npass < len(checks):
        print("FAILED:", [n for n, ok in checks if not ok])
    return 0 if npass == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
