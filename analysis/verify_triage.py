#!/usr/bin/env python3
"""verify_triage.py — assert the design-triage engine reproduces the model's verdict on every canonical
TIP-design regime (exit 0 iff all pass). Companion to triage.py. Lightweight (rule engine + the analytic
setpoint/NGM; no stochastic runs) — run under the 4 GB cap by habit.
"""
import sys
from triage import triage, VERDICTS

CASES = [
    ("decoupled chi=0 -> de Boer limit",         dict(chi=0.0, kappa=5.0),                          "NEUTRAL"),
    ("coupled + MARGINAL immunity -> helps",      dict(chi=1.0, kappa=5.0),                          "HELPS"),
    ("coupled + STRONG immunity -> no headroom",  dict(chi=1.0, kappa=9.0),                          "NEUTRAL"),
    ("sub-threshold psi -> cannot invade",        dict(chi=1.0, kappa=5.0, psi=5.0),                 "NO-BENEFIT"),
    ("immunity not maintained -> no control",     dict(chi=1.0, kappa=5.0, immune_maintained=False), "NO-BENEFIT"),
    ("bolus -> washes out before rebound",        dict(chi=1.0, kappa=5.0, delivery="bolus"),        "NEUTRAL"),
    ("partial coupling chi=0.5 -> conditional",   dict(chi=0.5, kappa=5.0),                          "CONDITIONAL"),
    ("coupled + WEAK immunity -> can't replace",  dict(chi=1.0, kappa=2.0),                          "NO-BENEFIT"),
    ("stealthy coupled + marginal -> helps",      dict(chi=1.0, kappa=5.0, nu=0.1),                  "HELPS"),
]


def main():
    print("verify_triage — design-triage reproduces the model's verdict on each regime:")
    npass = 0
    for name, kw, expect in CASES:
        got = triage(**kw)["verdict"]
        ok = got == expect
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {got}  (expect {expect})")
        npass += ok
    no_harm = all(triage(**kw)["verdict"] in VERDICTS for _, kw, _ in CASES)
    print(f"  [{'PASS' if no_harm else 'FAIL'}] no regime returns HARMS/CONTRAINDICATED "
          f"(the model finds no systematic backfire, P15)")
    npass += no_harm
    total = len(CASES) + 1
    print(f"\n==== triage verification: {npass}/{total} PASS ====")
    if npass < total:
        print("FAILED:", [n for n, kw, e in CASES if triage(**kw)["verdict"] != e])
    return 0 if npass == total else 1


if __name__ == "__main__":
    sys.exit(main())
