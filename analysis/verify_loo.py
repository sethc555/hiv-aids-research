#!/usr/bin/env python3
"""verify_loo.py — assert the (honest) out-of-sample findings reproduce (exit 0 iff all pass).

Companion to loo_validation.py. Verifies the reproducible facts — including the honest NEGATIVE: HIV's
clinical KM timing anchors do NOT transfer out-of-sample from a reservoir clock alone (the heavy
late-rebound tail is immune-heterogeneity-driven), and the reduced clock is a faithful proxy of the full
stochastic model (so the negative is real, not an artifact). HIV's genuine out-of-sample evidence is the
de Boer EXTERNAL reproduction (verify_claims `deboer`), stated narrowly. Lightweight; run under the 4 GB cap.
"""
import sys
import loo_validation as L


def main():
    print("verify_loo — out-of-sample findings re-asserted:")
    checks = [
        ("reduced clock is a faithful proxy of the full P6 KM curve (<8 pts)",
         L._fidelity is not None and L._fidelity < 0.08),
        ("HONEST NEGATIVE: rebound-timing anchors do NOT transfer from a clock alone (held-out err >25%)",
         L.MEDIAN_TIMING_ERR > 0.25),
    ]
    npass = 0
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        npass += bool(ok)
    print(f"\n==== out-of-sample verification: {npass}/{len(checks)} PASS ====")
    print("  (honest read: HIV's clinical KM anchors are fitted, not mutually-predictive; the real")
    print("   out-of-sample evidence is the de Boer external reproduction — narrower than T1D's clinical LOO.)")
    if npass < len(checks):
        print("FAILED:", [n for n, ok in checks if not ok])
    return 0 if npass == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
