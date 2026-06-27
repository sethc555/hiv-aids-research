#!/usr/bin/env python3
"""verify_corrections.py — re-derive the REVERSALS in the retraction trail (corrections.json), so each
withdrawn overclaim stays machine-verifiable: the OLD claim is reproducibly WRONG and the corrected claim
HOLDS. Pattern borrowed from the sibling T1D repo's `v1-withdrawn` correction-trail claim -- the old result
is kept reproducible precisely so the reversal cannot quietly disappear.

Lightweight (analytic NGM + the committed coupling grid; no stochastic runs). Run under the 4 GB cap by habit.
Exit 0 iff every reproduced reversal passes.
"""
import json
import os
import sys
import numpy as np
import analytic as A

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    cat = json.load(open(os.path.join(HERE, "corrections.json")))
    C = cat["corrections"]
    repro = [c["id"] for c in C if c["reproduced"]]
    docd = [c["id"] for c in C if not c["reproduced"]]
    print(f"verify_corrections — retraction trail: {len(C)} corrections "
          f"({len(repro)} machine-reproduced: {', '.join(repro)}; {len(docd)} documented: {', '.join(docd)})\n")

    eff = np.load(os.path.join(HERE, "p14_coupling_phase.npz"))["effect"]   # rows = chi, cols = kf

    checks = [
        # C1 -- the withdrawn linear 0.134*psi is reproducibly wrong; the corrected sqrt(psi) holds
        ("C1  withdrawn linear 0.134*psi WRONG; sqrt(psi) holds (ratio(20/7.5)~1.63; linear over-predicts >40%)",
         A.R0_tip(20.0) / A.R0_tip(7.5) < 2.0 and (0.134 * 20.0) / A.R0_tip(20.0) > 1.4),
        # C2 -- the withdrawn 'neutral' verdict was a chi=0 artifact; chi>0 helps (the reversal)
        ("C2  'neutral' holds ONLY at chi=0; coupled chi=1 HELPS (mean effect > +5 pts)",
         100 * float(np.mean(eff[-1])) > 5.0),
        # C3 -- withdrawn absolute 'never backfires'; corrected 'no SYSTEMATIC backfire' (min within noise)
        ("C3  no SYSTEMATIC backfire: grid min effect >= -2 pts (within noise; absolute 'never' was imprecise)",
         100 * float(eff.min()) >= -2.0),
    ]
    npass = 0
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        npass += bool(ok)
    print(f"\n==== correction-trail verification: {npass}/{len(checks)} reproduced reversals PASS ====")
    print("  (the documented framing/citation corrections C4/C5 are recorded in corrections.json; the trail is the credibility.)")
    if npass < len(checks):
        print("FAILED:", [n.split()[0] for n, ok in checks if not ok])
    return 0 if npass == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
