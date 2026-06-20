#!/usr/bin/env python3
"""Audit robustness sweep: is the P1.2 antagonism a tuning artifact or a structural
rule? Vary immune strength (k) and TIP design (entry bt, diversion rho); for each,
measure the audit metric via EARLY co-inoculation (single timing -> avoids the
oscillation-phase/timing confound). Control = psi=1 suppression (if a no-advantage
TIP 'suppresses' WT, that setting is immune/phase-dominated, not a TIP effect).
"""
import numpy as np
from tip_model import P
import tip_model_p12 as m

b0 = P["b"]
PSI = np.linspace(1.0, 20.0, 9)
SETTINGS = [(k, rho, bm) for k in (1e-4, 2e-4, 4e-4)
            for rho in (0.85, 0.95) for bm in (1.0, 3.0)]

print(f"{'k':>6} {'rho':>5} {'bt/b':>5} | {'base_vw':>9} {'psi1_wt':>7} "
      f"{'max_wt':>7} {'wt@cd8>=.5':>10} {'compat?':>8}")
print("-"*72)
any_compat = False
for k, rho, bm in SETTINGS:
    m.IM["k"] = k; P["rho"] = rho; P["bt"] = b0*bm
    base_vw, base_em, _ = m.tail_avg(m.ACUTE, 1.0, 0.0)
    wts, cd8s = [], []
    for psi in PSI:
        vw, em, _ = m.run_tip(psi, 0.0)            # early co-inoculation
        wts.append(np.log10(base_vw) - np.log10(vw))
        cd8s.append(em / max(base_em, 1e-9))
    wts, cd8s = np.array(wts), np.array(cd8s)
    psi1_wt = wts[0]                               # control: no-advantage TIP
    kept = cd8s >= 0.5
    max_wt_cd8kept = wts[kept].max() if kept.any() else 0.0
    compat = ((wts >= 1.0) & kept & (psi1_wt < 0.5)).any()   # genuine TIP, not immune
    any_compat = any_compat or compat
    print(f"{k:>6.0e} {rho:>5.2f} {bm:>5.1f} | {base_vw:>9.2e} {psi1_wt:>7.2f} "
          f"{wts.max():>7.2f} {max_wt_cd8kept:>10.2f} {str(compat):>8}")

print("-"*72)
print(f"ANY immune-compatible TIP across {len(SETTINGS)} parameterizations: "
      f"{'YES — antagonism is tuning-dependent' if any_compat else 'NO — antagonism robust to (k, rho, bt)'}")
print("(compat = some psi with WT down >=1 log AND CD8 >=50% AND psi=1 control <0.5 log)")
