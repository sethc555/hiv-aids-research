#!/usr/bin/env python3
"""Decisive control: run the P1.2 machinery with immunity fully OFF (k=0). If the TIP
then achieves multi-log suppression (like the P1 static model at kap=0), the P1.2 TIP
implementation is sound and the suppression-failure under immunity is genuinely immune-
driven. If the TIP still fails with k=0, P1.2 is broken (oscillation/averaging or a bug)
and its antagonism claim must be retracted.
"""
import numpy as np
from tip_model import P
import tip_model_p12 as m

m.IM["k"] = 0.0; P["rho"] = 0.9; P["bt"] = P["b"]          # immunity OFF
base_vw, base_em, _ = m.tail_avg(m.ACUTE, 1.0, 0.0)
print(f"immunity OFF (k=0): baseline time-avg Vw={base_vw:.2e} (expect ~7.7e5, uncontrolled)")
print(f"{'psi':>5} {'late_wt':>9}   (log10 WT reduction; P1 static gave ~1.4 @psi12, ~3.8 @psi20)")
for psi in [8, 12, 16, 20]:
    vl, _, _ = m.run_tip(psi, 400.0)
    print(f"{psi:>5} {np.log10(base_vw)-np.log10(vl):>9.2f}")
