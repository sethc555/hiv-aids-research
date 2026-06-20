#!/usr/bin/env python3
"""Audit check: is the P1.2 TIP underperforming because of IMMUNITY or because of the
early-co-inoculation PROTOCOL? At weak immunity (k=1e-4, immunity barely controls WT),
compare early co-inoc vs late introduction at the established setpoint (the P1 protocol).
If late >> early, the robustness sweep's 'TIP barely works' is a protocol artifact, and
P1.2 cannot cleanly attribute TIP failure to immune antagonism.
"""
import numpy as np
from tip_model import P
import tip_model_p12 as m

m.IM["k"] = 1e-4; P["rho"] = 0.9; P["bt"] = P["b"]
base_vw, base_em, sb = m.tail_avg(m.ACUTE, 1.0, 0.0)
# is immunity even active at k=1e-4?
print(f"k=1e-4 baseline: time-avg Vw={base_vw:.2e} (no-immunity ref 7.7e5), "
      f"CD8 E+M={base_em:.2e}")
print(f"{'psi':>5} {'early_wt':>9} {'late_wt':>9}   (log10 WT reduction)")
for psi in [3, 8, 12, 16, 20]:
    ve, _, _ = m.run_tip(psi, 0.0)              # early co-inoculation
    vl, _, _ = m.run_tip(psi, 400.0)            # late, at established setpoint (P1 protocol)
    print(f"{psi:>5} {np.log10(base_vw)-np.log10(ve):>9.2f} "
          f"{np.log10(base_vw)-np.log10(vl):>9.2f}")
