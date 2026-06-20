#!/usr/bin/env python3
"""Verify the spatial 'escape' is a real fixed point, not a transient/oscillation
artifact. Check convergence (Vw_P at increasing times) and tail oscillation for the
claimed escape cell and neighbors.
"""
import numpy as np
from scipy.integrate import solve_ivp
from tip_model import _san
import tip_model_p13_spatial as sp

for phi, g in [(0.5, 0.03), (0.0, 0.03), (0.25, 0.03), (0.5, 0.10)]:
    # TIP run: to setpoint then inject, long integration with tail sampling
    s1 = solve_ivp(sp.rhs, (0, 400), _san(sp.ACUTE), args=(22.0, phi, g),
                   method="LSODA", rtol=1e-7, atol=1e-2, max_step=4.0)
    y0 = _san(s1.y[:, -1]); y0[5] = 1e2
    s = solve_ivp(sp.rhs, (0, 4000), y0, args=(22.0, phi, g), method="LSODA",
                  rtol=1e-7, atol=1e-2, max_step=4.0, t_eval=np.linspace(0, 4000, 2000))
    VwP = np.clip(s.y[4], 1e-3, None)
    snap = {t: VwP[np.argmin(np.abs(s.t-t))] for t in (500, 1000, 2000, 3000, 4000)}
    tail = VwP[s.t > 3000]
    print(f"phi={phi}, g={g}: Vw_P @ "
          + ", ".join(f"{t}d={v:.2e}" for t, v in snap.items())
          + f" | tail max/min={tail.max()/max(tail.min(),1e-3):.2f} "
          + ("STABLE" if tail.max()/max(tail.min(), 1e-3) < 1.1 else "OSCILLATING/UNCONVERGED"))
