#!/usr/bin/env python3
"""Audit gate for P1.4: is the 'escape' a converged fixed point or a transient/
oscillation (the trap that killed the spatial model)? Check convergence of Vw and Id."""
import numpy as np
from scipy.integrate import solve_ivp
from tip_model import T0, _san
import tip_model_p14 as m

for nu, Adef, tag in [(0.42, 11250, "best escape"), (0.5, 5214, "realistic A_def~1x"),
                      (1.0, 0, "no refinements (P1.3 baseline)")]:
    base = m.integ([T0, 0, 0, 0, 1e-3, 0], 1.0, Adef)
    y0 = list(base); y0[5] = 1e2
    s = solve_ivp(m.rhs, (0, 5000), _san(y0), args=(nu, Adef), method="LSODA",
                  rtol=1e-7, atol=1e-2, max_step=4.0, t_eval=np.linspace(0, 5000, 2000))
    Vw = np.clip(s.y[4], 1e-3, None)
    snap = {t: Vw[np.argmin(np.abs(s.t-t))] for t in (1000, 2000, 3000, 5000)}
    tail = Vw[s.t > 4000]
    print(f"{tag} (nu={nu}, Adef={Adef}): Vw @ "
          + ", ".join(f"{t}d={v:.2e}" for t, v in snap.items())
          + f" | tail max/min={tail.max()/max(tail.min(),1e-3):.3f} "
          + ("STABLE" if tail.max()/max(tail.min(), 1e-3) < 1.05 else "OSCILLATING"))
