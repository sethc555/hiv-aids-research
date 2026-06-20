#!/usr/bin/env python3
"""P2b — the last mean-field variation: saturating (Holling-II) immune killing.

P2 closed the antigen-pool question (decoupling doesn't rescue the TIP) but flagged one
remaining mean-field assumption: killing is MASS-ACTION in the effector, kw = k*E, so a
reservoir-maintained CD8 pool kills the TIP's dual factories arbitrarily fast (that IS the
obstruction). P1.5 also listed Holling-II / handling-time saturation as untried.

This variant saturates the PER-CELL killing rate in effector abundance:
    kill(E) = k*E / (1 + E/Esat)
so once CD8 is abundant (reservoir-primed), adding more does NOT increase per-cell killing
of the TIP's factories. If the antagonism is purely the unbounded-predation term, saturation
should OPEN a stable immune-compatible escape. Everything else is P2 (two-pool, slow
reservoir, QSS virus, same audit discipline).
"""
import os
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tip_model import P, T0, _san
from tip_model_p13_wm import QS, Estar

HERE = os.path.dirname(os.path.abspath(__file__))
PSI = 22.0
pc = P["p"] / P["c"]
S_SEED = 1e-4
ESAT = 1.5e4          # effector half-saturation for killing (baseline E* ~ 1.9e4)


def kill(E):
    return QS["k"] * E / (1.0 + E / ESAT)          # Holling-II per-cell killing rate


def rhs(t, y, nu, R, g):
    T, Iw, It, Id, L = np.clip(y, 0.0, 1e13)
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]
    Vw = pc * (Iw + (1 - rho) * Id)
    Vt = pc * (PSI * rho * Id)
    E = Estar(Iw + nu * Id + L)
    kw, kd = kill(E), nu * kill(E)
    return [
        lam - dT * T - b * T * Vw - b * T * Vt,
        b * T * Vw - (d + kw) * Iw - b * Iw * Vt,
        b * T * Vt - d * It - b * It * Vw,
        b * Iw * Vt + b * It * Vw - (d + kd) * Id,
        g * L * (1.0 - L / R) + S_SEED * (Iw + Id),
    ]


def integ(y0, nu, R, g, tmax=3000):
    return solve_ivp(rhs, (0, tmax), _san(y0), args=(nu, R, g), method="LSODA",
                     rtol=1e-8, atol=1e-3, max_step=5.0)


def vload(y):
    return pc * (y[1] + (1 - P["rho"]) * y[3])


def tail_osc(s, win=400):
    vw = np.clip(pc * (s.y[1] + (1 - P["rho"]) * s.y[3]), 1e-6, None)
    v = vw[s.t > s.t[-1] - win]
    return v.max() / max(v.min(), 1e-9)


def main():
    print(f"P2b saturating-kill (Esat={ESAT:.0e}); baseline E* ~ {Estar(5000):.2e}")
    g = 0.02
    NU = np.linspace(0.0, 1.0, 13)
    RES = np.linspace(0.0, 15000.0, 13)
    wtred = np.zeros((len(NU), len(RES)))
    cd8ret = np.zeros((len(NU), len(RES)))
    stable = np.zeros((len(NU), len(RES)))
    A0 = None
    for j, R in enumerate(RES):
        Rs = max(R, 1.0)
        base = integ([T0, 1e-3, 0, 0, 1.0], 1.0, Rs, g)
        bend = base.y[:, -1]; vw0 = max(vload(bend), 1e-9)
        E0 = Estar(bend[1] + bend[3] + bend[4])
        if j == 0:
            A0 = bend[1] + bend[3]
        for i, nu in enumerate(NU):
            st = list(bend); st[2] += 50; st[3] += 50
            s = integ(st, nu, Rs, g)
            end = s.y[:, -1]
            wtred[i, j] = np.log10(vw0) - np.log10(max(vload(end), 1e-9))
            cd8ret[i, j] = Estar(end[1] + nu * end[3] + end[4]) / max(E0, 1e-9)
            stable[i, j] = 1.0 if tail_osc(s) < 1.05 else 0.0
    escape = (wtred >= 0.5) & (cd8ret >= 0.7) & (stable > 0.5)
    st = stable > 0.5
    print(f"baseline active antigen A0 = {A0:.0f}")
    print(f"STABLE fraction: {stable.mean():.0%}")
    print(f"max WT reduction (stable cells): {wtred[st].max() if st.any() else 0:.2f} log")
    print(f"IMMUNE-COMPATIBLE & STABLE: {escape.sum()}/{escape.size} -> "
          f"{'EXISTS' if escape.any() else 'NONE'}")
    if escape.any():
        ii, jj = np.where(escape); b = np.argmax(wtred[ii, jj])
        print(f"  best: nu={NU[ii[b]]:.2f}, R={RES[jj[b]]/A0:.1f}xA0 -> "
              f"WT down {wtred[ii[b],jj[b]]:.2f} log, CD8 {cd8ret[ii[b],jj[b]]:.0%}, STABLE")
    np.savez(os.path.join(HERE, "p2b_sweep.npz"), NU=NU, RES=RES, A0=A0,
             wtred=wtred, cd8ret=cd8ret, stable=stable)
    print("wrote p2b_sweep.npz")


if __name__ == "__main__":
    main()
