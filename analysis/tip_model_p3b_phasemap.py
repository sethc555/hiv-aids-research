#!/usr/bin/env python3
"""P3b — the stochastic phase map: P(cure)/P(TIP-loss) over (reservoir R, visibility nu).

P3 profiled ONE point (R=5000, nu=0.2) in depth and found the immune-compatible TIP is a
~16% outcome. P3b turns that into a SURFACE: the deterministic phase diagrams (P1, P1.5, P2)
plotted a log-reduction or a yes/no escape; the honest stochastic analogue is a PROBABILITY
of each absorbing outcome. The design objective named in P3 -- "raise P(cure)/P(TIP-loss)"
-- is exactly this surface.

Efficiency: the entire grid (nu x R x replicates) is flattened into ONE vectorized
tau-leaping call (per-replicate nu and R arrays), so it is a single step-loop over big
arrays, not a Python loop over grid cells. Reuses P3's tau_leap/classify unchanged (they
already broadcast nu and R).

AUDIT DISCIPLINE: probabilities over replicates; coexistence reported separately (P3 found
it ~0 -- if the map shows a coexistence island that would be a genuine new finding).
"""
import os
import numpy as np
from tip_model import P, T0
from tip_model_p13_wm import Estar
from tip_model_p2_twopool import integ as det_integ
from tip_model_p3_stochastic import tau_leap, classify, pc

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    g = 0.02
    NU = np.linspace(0.1, 1.0, 7)
    RES = np.array([0.0, 2500.0, 5000.0, 7500.0, 10000.0])
    REPS = 80
    tmax, dt = 1200.0, 0.04

    # per-R deterministic burn-in (baseline setpoint + E0 depend on reservoir)
    seeds, E0s, A0 = {}, {}, None
    for R in RES:
        Rs = max(R, 1.0)
        b = det_integ([T0, 1e-3, 0, 0, 1.0], 1.0, Rs, g).y[:, -1]
        s = np.array(b); s[2] = 100; s[3] = 100
        seeds[R] = s
        E0s[R] = Estar(b[1] + b[3] + b[4])
        if R == 0.0:
            A0 = b[1] + b[3]

    # flatten (nu, R, rep) -> one big replicate array with per-replicate nu, R, E0
    state0, nu_arr, R_arr, E0_arr, idx = [], [], [], [], []
    for i, nu in enumerate(NU):
        for j, R in enumerate(RES):
            Rs = max(R, 1.0)
            state0.append(np.tile(seeds[R], (REPS, 1)))
            nu_arr.append(np.full(REPS, nu)); R_arr.append(np.full(REPS, Rs))
            E0_arr.append(np.full(REPS, E0s[R])); idx.append(np.full(REPS, i * len(RES) + j))
    state0 = np.vstack(state0)
    nu_arr = np.concatenate(nu_arr); R_arr = np.concatenate(R_arr)
    E0_arr = np.concatenate(E0_arr); idx = np.concatenate(idx)
    print(f"phase map: {len(NU)}x{len(RES)} grid x {REPS} reps = {state0.shape[0]} replicates, "
          f"tmax={tmax}, dt={dt}, A0={A0:.0f}")

    final = tau_leap(state0, nu_arr, R_arr, g, tmax=tmax, dt=dt, seed=2024)
    c = classify(final, nu_arr, E0_arr)

    ncell = len(NU) * len(RES)
    p_cure = np.zeros(ncell); p_tiploss = np.zeros(ncell)
    p_coex = np.zeros(ncell); cd8_ctrl = np.full(ncell, np.nan)
    for cell in range(ncell):
        m = idx == cell
        controlled = (c["wt_gone"] | c["coexist"]) & m
        p_cure[cell] = c["wt_gone"][m].mean()
        p_tiploss[cell] = c["wt_persist"][m].mean()
        p_coex[cell] = c["coexist"][m].mean()
        if controlled.any():
            cd8_ctrl[cell] = c["cd8"][controlled].mean()
    shape = (len(NU), len(RES))
    p_cure, p_tiploss, p_coex, cd8_ctrl = (a.reshape(shape) for a in (p_cure, p_tiploss, p_coex, cd8_ctrl))

    np.savez(os.path.join(HERE, "p3b_phasemap.npz"),
             NU=NU, RES=RES, A0=A0, p_cure=p_cure, p_tiploss=p_tiploss,
             p_coex=p_coex, cd8_ctrl=cd8_ctrl, REPS=REPS)

    print("\nP(cure) = P(WT-producing cells extinct), rows=nu, cols=R/A0:")
    print("        " + "  ".join(f"{r/A0:4.1f}x" for r in RES))
    for i, nu in enumerate(NU):
        print(f"nu={nu:.2f} " + "  ".join(f"{p_cure[i,j]:4.0%}" for j in range(len(RES))))
    print(f"\nmax P(cure) on grid: {p_cure.max():.0%} at "
          f"nu={NU[np.unravel_index(p_cure.argmax(),shape)[0]]:.2f}, "
          f"R={RES[np.unravel_index(p_cure.argmax(),shape)[1]]/A0:.1f}xA0")
    print(f"stable-coexistence cells (P>5%): {(p_coex>0.05).sum()}  (P3 expected ~0)")

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(17, 5)); x = RES / A0
    for a, M, ttl, cm, vlim in [
        (ax[0], p_cure, "P(cure): WT controlled", "viridis", (0, 1)),
        (ax[1], p_tiploss, "P(TIP lost, WT rebounds)", "magma", (0, 1)),
        (ax[2], cd8_ctrl, "CD8 retained in cures", "RdYlGn", (0, 1.2))]:
        im = a.pcolormesh(x, NU, M, shading="auto", cmap=cm, vmin=vlim[0], vmax=vlim[1])
        a.set_xlabel("reservoir antigen R (x active A0)"); a.set_ylabel("dual-cell visibility nu")
        a.set_title(ttl); fig.colorbar(im, ax=a)
    fig.suptitle("P3b stochastic phase map: the escape is a probability surface "
                 "(cure needs strong evasion AND a reservoir)")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "p3b_phasemap.png"), dpi=130)
    print("wrote p3b_phasemap.npz, p3b_phasemap.png")


if __name__ == "__main__":
    main()
