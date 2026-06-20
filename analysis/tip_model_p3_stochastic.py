#!/usr/bin/env python3
"""P3 — the stochastic successor: does discreteness resolve the escape the ODEs couldn't?

The whole ODE program (P1-P1.5, P2 antigen-decoupling, P2b Holling-II killing) reached the
same wall: the ONLY regime where the TIP suppresses WT is dual-cell immune evasion (nu<1),
and there the mean-field equations OSCILLATE -- never a stable fixed point, so "does an
immune-compatible TIP exist?" was undecidable in the continuum. P1.5/P2 handed that question
to the right tool: a STOCHASTIC model, where (a) the limit-cycle troughs that the ODE rides
through become real EXTINCTION events, and (b) the Simonetti defective clones are discrete
individuals, not a mean field.

This is a demographic (tau-leaping) version of the P2 cell model: the same 5 integer
populations (T, Iw, It, Id, L), QSS free virus, two-pool antigen (Aeff = Iw + nu*Id + L),
slow self-maintaining reservoir. We run many replicates in the nu<1 regime that oscillated
in mean-field and ask: what does the oscillation BECOME when populations are discrete?

Hypothesis: the mean-field limit cycle is replaced by a LOTTERY between two absorbing
outcomes -- (i) WT-producing cells go extinct in a trough -> functional control ("cure"),
or (ii) the TIP lineage (It,Id) goes extinct -> WT rebounds ("TIP lost") -- with little
stable coexistence. If so, the "immune-compatible TIP" is not a fixed point you design for
but a PROBABILITY you'd have to win -- a fundamentally different (and more honest) framing.

AUDIT DISCIPLINE: report the full outcome DISTRIBUTION over replicates (not a single run);
classify by absorbing state; report CD8 retention only among WT-controlled survivors.
"""
import os
import numpy as np
from tip_model import P, T0, _san
from tip_model_p13_wm import QS, Estar
from tip_model_p2_twopool import integ as det_integ          # reuse P2 deterministic burn-in

HERE = os.path.dirname(os.path.abspath(__file__))
PSI = 22.0
pc = P["p"] / P["c"]
S_SEED = 1e-4


def tau_leap(state0, nu, R, g, tmax=1500.0, dt=0.03, seed=0, record=False):
    """Vectorized tau-leaping over many replicates. state0: (5,) -> broadcast to (Nrep,5)."""
    rng = np.random.default_rng(seed)
    Nrep = state0.shape[0]
    T, Iw, It, Id, L = (state0[:, k].astype(np.float64).copy() for k in range(5))
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]
    k = QS["k"]
    nsteps = int(tmax / dt)
    rec_t, rec = ([], []) if record else (None, None)
    for n in range(nsteps):
        Vw = pc * (Iw + (1 - rho) * Id)
        Vt = pc * (PSI * rho * Id)
        E = Estar(Iw + nu * Id + L)
        kw, kd = k * E, nu * k * E
        # propensities (per day) -> expected events over dt
        def pois(rate):
            return rng.poisson(np.clip(rate, 0, None) * dt)
        Tprod = pois(np.full(Nrep, lam)); Tdeath = pois(dT * T)
        infW = pois(b * T * Vw); infT = pois(b * T * Vt)
        IwD = pois((d + kw) * Iw); IwS = pois(b * Iw * Vt)
        ItD = pois(d * It); ItS = pois(b * It * Vw)
        IdD = pois((d + kd) * Id)
        Lb = pois(g * L + S_SEED * (Iw + Id)); Ld = pois(g * L * L / R)
        # apply (clip at 0; consume superinfection only up to available cells)
        IwS = np.minimum(IwS, Iw); ItS = np.minimum(ItS, It)
        T = np.clip(T + Tprod - Tdeath - infW - infT, 0, None)
        Iw = np.clip(Iw + infW - IwD - IwS, 0, None)
        It = np.clip(It + infT - ItD - ItS, 0, None)
        Id = np.clip(Id + IwS + ItS - IdD, 0, None)
        L = np.clip(L + Lb - Ld, 0, None)
        if record and n % int(2 / dt) == 0:                   # ~every 2 days
            rec_t.append(n * dt); rec.append(pc * (Iw + (1 - rho) * Id).copy())
    final = np.stack([T, Iw, It, Id, L], axis=1)
    return (final, (np.array(rec_t), np.array(rec))) if record else final


def classify(final, nu, E0):
    """Per-replicate absorbing-state classification + CD8 retention."""
    T, Iw, It, Id, L = (final[:, k] for k in range(5))
    Vw = pc * (Iw + (1 - P["rho"]) * Id)
    cd8 = Estar(Iw + nu * Id + L) / np.maximum(E0, 1e-9)
    wt_gone = (Iw + Id) < 1                       # no WT-producing cells -> functional control
    tip_gone = (It + Id) < 1                       # TIP lineage extinct
    coexist = (~wt_gone) & (~tip_gone) & (Vw < 1e4)   # both present AND WT suppressed
    wt_persist = (~wt_gone) & tip_gone            # TIP lost, WT back
    return dict(wt_gone=wt_gone, tip_gone=tip_gone, coexist=coexist, wt_persist=wt_persist,
                Vw=Vw, cd8=cd8)


def main():
    g, R, Nrep = 0.02, 5000.0, 200
    # burn-in: WT + reservoir to setpoint (deterministic), then seed a TIP bolus
    base = det_integ([T0, 1e-3, 0, 0, 1.0], 1.0, R, g).y[:, -1]
    vw0 = pc * (base[1] + (1 - P["rho"]) * base[3])
    E0 = Estar(base[1] + base[3] + base[4])
    print(f"burn-in WT setpoint Vw0={vw0:.2e}, Iw={base[1]:.0f}, reservoir L={base[4]:.0f}, "
          f"E*0={E0:.2e}  (R={R:.0f}, g={g})")
    seed_state = np.array(base); seed_state[2] = 100; seed_state[3] = 100   # TIP carriers + dual

    print(f"\n{Nrep} stochastic replicates per nu, tmax=1500d "
          f"(mean-field OSCILLATES for every nu<1 here):")
    print(f"{'nu':>5} | {'WT-controlled':>13} {'TIP-lost':>9} {'coexist':>8} | "
          f"{'CD8 kept (controlled)':>22}")
    rows = []
    for ni, nu in enumerate([0.2, 0.4, 0.6, 0.8, 1.0]):
        s0 = np.tile(seed_state, (Nrep, 1))
        final = tau_leap(s0, nu, R, g, seed=1000 + ni)
        c = classify(final, nu, E0)
        controlled = c["wt_gone"] | c["coexist"]
        cd8_ctrl = c["cd8"][controlled].mean() if controlled.any() else float("nan")
        print(f"{nu:>5.1f} | {c['wt_gone'].mean():>12.0%} {c['wt_persist'].mean():>8.0%} "
              f"{c['coexist'].mean():>7.0%} | {cd8_ctrl:>21.0%}")
        rows.append((nu, c["wt_gone"].mean(), c["wt_persist"].mean(), c["coexist"].mean(), cd8_ctrl))

    arr = np.array(rows)
    np.savez(os.path.join(HERE, "p3_outcomes.npz"),
             nu=arr[:, 0], wt_controlled=arr[:, 1], tip_lost=arr[:, 2],
             coexist=arr[:, 3], cd8_kept=arr[:, 4], R=R, g=g, Nrep=Nrep)
    print("\nKey question: in the controlled outcomes, is CD8 still primed (immune-compatible),")
    print("and is control a STABLE design point or a stochastic coin-flip? (see distribution above)")
    print("wrote p3_outcomes.npz")


if __name__ == "__main__":
    main()
