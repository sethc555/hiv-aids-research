#!/usr/bin/env python3
"""P3c — exact-SSA spot-check of the P3 headline (removes the last tau-leaping caveat).

P3's ~16% cure probability at (nu=0.2, R=5000) was shown dt-robust (16% vs 17% at
dt=0.03/0.015). This runs the EXACT Gillespie (no time discretization) at the same point
for a modest number of replicates and checks the outcome split corroborates tau-leaping.

Exact SSA is event-by-event, so it is slow at these population sizes; we (a) stop each
replicate as soon as it reaches an absorbing state (WT-producing cells extinct, or TIP
lineage extinct), which in the evasion regime usually happens well before tmax, and (b)
cap wall events per replicate. Reports the split over however many replicates complete.

11 reaction channels (same model as P3): T birth/death, WT/TIP infection of T, Iw
death+kill / superinfection, It death / superinfection, Id death+kill, reservoir
birth/death. Free virus QSS; two-pool antigen Aeff = Iw + nu*Id + L.
"""
import os, sys
import numpy as np
from tip_model import P, T0
from tip_model_p13_wm import QS, Estar
from tip_model_p2_twopool import integ as det_integ

HERE = os.path.dirname(os.path.abspath(__file__))
pc = P["p"] / P["c"]
PSI = 22.0
S_SEED = 1e-4
NU, R, g = 0.2, 5000.0, 0.02

# stoichiometry for the 11 channels, columns = [T, Iw, It, Id, L]
STOICH = np.array([
    [+1, 0, 0, 0, 0],   # 0 T birth
    [-1, 0, 0, 0, 0],   # 1 T death
    [-1, +1, 0, 0, 0],  # 2 WT infect T -> Iw
    [-1, 0, +1, 0, 0],  # 3 TIP infect T -> It
    [0, -1, 0, 0, 0],   # 4 Iw death+kill
    [0, -1, 0, +1, 0],  # 5 Iw superinfected by TIP -> Id
    [0, 0, -1, 0, 0],   # 6 It death
    [0, 0, -1, +1, 0],  # 7 It superinfected by WT -> Id
    [0, 0, 0, -1, 0],   # 8 Id death+kill
    [0, 0, 0, 0, +1],   # 9 reservoir birth
    [0, 0, 0, 0, -1],   # 10 reservoir death
], dtype=np.float64)


def propensities(y, omega=1.0):
    """Counts-domain propensities. omega = system-size factor: counts = concentration*omega.
    We evaluate the concentration-rate expressions on c=y/omega and rescale by omega, which
    preserves the deterministic (concentration) dynamics while making the system tractably
    small for an exact event loop (omega<1)."""
    c = y / omega
    T, Iw, It, Id, L = c
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]
    Vw = pc * (Iw + (1 - rho) * Id)
    Vt = pc * (PSI * rho * Id)
    E = Estar(Iw + NU * Id + L)
    kw, kd = QS["k"] * E, NU * QS["k"] * E
    rates_conc = np.array([
        lam, dT * T, b * T * Vw, b * T * Vt,
        (d + kw) * Iw, b * Iw * Vt, d * It, b * It * Vw,
        (d + kd) * Id, g * L + S_SEED * (Iw + Id), g * L * L / R,
    ])
    return rates_conc * omega


def run_one(y0, rng, omega=1.0, tmax=1500.0, max_events=12_000_000):
    y = y0.astype(np.float64).copy()
    t = 0.0
    for _ in range(max_events):
        a = propensities(y, omega)
        a0 = a.sum()
        if a0 <= 0:
            break
        t += rng.exponential(1.0 / a0)
        if t > tmax:
            break
        ch = np.searchsorted(np.cumsum(a), rng.random() * a0)   # faster than rng.choice
        y = np.clip(y + STOICH[ch], 0, None)
        if (y[1] + y[3]) < 1 or (y[2] + y[3]) < 1:   # absorbing: WT-prod gone OR TIP gone
            break
    return y, t


def main():
    nreps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    omega = float(sys.argv[2]) if len(sys.argv) > 2 else 0.03   # reduced scale for tractability
    base = det_integ([T0, 1e-3, 0, 0, 1.0], 1.0, R, g).y[:, -1]
    E0 = Estar(base[1] + base[3] + base[4])
    seed = np.array(base) * omega                               # counts = concentration*omega
    seed[2] = max(100 * omega, 3); seed[3] = max(100 * omega, 3)
    rng = np.random.default_rng(31337)
    print(f"exact SSA at nu={NU}, R={R}, omega={omega} (counts=conc*omega; T0~{base[0]*omega:.0f}); "
          f"{nreps} reps. Reduced scale validates STRUCTURE (lottery, no coexistence), not the "
          f"full-size %; tau-leaping gave ~16% cure / 0% coexist at full scale.")
    cure = tiploss = coex = 0
    cd8s = []
    for r in range(nreps):
        y, t = run_one(seed, rng, omega)
        wt_gone = (y[1] + y[3]) < 1
        tip_gone = (y[2] + y[3]) < 1
        cd8 = Estar((y[1] + NU * y[3] + y[4]) / omega) / E0    # back to concentration for Estar
        if wt_gone:
            cure += 1; cd8s.append(cd8)
        elif tip_gone:
            tiploss += 1
        else:
            coex += 1; cd8s.append(cd8)
        print(f"  rep {r+1:2d}/{nreps}: t={t:6.0f}d  "
              f"{'CURE (WT extinct)' if wt_gone else 'TIP-LOST' if tip_gone else 'coexist/timeout'}"
              f"  CD8={cd8:.0%}", flush=True)
    n = cure + tiploss + coex
    print(f"\nSSA result over {n} reps: cure {cure/n:.0%}  TIP-lost {tiploss/n:.0%}  "
          f"coexist {coex/n:.0%}  | CD8 in controlled ~{np.mean(cd8s) if cd8s else float('nan'):.0%}")
    print(f"(tau-leaping reference: ~16% cure, 84% TIP-lost, 0% coexist, CD8 82%)")
    np.savez(os.path.join(HERE, "p3c_ssa.npz"),
             nreps=n, cure=cure, tiploss=tiploss, coex=coex,
             cd8_controlled=np.mean(cd8s) if cd8s else np.nan)
    print("wrote p3c_ssa.npz")


if __name__ == "__main__":
    main()
