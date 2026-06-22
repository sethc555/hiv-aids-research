#!/usr/bin/env python3
"""P13 — the harshest test: ACTIVE immune damage (true exhaustion) + a very STEALTHY coupled TIP.

P12 found no backfire under a merely-*waning* immune memory. The remaining corner (P12's own
caveat): an immune system that is ACTIVELY DEGRADED by viral burden (true T-cell exhaustion, not
just under-stimulation), combined with a TIP whose dual cells EVADE killing (very low visibility).
The feared mechanism: a stealthy coupled TIP converts a clearable rebound into a PERSISTENT
smolder of evasive dual cells that the immune system can't clear yet that keeps exhausting it ->
CD8 degrades -> control is lost. If a backfire exists anywhere, it should be here.

Two changes vs P12:
  - maintenance is epitope-specific (nu-weighted):  Aeff = Iw + nu*Id  feeds E
  - exhaustion is driven by TOTAL burden (evadable killing does NOT spare you from inflammation):
        dE includes  - exhaust * burden/(burden+Kx) * E ,   burden = Iw + Id
  so evasive dual cells (low nu) EXHAUST immunity (via burden) without FEEDING it (low Aeff).
Sweep exhaustion strength x very-low visibility nu; compare no-TIP vs fully-coupled TIP.
"""
import os
import numpy as np
from tip_model import T0, P
from tip_model_p4_reservoir import pc, PSI, F_LAT, A_REACT, DL, PL

HERE = os.path.dirname(os.path.abspath(__file__))
kclear, rhoE, deltaE, KA = 14.0, 0.05, 0.004, 1000.0
Kx = 2000.0                                   # burden at half-max exhaustion


def simulate_dmg(N, chi, nu, exhaust, e_boost=0.5, t_chr=500.0, t_art=300.0, t_ati=900.0,
                 dt=0.05, seed=0):
    """Stochastic CHRONIC->ART->ATI with ACTIVE exhaustion. Returns (P(control), mean final E)."""
    rng = np.random.default_rng(seed)
    T = np.full(N, T0); Iw = np.full(N, 10.0); It = np.zeros(N); Id = np.zeros(N); Llat = np.zeros(N)
    E = np.full(N, 0.05)
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]
    n_chr, n_art = int(t_chr / dt), int(t_art / dt)
    for n in range(n_chr + n_art + int(t_ati / dt)):
        art = 1.0 if n < n_chr else (1e-3 if n < n_chr + n_art else 1.0)
        bt = b * art
        if n == n_chr + n_art and e_boost:
            E = np.maximum(E, e_boost)
        Aeff = Iw + nu * Id                          # feeds CD8 (epitope-specific)
        burden = Iw + Id                             # exhausts CD8 (inflammation, unevadable)
        E = np.clip(E + dt * (rhoE * (Aeff / (Aeff + KA)) * (1 - E)
                              - deltaE * E
                              - exhaust * (burden / (burden + Kx)) * E), 0.0, 1.0)
        Vw = pc * (Iw + (1 - rho) * Id); Vt = pc * (PSI * rho * Id)
        kw, kd = kclear * E, nu * kclear * E

        def po(r):
            return rng.poisson(np.clip(r, 0, None) * dt)
        Tp = po(np.full(N, lam)); Td = po(dT * T)
        iW = po(bt * T * Vw); iT = po(bt * T * Vt)
        iWs = np.minimum(po(bt * Iw * Vt), Iw.astype(np.int64)); iTs = np.minimum(po(bt * It * Vw), It.astype(np.int64))
        iWd = po((d + kw) * Iw); iTd = po(d * It); idd = po((d + kd) * Id)
        react = np.minimum(po(A_REACT * Llat), Llat.astype(np.int64))
        to_Id = rng.binomial(react.astype(np.int64), chi); to_Iw = react - to_Id
        ldd = po(DL * Llat); lp = po(PL * Llat)
        lat = rng.binomial(np.maximum(iW.astype(np.int64), 0), F_LAT)
        T = np.clip(T + Tp - Td - iW - iT, 0, None)
        Iw = np.clip(Iw + (iW - lat) - iWd - iWs + to_Iw, 0, None)
        It = np.clip(It + iT - iTd - iTs, 0, None)
        Id = np.clip(Id + (iWs + iTs) - idd + to_Id, 0, None)
        Llat = np.clip(Llat + lat - react - ldd + lp, 0, None)
    ctrl = (Iw + Id) < 50
    return float(ctrl.mean()), float(E.mean())


def main():
    N = 200
    print("P13 harshest test: ACTIVE immune exhaustion (burden-driven) + STEALTHY coupled TIP.")
    print("Looking for a BACKFIRE: does a coupled TIP (chi=1) LOWER control vs no-TIP (chi=0)?\n")
    print(f"{'exhaust':>8} {'nu':>5} {'ctrl noTIP':>11} {'ctrl +cpldTIP':>14} {'dControl':>9} {'verdict':>9}")
    rows = []
    for exhaust in [0.05, 0.15, 0.30]:
        for nu in [0.2, 0.1, 0.05]:
            c0 = simulate_dmg(N, 0.0, nu, exhaust, seed=int(exhaust * 200 + nu * 17))[0]
            c1 = simulate_dmg(N, 1.0, nu, exhaust, seed=int(exhaust * 200 + nu * 17))[0]
            dd = 100 * (c1 - c0)
            v = "HELPS" if dd > 8 else ("BACKFIRE" if dd < -8 else "~neutral")
            rows.append((exhaust, nu, c0, c1, dd))
            print(f"{exhaust:>8.2f} {nu:>5.2f} {100*c0:>10.0f}% {100*c1:>13.0f}% {dd:>+8.0f} {v:>9}")
    np.savez(os.path.join(HERE, "p13_exhaustion_damage.npz"), rows=np.array(rows))
    anyback = any(r[4] < -8 for r in rows)
    print(f"\n=> BACKFIRE found somewhere: {anyback}. "
          f"(verdict per row: HELPS / BACKFIRE / ~neutral). wrote p13_exhaustion_damage.npz")


if __name__ == "__main__":
    main()
