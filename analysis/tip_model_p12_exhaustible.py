#!/usr/bin/env python3
"""P12 — coupled TIP under an EXHAUSTIBLE (waning, no-floor) immune memory: help or backfire?

P11 (coupled TIP) HELPED, but used quasi-steady immunity with a constant antigen FLOOR, so CD8
could not be starved. AUDIT2's last open finding: with a waning immune MEMORY and no floor, a
coupled TIP that suppresses the rebound might starve the antigen CD8 needs (the P1 antagonism)
and BACKFIRE. The control-vs-rebound distinction is genuinely STOCHASTIC (reservoir dies out or
re-establishes) — a deterministic ODE only finds one blurred setpoint — so this is built in the
tau-leaping engine, with the QSS+floor immunity replaced by a dynamic effector pool E:

    E earned from ACTIVE antigen only (no floor), with memory:
        E += dt*( rhoE * A/(A+KA) * (1-E) - deltaE*E ),   A = Iw + nu*Id
    immune killing = kclear*E (on Iw),  nu*kclear*E (on the TIP's dual cells Id)

Knobs: chi (TIP-reservoir coupling, as P11) and nu (how VISIBLE the TIP's dual cells are).
Hypothesis: coupled TIP HELPS when nu high (its cells keep feeding CD8) but HARMS when nu low
(it handles the rebound yet starves the memory -> E wanes -> control lost).
"""
import os
import numpy as np
from tip_model import T0, P
from tip_model_p4_reservoir import pc, PSI, F_LAT, A_REACT, DL, PL

HERE = os.path.dirname(os.path.abspath(__file__))
# rhoE small = SLOW immune recovery: puts the cohort in the marginal regime (partial baseline
# control) where a TIP effect is actually visible; with fast recovery immunity self-rescues to ~98%
# control and nothing is distinguishable. (illustrative, like all params here.)
kclear, rhoE, deltaE, KA = 14.0, 0.05, 0.004, 1000.0    # kill scale; earn; memory decay (~170d); half-sat


def simulate_exh(N, chi, nu, e_boost=0.0, t_chr=500.0, t_art=300.0, t_ati=900.0, dt=0.05, seed=0):
    """Stochastic CHRONIC->ART->ATI with dynamic waning immune memory E (no floor).
    e_boost = immune memory level set at ATI onset (a 'vaccine'); 0 = whatever survived ART.
    Returns (P(control = active<50 at end), mean final E in controlled reps)."""
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
        A = Iw + nu * Id
        E = np.clip(E + dt * (rhoE * (A / (A + KA)) * (1 - E) - deltaE * E), 0.0, 1.0)
        Vw = pc * (Iw + (1 - rho) * Id); Vt = pc * (PSI * rho * Id)
        kw, kd = kclear * E, nu * kclear * E

        def po(r):
            return rng.poisson(np.clip(r, 0, None) * dt)
        Tp = po(np.full(N, lam)); Td = po(dT * T)
        iW = po(bt * T * Vw); iT = po(bt * T * Vt)
        iWs = np.minimum(po(bt * Iw * Vt), Iw.astype(np.int64)); iTs = np.minimum(po(bt * It * Vw), It.astype(np.int64))
        iWd = po((d + kw) * Iw); iTd = po(d * It); idd = po((d + kd) * Id)
        react = np.minimum(po(A_REACT * Llat), Llat.astype(np.int64))
        to_Id = rng.binomial(react.astype(np.int64), chi); to_Iw = react - to_Id    # coupling
        ldd = po(DL * Llat); lp = po(PL * Llat)
        lat = rng.binomial(np.maximum(iW.astype(np.int64), 0), F_LAT)
        T = np.clip(T + Tp - Td - iW - iT, 0, None)
        Iw = np.clip(Iw + (iW - lat) - iWd - iWs + to_Iw, 0, None)
        It = np.clip(It + iT - iTd - iTs, 0, None)
        Id = np.clip(Id + (iWs + iTs) - idd + to_Id, 0, None)
        Llat = np.clip(Llat + lat - react - ldd + lp, 0, None)
    ctrl = (Iw + Id) < 50
    return float(ctrl.mean()), float(E[ctrl].mean() if ctrl.any() else np.nan)


def main():
    N = 200
    print("P12 coupled TIP under WANING immune memory (no floor). A 'vaccine' sets E at ATI;")
    print("we compare no-TIP (chi=0) vs coupled-TIP (chi=1) across dual-cell visibility nu.\n")
    print(f"{'e_boost':>8} {'nu':>5} {'ctrl noTIP':>11} {'ctrl +cpldTIP':>14} {'verdict':>10}")
    rows = []
    for eb in [0.6, 0.45, 0.35]:                # post-ART immune memory levels (marginal regime)
        for nu in [1.0, 0.5, 0.2]:
            c0 = simulate_exh(N, 0.0, nu, e_boost=eb, t_ati=700.0, seed=int(eb * 100 + nu * 13))[0]
            c1 = simulate_exh(N, 1.0, nu, e_boost=eb, t_ati=700.0, seed=int(eb * 100 + nu * 13))[0]
            dd = 100 * (c1 - c0)
            v = "HELPS" if dd > 8 else ("HARMS" if dd < -8 else "~neutral")
            rows.append((eb, nu, c0, c1, dd)); print(f"{eb:>8.1f} {nu:>5.1f} {100*c0:>10.0f}% {100*c1:>13.0f}% {v:>10}")
    np.savez(os.path.join(HERE, "p12_exhaustible.npz"), rows=np.array(rows))
    print("\nverdict per row: HELPS = coupled TIP raises control; HARMS = lowers it (antagonism).")
    print("wrote p12_exhaustible.npz")


if __name__ == "__main__":
    main()
