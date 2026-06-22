#!/usr/bin/env python3
"""P11 — the audit's decisive test: COUPLE the TIP to the reservoir (and weaken the antigen floor).

AUDIT2's #1/#2 findings: the "TIP neutral to cure" verdict is forced by two structural choices —
(1) the reactivating reservoir virus never transits a TIP-accessible co-infected state (L_lat -> Iw
only; the TIP's pool is disjoint), and (2) a constant antigen FLOOR (Ldef) keeps CD8 primed under
suppression. A real TIP is *designed* to mobilize on WT co-infection, so a reactivating provirus
SHOULD be able to carry/acquire the TIP. P11 adds a single coupling knob:

  chi = fraction of reservoir reactivations that emerge as TIP-carrying DUAL cells (Id) instead of
        pure WT cells (Iw).   chi=0 = the old decoupled model;  chi=1 = the TIP rides every rebound.

We then ask the project's question again on the calibrated cohort: as chi rises (and at a weaker
antigen floor RDEF), does the TIP-neutral verdict FLIP — to help (TIP intercepts rebound) or to
harm (TIP suppresses WT but starves the CD8 that maintains control = the P1 antagonism)? We report
BOTH post-treatment control AND CD8 retention, so we can tell help from harm.

Stochastic tau-leaping, 6 states (T, Iw, It, Id, Ldef, Llat); run memory-capped.
"""
import os
import numpy as np
from tip_model import T0, P
from tip_model_p4_reservoir import pc, PSI, F_LAT, A_REACT, DL, PL, G, RDEF, S_SEED
from tip_model_p13_wm import QS, Estar

HERE = os.path.dirname(os.path.abspath(__file__))
DETECT = 1000.0


def simulate_coupled(N, kf, chi, rdef=RDEF, nu=0.9, psi=60.0, f_lat=F_LAT,
                     t_chr=500.0, t_art=300.0, t_ati=500.0, dt=0.05, seed=0):
    """CHRONIC->ART->ATI. chi = fraction of reservoir reactivations emerging as TIP-carrying duals.
    Returns (P(control), mean CD8 retention in controlled reps)."""
    rng = np.random.default_rng(seed)
    kf = np.full(N, float(kf)) if np.isscalar(kf) else kf
    T = np.full(N, T0); Iw = np.full(N, 10.0)
    It = np.zeros(N); Id = np.zeros(N); Ldef = np.zeros(N); Llat = np.zeros(N)
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]
    k = QS["k"]
    E0 = None
    n_chr, n_art = int(t_chr / dt), int(t_art / dt)
    nsteps = n_chr + n_art + int(t_ati / dt)
    for n in range(nsteps):
        art = 1.0 if n < n_chr else (1e-3 if n < n_chr + n_art else 1.0)
        bt = b * art
        if n == n_chr + n_art:                       # record CD8 baseline at ATI onset
            E0 = Estar(Iw + nu * Id + Ldef)
        Vw = pc * (Iw + (1 - rho) * Id)
        Vt = pc * (psi * rho * Id)
        E = Estar(Iw + nu * Id + Ldef)
        kw, kd = kf * k * E, nu * kf * k * E

        def po(r):
            return rng.poisson(np.clip(r, 0, None) * dt)
        Tp = po(np.full(N, lam)); Td = po(dT * T)
        iW = po(bt * T * Vw); iT = po(bt * T * Vt)
        iWs = np.minimum(po(bt * Iw * Vt), Iw.astype(np.int64))
        iTs = np.minimum(po(bt * It * Vw), It.astype(np.int64))
        iWd = po((d + kw) * Iw); iTd = po(d * It); idd = po((d + kd) * Id)
        Lb = po(G * Ldef + S_SEED * (Iw + Id)); Ld = po(G * Ldef * Ldef / max(rdef, 1.0))
        react = np.minimum(po(A_REACT * Llat), Llat.astype(np.int64))
        to_Id = rng.binomial(react.astype(np.int64), chi)      # *** COUPLING ***
        to_Iw = react - to_Id
        ldd = po(DL * Llat); lp = po(PL * Llat)
        lat = rng.binomial(np.maximum(iW.astype(np.int64), 0), f_lat)
        T = np.clip(T + Tp - Td - iW - iT, 0, None)
        Iw = np.clip(Iw + (iW - lat) - iWd - iWs + to_Iw, 0, None)
        It = np.clip(It + iT - iTd - iTs, 0, None)
        Id = np.clip(Id + (iWs + iTs) - idd + to_Id, 0, None)
        Ldef = np.clip(Ldef + Lb - Ld, 0, None)
        Llat = np.clip(Llat + lat - react - ldd + lp, 0, None)
    prod = Iw + Id
    ctrl = prod < 50
    cd8 = Estar(Iw + nu * Id + Ldef) / np.maximum(E0, 1e-9)
    return float(ctrl.mean()), float(cd8[ctrl].mean() if ctrl.any() else np.nan)


def main():
    N = 200
    print(f"P11 coupled-TIP test (N={N}). chi = fraction of reservoir reactivations that carry the TIP.")
    print("Question: does coupling the TIP to the rebound FLIP the neutral verdict?\n")
    print("--- baseline antigen floor (RDEF=5000) ---")
    print(f"{'kf':>4} {'chi':>5} {'P(control)':>11} {'CD8 kept(ctrl)':>15}")
    rows = []
    for kf in [11.0, 13.0]:
        for chi in [0.0, 0.25, 0.5, 1.0]:
            c, cd = simulate_coupled(N, kf, chi, seed=int(kf * 10 + chi * 7))
            rows.append((5000, kf, chi, c, cd))
            print(f"{kf:>4.0f} {chi:>5.2f} {100*c:>10.0f}% {100*cd:>14.0f}%")
    print("\n--- WEAK antigen floor (RDEF=500: immunity must be earned from active antigen) ---")
    for kf in [13.0, 16.0]:
        for chi in [0.0, 0.5, 1.0]:
            c, cd = simulate_coupled(N, kf, chi, rdef=500.0, seed=int(kf * 10 + chi * 7) + 1)
            rows.append((500, kf, chi, c, cd))
            print(f"{kf:>4.0f} {chi:>5.2f} {100*c:>10.0f}% {100*cd:>14.0f}%")

    # verdict
    arr = np.array(rows)
    base = arr[arr[:, 0] == 5000]
    for kf in [11.0, 13.0]:
        sub = base[base[:, 1] == kf]
        d_ctrl = 100 * (sub[-1, 3] - sub[0, 3])      # chi=1 vs chi=0
        verdict = "HELPS" if d_ctrl > 8 else ("HARMS" if d_ctrl < -8 else "still ~neutral")
        print(f"\nkf={kf:.0f}, floor=5000: coupling chi 0->1 changes control by {d_ctrl:+.0f} pts -> {verdict}")
    np.savez(os.path.join(HERE, "p11_coupled.npz"), rows=arr)
    print("\nwrote p11_coupled.npz")


if __name__ == "__main__":
    main()
