#!/usr/bin/env python3
"""P8 — booster vaccines / sustained immunotherapy: the positive mirror of the TIP result.

The TIP was neutral to cure because it is PARASITIC on active virus -- it needs the very
replication the cure removes. A booster vaccine / therapeutic antibody is the OPPOSITE: it
raises immune control (kf) DIRECTLY, on the axis the cure actually turns on. So it can help
where the TIP could not. This phase quantifies HOW, on the P6-calibrated heterogeneous cohort:

  BOOSTERS   single vs repeated vs continuous dosing -> delay-only vs durable cure
  TITRATION  (1) maintenance threshold: constant kf level -> P(durable control)
             (2) cure = level x duration: hold then STOP -> bankable durable cure
             (3) trough rule: same peak, vary decay -> control tracks the TROUGH not the peak
  MORE ROOM  reservoir size x maintained immunity -> synergy (the CML 'deep response x
             duration' lever: shrink the reservoir AND boost -> where the big cures live)

A 'booster' raises immune strength kf(t) transiently; immunity wanes (decay half-life), so the
relevant quantity is the TROUGH between doses, exactly like a bNAb trough level. kf is a
dimensionless model knob (not a clinical dose): the transferable findings are the STRUCTURE
(trough-not-peak, frequency-over-amplitude, level-sets-ceiling, duration-fills-it, reservoir
synergy), not the specific percentages. Run a section: `python3 tip_model_p8_boosters.py {boost|titr|synergy|all}`.

AUDIT DISCIPLINE: stochastic, P over the cohort; 'durable control' = never rebounded through
the full post-ATI window INCLUDING after dosing stops (true treatment-free remission).
"""
import os, sys
import numpy as np
from tip_model import T0, P
from tip_model_p4_reservoir import simulate, pc, DL, PL, G, RDEF, S_SEED, A_REACT
from tip_model_p13_wm import QS, Estar
from tip_model_p6_heterogeneous import draw_cohort

HERE = os.path.dirname(os.path.abspath(__file__))
DETECT, dt, T_ATI = 1000.0, 0.05, 1200.0
kf_base, f_lat = draw_cohort(10.0, 2.5, np.random.default_rng(7))   # the P6-fit cohort
N = len(kf_base)
_art = None


def art_state():
    global _art
    if _art is None:
        s0 = np.zeros((N, 7)); s0[:, 0] = T0; s0[:, 1] = 10
        _art = simulate(s0, 1.0, kf_base, 500.0, 300.0, 0.0, f_lat=f_lat, seed=1)
    return _art


def run(level_func, resv_scale=1.0, t_ati=T_ATI, seed=2):
    """ATI from the ART state; kf(t)=max(kf_base, level_func(t)); reservoir scaled by resv_scale.
    Returns P(never rebounded through t_ati) = durable treatment-free control."""
    art = art_state()
    rng = np.random.default_rng(seed)
    T, Iw, It, Id, Ldef, Llat, Llatd = (art[:, j].astype(float).copy() for j in range(7))
    Llat = Llat * resv_scale                                      # reservoir reduction (early ART / LRA)
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]; k = QS["k"]
    treb = np.full(N, np.nan)
    for n in range(int(t_ati / dt)):
        t = n * dt; Vw = pc * (Iw + (1 - rho) * Id)
        treb[np.isnan(treb) & (Vw > DETECT)] = t
        kf = np.maximum(kf_base, level_func(t))
        E = Estar(Iw + Id + Ldef); kw = kf * k * E

        def po(r):
            return rng.poisson(np.clip(r, 0, None) * dt)
        Tp = po(np.full(N, lam)); Td = po(dT * T); iW = po(b * T * Vw); iWd = po((d + kw) * Iw)
        Lb = po(G * Ldef + S_SEED * (Iw + Id)); Ld = po(G * Ldef * Ldef / RDEF)
        lr = np.minimum(po(A_REACT * Llat), Llat.astype(np.int64)); ldd = po(DL * Llat); lp = po(PL * Llat)
        lat = rng.binomial(np.maximum(iW.astype(np.int64), 0), f_lat)
        T = np.clip(T + Tp - Td - iW, 0, None); Iw = np.clip(Iw + (iW - lat) - iWd + lr, 0, None)
        Ldef = np.clip(Ldef + Lb - Ld, 0, None); Llat = np.clip(Llat + lat - lr - ldd + lp, 0, None)
    return float(np.mean(np.isnan(treb)))


def pulses(times, amp, thalf, base=0.0):
    def f(t):
        b = base
        for tk in times:
            if t >= tk:
                b += amp * 0.5 ** ((t - tk) / thalf)
        return b
    return f


def section_boost():
    print("BOOSTERS (cohort kf~N(10,2.5); 'cure' = durable TFR through 1200d):")
    rows = [("no booster", run(lambda t: 0.0)),
            ("single boost +6 (t1/2 60d)", run(pulses([0], 6, 60))),
            ("repeat +6 q90d x2yr then STOP", run(pulses(list(np.arange(0, 720, 90)), 6, 60))),
            ("repeat +6 q45d x2yr then STOP", run(pulses(list(np.arange(0, 720, 45)), 6, 60))),
            ("repeat +6 q90d FOREVER", run(pulses(list(np.arange(0, 1200, 90)), 6, 60)))]
    for name, c in rows:
        print(f"  {name:34s} {100*c:4.0f}%")
    return rows


def section_titr():
    print("TITRATION 1 -- maintenance threshold (hold kf=L forever):")
    L1 = [(L, run(lambda t, L=L: L)) for L in [10, 12, 13, 14, 15, 16, 18]]
    for L, c in L1:
        print(f"  hold kf={L:<3} {100*c:4.0f}%")
    print("TITRATION 2 -- cure = level x duration (hold kf=15 for D then STOP):")
    L2 = [(D, run(lambda t, D=D: 15.0 if t < D else 0.0)) for D in [365, 730, 1095, 1460]]
    for D, c in L2:
        print(f"  hold 15 for {D:>5}d {100*c:4.0f}%")
    print("TITRATION 3 -- trough rule (peak 18 q120d, vary decay -> trough):")
    L3 = []
    for thalf in [20, 40, 60, 120]:
        trough = 10 + 8 * 0.5 ** (120 / thalf)
        c = run(pulses(list(np.arange(0, 1200, 120)), 8, thalf, base=10))
        L3.append((thalf, trough, c)); print(f"  t1/2={thalf:>4}d trough_kf={trough:4.1f} {100*c:4.0f}%")
    return L1, L2, L3


def section_synergy():
    print("MORE ROOM -- reservoir reduction x maintained immunity (hold kf=L forever):")
    RS = [1.0, 0.3, 0.1, 0.03]      # reservoir scale (1x, 3x, 10x, 33x reduction)
    LV = [11, 13, 15, 17]
    grid = np.zeros((len(RS), len(LV)))
    print("   reservoir\\kf   " + "  ".join(f"{L:>4}" for L in LV))
    for i, rs in enumerate(RS):
        for j, L in enumerate(LV):
            grid[i, j] = run(lambda t, L=L: L, resv_scale=rs)
        print(f"   {1/rs:>5.0f}x smaller  " + "  ".join(f"{100*grid[i,j]:3.0f}%" for j in range(len(LV))))
    return RS, LV, grid


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = {}
    if which in ("boost", "all"):
        out["boost"] = section_boost()
    if which in ("titr", "all"):
        out["titr"] = section_titr()
    if which in ("synergy", "all"):
        rs, lv, grid = section_synergy()
        np.savez(os.path.join(HERE, "p8_synergy.npz"), resv_scale=np.array(rs), levels=np.array(lv), cure=grid)
        # figure for the synergy (the headline 'more room' result)
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        im = ax.imshow(100 * grid, origin="lower", aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(lv))); ax.set_xticklabels(lv)
        ax.set_yticks(range(len(rs))); ax.set_yticklabels([f"{1/r:.0f}x" for r in rs])
        ax.set_xlabel("maintained immune strength  kf"); ax.set_ylabel("reservoir reduction (x smaller)")
        ax.set_title("P8 'more room': reservoir reduction x boosting -> durable cure %\n"
                     "(shrinking the reservoir lowers the immune level needed AND raises the ceiling)")
        for i in range(len(rs)):
            for j in range(len(lv)):
                ax.text(j, i, f"{100*grid[i,j]:.0f}%", ha="center", va="center",
                        color="w" if grid[i, j] < 0.6 else "k", fontsize=10)
        fig.colorbar(im, label="durable cure %"); fig.tight_layout()
        fig.savefig(os.path.join(HERE, "p8_synergy.png"), dpi=130)
        print("wrote p8_synergy.npz, p8_synergy.png")


if __name__ == "__main__":
    main()
