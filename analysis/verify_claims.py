#!/usr/bin/env python3
"""verify_claims.py — machine-checkable verification of the whole claim chain (P1 -> P8).

Each load-bearing claim in the README / findings docs is INDEPENDENTLY re-derived here from
the phase modules and asserted within a tolerance. Stochastic claims use reduced replicate
counts (faster, wider tolerances) -- they check the DIRECTION / ORDER of magnitude that the
findings rest on, not the last digit. Exit code 0 iff every claim passes.

Run: `python3 verify_claims.py`   (a few minutes; deterministic via fixed seeds)

This is the backbone for "verify everything along the claim chain": every headline number maps
to a check here, and a regression in any phase script trips a FAIL.
"""
import sys
import numpy as np

CHECKS = []


def check(name, got, lo, hi, detail=""):
    ok = (got >= lo) and (got <= hi)
    CHECKS.append((ok, name, got, lo, hi, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: got={got:.4g} expect[{lo:g},{hi:g}] {detail}")
    return ok


# ---------- P1: deterministic core ----------
def verify_p1():
    print("P1 — within-host core (deterministic, exact):")
    from tip_model import P, T0, wt_setpoint, tip_arm
    R0 = P["b"] * T0 * P["p"] / (P["d"] * P["c"])
    check("R0 (WT basic reproductive ratio) = 8.70", R0, 8.65, 8.75, "[tip_model.py]")

    # TIP invasion threshold psi*: below it the TIP cannot suppress WT, above it it can.
    base = wt_setpoint(0.0)
    drops = {}
    for psi in [5, 6, 7, 8, 9, 10]:
        end, _ = tip_arm(base, 0.0, float(psi), tmax=500)
        drops[psi] = np.log10(base[4] / max(end[4], 1e-9))
    thr = next((psi for psi in [5, 6, 7, 8, 9, 10] if drops[psi] > 0.1), 99)
    check("TIP invasion threshold psi* ~ 7.5", thr, 6, 9, "(ODE crossing; AUDIT NGM root 7.49)")

    # immunity collapses the TIP benefit (P1 headline). NB: P1 endpoints oscillate and their
    # precise magnitudes were RETRACTED by the audit -- so verify the ROBUST claims: the
    # suppression a supra-threshold TIP ACHIEVES (tail-minimum, phase-robust), and that immunity
    # collapses it. Not a fixed endpoint magnitude (that would re-assert a retracted number).
    def supp(psi, kap):
        sp = wt_setpoint(kap)
        _, s = tip_arm(sp, kap, float(psi), tmax=700)
        vw_min = np.nanmin(np.clip(s.y[4][s.t > s.t[-1] - 300], 1e-9, None))
        return np.log10(max(sp[4], 1e-9) / vw_min)
    s_lo, s_hi = supp(20, 0.0), supp(20, 0.8)
    check("supra-threshold TIP suppresses WT at kap=0 (>1 log, tail-min)", s_lo, 1.0, 6.0)
    check("immunity collapses TIP suppression (kap0 - kap0.8 > 0.5)", s_lo - s_hi, 0.5, 6.0,
          f"(kap0 {s_lo:.2f} vs kap0.8 {s_hi:.2f} log)")


# ---------- P3: stochastic resolution ----------
def verify_p3():
    print("P3 — stochastic successor (reduced reps; nu=0.2, R=5000):")
    import numpy as np
    from tip_model_p3_stochastic import tau_leap, classify, det_integ
    from tip_model import T0
    from tip_model_p13_wm import Estar
    g, R, REPS = 0.02, 5000.0, 120
    base = det_integ([T0, 1e-3, 0, 0, 1.0], 1.0, R, g).y[:, -1]
    E0 = Estar(base[1] + base[3] + base[4])
    sd = np.array(base); sd[2] = 100; sd[3] = 100
    f = tau_leap(np.tile(sd, (REPS, 1)), 0.2, R, g, tmax=1500, dt=0.04, seed=42)
    c = classify(f, 0.2, E0)
    ctrl = c["wt_gone"] | c["coexist"]
    check("P3 stable coexistence = 0%", 100 * c["coexist"].mean(), -0.1, 0.1, "(continuum artifact)")
    check("P3 cure probability ~order-15%", 100 * c["wt_gone"].mean(), 5, 30)
    check("P3 CD8 retained in cures ~82%", 100 * (c["cd8"][ctrl].mean() if ctrl.any() else 0), 70, 95)


# ---------- shared cohort burn-in for P4/P5/P6/P8 ----------
_BURN = {}


def cohort_burn(N=120, m_kf=10.0, s_kf=2.5, seed=7):
    key = (N, m_kf, s_kf, seed)
    if key not in _BURN:
        from tip_model import T0
        from tip_model_p4_reservoir import simulate
        rng = np.random.default_rng(seed)
        kf = np.clip(rng.normal(m_kf, s_kf, N), 4, 40)
        f_lat = 10.0 ** rng.normal(np.log10(8e-5), 0.35, N)
        s0 = np.zeros((N, 7)); s0[:, 0] = T0; s0[:, 1] = 10
        art = simulate(s0, 1.0, kf, 500.0, 300.0, 0.0, f_lat=f_lat, seed=1)
        _BURN[key] = (kf, f_lat, art)
    return _BURN[key]


def _ati_control(art, kf_func, resv_scale=1.0, t_ati=900.0, dt=0.05, tip_sustained=0.0,
                 psi=60.0, nu=1.0, seed=2, detect=1000.0):
    """Generic ATI from a burn-in tuple (kf_base, f_lat, ART_state).
    Returns (P(durable control), median rebound day). kf_func(t, kf_base) -> kf array."""
    from tip_model import P
    from tip_model_p4_reservoir import pc, DL, PL, G, RDEF, S_SEED, A_REACT
    from tip_model_p13_wm import QS, Estar
    kf_base, f_lat, A = art
    T, Iw, It, Id, Ldef, Llat, Llatd = (A[:, j].astype(float).copy() for j in range(7))
    Llat = Llat * resv_scale
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]; k = QS["k"]
    N = len(kf_base); rng = np.random.default_rng(seed); treb = np.full(N, np.nan)
    for n in range(int(t_ati / dt)):
        t = n * dt; Vw = pc * (Iw + (1 - rho) * Id); treb[np.isnan(treb) & (Vw > detect)] = t
        if tip_sustained:
            It = It + rng.poisson(np.full(N, tip_sustained * dt))
        kf = np.maximum(kf_base, kf_func(t, kf_base))
        Vt = pc * (psi * rho * Id)
        E = Estar(Iw + nu * Id + Ldef); kw, kd = kf * k * E, nu * kf * k * E

        def po(r):
            return rng.poisson(np.clip(r, 0, None) * dt)
        Tp = po(np.full(N, lam)); Td = po(dT * T)
        iW = po(b * T * Vw); iT = po(b * T * Vt)
        iWs = np.minimum(po(b * Iw * Vt), Iw.astype(np.int64)); iTs = np.minimum(po(b * It * Vw), It.astype(np.int64))
        iWd = po((d + kw) * Iw); iTd = po(d * It); idd = po((d + kd) * Id)
        Lb = po(G * Ldef + S_SEED * (Iw + Id)); Ld = po(G * Ldef * Ldef / RDEF)
        lr = np.minimum(po(A_REACT * Llat), Llat.astype(np.int64)); ldd = po(DL * Llat); lp = po(PL * Llat)
        ldr = np.minimum(po(A_REACT * Llatd), Llatd.astype(np.int64)); lddd = po(DL * Llatd); ldp = po(PL * Llatd)
        lat = rng.binomial(np.maximum(iW.astype(np.int64), 0), f_lat)
        T = np.clip(T + Tp - Td - iW - iT, 0, None)
        Iw = np.clip(Iw + (iW - lat) - iWd - iWs + lr, 0, None)
        It = np.clip(It + iT - iTd - iTs, 0, None)
        Id = np.clip(Id + (iWs + iTs) - idd + ldr, 0, None)
        Ldef = np.clip(Ldef + Lb - Ld, 0, None)
        Llat = np.clip(Llat + lat - lr - ldd + lp, 0, None)
        Llatd = np.clip(Llatd - ldr - lddd + ldp, 0, None)
    reb = treb[~np.isnan(treb)]
    med = np.median(reb) if len(reb) else np.inf
    return float(np.mean(np.isnan(treb))), med


def verify_p4_p6():
    print("P4/P5/P6 — reservoir, calibration, cohort (reduced cohort N=120):")
    art = cohort_burn()
    # P4: control is a SHARP threshold in kf (homogeneous-ish via constant high/low)
    c_lo, _ = _ati_control(art, lambda t, kb: 8.0, seed=11)
    c_hi, _ = _ati_control(art, lambda t, kb: 18.0, seed=12)
    check("P4 control sharp: low immunity -> low control", 100 * c_lo, 0, 25)
    check("P4 control sharp: high immunity -> high control", 100 * c_hi, 35, 100, "(threshold behaviour)")
    # P4 TIP neutrality: no-TIP vs sustained TIP at a mid level -> within tolerance
    c_no, _ = _ati_control(art, lambda t, kb: kb, seed=21)
    c_tip, _ = _ati_control(art, lambda t, kb: kb, tip_sustained=2000.0, psi=60.0, nu=0.9, seed=22)
    check("P4 TIP ~neutral (|dControl| small)", 100 * abs(c_tip - c_no), 0, 20, f"(noTIP {100*c_no:.0f}%, +TIP {100*c_tip:.0f}%)")
    # P6: fitted cohort PTC ~5% and rebound median in clinical range (natural immunity, no boost)
    ptc, med = _ati_control(art, lambda t, kb: kb, seed=31)
    check("P6 spontaneous PTC ~5%", 100 * ptc, 1, 12, "(CHAMP/Gunst ~4-5%)")
    check("P6 rebound median in clinical range 12-30d", med, 12, 32, "(A5345 ~22d)")


def verify_p8():
    print("P8 — boosters: trough rule + reservoir x immunity synergy (reduced):")
    art = cohort_burn()
    # trough rule: same PEAK (18), two decay rates -> higher trough must give >= control
    def boosted(thalf):
        return lambda t, kb: 10.0 + 8.0 * 0.5 ** ((t - (t // 120) * 120) / thalf)
    c_lowtrough, _ = _ati_control(art, boosted(20), seed=41)   # trough ~10
    c_hitrough, _ = _ati_control(art, boosted(120), seed=42)   # trough ~14
    check("P8 trough rule: higher trough -> more control (same peak)",
          100 * (c_hitrough - c_lowtrough), 3, 100, f"(trough10 {100*c_lowtrough:.0f}% vs trough14 {100*c_hitrough:.0f}%)")
    # synergy: reservoir reduction raises cure at fixed immunity
    c_full, _ = _ati_control(art, lambda t, kb: 13.0, resv_scale=1.0, seed=51)
    c_small, _ = _ati_control(art, lambda t, kb: 13.0, resv_scale=0.1, seed=52)
    check("P8 synergy: 10x reservoir reduction raises cure (kf=13)",
          100 * (c_small - c_full), 10, 100, f"(1x {100*c_full:.0f}% vs 10x {100*c_small:.0f}%)")


def main():
    for fn in (verify_p1, verify_p3, verify_p4_p6, verify_p8):
        fn(); print()
    npass = sum(1 for c in CHECKS if c[0]); ntot = len(CHECKS)
    print(f"==== claim-chain verification: {npass}/{ntot} checks PASS ====")
    if npass < ntot:
        print("FAILED:", [c[1] for c in CHECKS if not c[0]])
    sys.exit(0 if npass == ntot else 1)


if __name__ == "__main__":
    main()
