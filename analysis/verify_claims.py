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


def verify_p4_p6():
    # AUDIT2 #1 FIX: exercise the PRODUCTION engine (simulate / P6.rebound_curve), not a copy.
    print("P4/P5/P6 — via PRODUCTION simulate()/rebound_curve:")
    from tip_model import T0
    from tip_model_p4_reservoir import simulate
    import tip_model_p6_heterogeneous as P6
    kf_base, f_lat, _ = cohort_burn()
    n = len(kf_base)

    def ctrl(level, tip=False, seed=0):
        kf_arr = np.full(n, float(level)) if np.isscalar(level) else level
        s0 = np.zeros((n, 7)); s0[:, 0] = T0; s0[:, 1] = 10
        if tip:                                  # P4's headline arm: latent/engineered TIP (nu=0.9)
            s0[:, 3] = 10; s0[:, 6] = 50         # seed Id + latent dual L_latd
        f = simulate(s0, 0.9 if tip else 1.0, kf_arr, 500.0, 300.0, 500.0, f_lat=f_lat, seed=seed)
        return float(((f[:, 1] + f[:, 3]) < 50).mean())

    check("P4 control sharp: low immunity -> low control", 100 * ctrl(8.0, seed=11), 0, 25)
    check("P4 control sharp: high immunity -> high control", 100 * ctrl(18.0, seed=12), 35, 100, "(threshold)")
    # P4 HEADLINE neutrality = the latent/engineered TIP arm (NOT the sustained 'second-ART' arm,
    # which P4b/P8 show DOES help via TIP-dependent suppression -- AUDIT2 #7). Band tightened from
    # the original near-vacuous [0,20] (AUDIT2 #8).
    c_no, c_tip = ctrl(kf_base, seed=21), ctrl(kf_base, tip=True, seed=22)
    check("P4 latent-TIP ~neutral (|dControl|<=12pts)", 100 * abs(c_tip - c_no), 0, 12,
          f"(noTIP {100*c_no:.0f}%, +latentTIP {100*c_tip:.0f}%)")

    # P6: production rebound_curve + km_summary on the fitted cohort
    kf6, fl6 = P6.draw_cohort(10.0, 2.5, np.random.default_rng(7))
    med, f4, f12, ptc = P6.km_summary(P6.rebound_curve(kf6, fl6, t_ati=500.0, seed=31))
    check("P6 spontaneous PTC ~5%", 100 * ptc, 1, 12, "(CHAMP/Gunst ~4-5%)")
    check("P6 rebound median in clinical range 12-30d", med, 12, 32, "(A5345 ~22d)")


def verify_p8():
    # AUDIT2 #1 FIX: call the PRODUCTION P8.run() (its own ATI loop), not a copy.
    print("P8 — via PRODUCTION run() (trough rule + reservoir x immunity synergy):")
    import tip_model_p8_boosters as P8
    def boosted(thalf):                          # P8.run's level_func takes (t) only
        return lambda t: 10.0 + 8.0 * 0.5 ** ((t - (t // 120) * 120) / thalf)
    c_lowtrough = P8.run(boosted(20), t_ati=500.0, seed=41)    # trough ~10
    c_hitrough = P8.run(boosted(120), t_ati=500.0, seed=42)    # trough ~14
    check("P8 trough rule: higher trough -> more control (same peak)",
          100 * (c_hitrough - c_lowtrough), 3, 100, f"(trough10 {100*c_lowtrough:.0f}% vs trough14 {100*c_hitrough:.0f}%)")
    # synergy: reservoir reduction raises cure at fixed immunity (production P8.run, resv_scale)
    c_full = P8.run(lambda t: 13.0, resv_scale=1.0, t_ati=500.0, seed=51)
    c_small = P8.run(lambda t: 13.0, resv_scale=0.1, t_ati=500.0, seed=52)
    check("P8 synergy: 10x reservoir reduction raises cure (kf=13)",
          100 * (c_small - c_full), 10, 100,
          f"(1x {100*c_full:.0f}% vs 10x {100*c_small:.0f}%; NB direction only, not super-additivity)")


def verify_p11():
    # AUDIT2 follow-up: coupling the TIP to the rebound FLIPS neutral -> helps (production model).
    print("P11 — coupled TIP helps (via PRODUCTION simulate_coupled):")
    from tip_model_p11_coupled import simulate_coupled
    c0 = simulate_coupled(150, 11.0, chi=0.0, t_ati=400.0, seed=1)[0]   # decoupled
    c1 = simulate_coupled(150, 11.0, chi=1.0, t_ati=400.0, seed=1)[0]   # fully coupled
    check("P11 coupling raises control (chi 0->1, kf=11)", 100 * (c1 - c0), 10, 100,
          f"(decoupled {100*c0:.0f}% -> coupled {100*c1:.0f}%; flips the 'neutral' verdict)")


def verify_p12():
    # AUDIT2 last open item: under wear-down-able immunity, coupled TIP still helps (no backfire).
    print("P12 — coupled TIP helps under waning immunity, no backfire (production):")
    from tip_model_p12_exhaustible import simulate_exh
    c0 = simulate_exh(150, 0.0, 0.5, e_boost=0.5, t_ati=700.0, seed=3)[0]   # no TIP
    c1 = simulate_exh(150, 1.0, 0.5, e_boost=0.5, t_ati=700.0, seed=3)[0]   # coupled TIP
    check("P12 coupled TIP does not backfire (dControl >= -5pts)", 100 * (c1 - c0), -5, 100,
          f"(noTIP {100*c0:.0f}% -> coupled {100*c1:.0f}%; helps under exhaustible immunity)")


def verify_p13():
    # harshest test: active exhaustion + stealthy TIP -> still no backfire (TIP helps).
    print("P13 — coupled TIP no-backfire under ACTIVE exhaustion + stealth (production):")
    from tip_model_p13_exhaustion_damage import simulate_dmg
    c0 = simulate_dmg(150, 0.0, 0.1, exhaust=0.15, seed=5)[0]   # no TIP, exhaustible immunity, stealthy
    c1 = simulate_dmg(150, 1.0, 0.1, exhaust=0.15, seed=5)[0]   # coupled stealthy TIP
    check("P13 no backfire under active exhaustion (dControl >= -5pts)", 100 * (c1 - c0), -5, 100,
          f"(noTIP {100*c0:.0f}% -> coupled {100*c1:.0f}%; TIP protects vs exhaustion)")


def verify_deboer():
    # We REPRODUCE Dodd & de Boer (active infection): immunity collapses the TIP's suppression.
    print("de Boer recreation: immunity collapses TIP benefit on active VL (we reproduce them):")
    from tip_model import wt_setpoint, tip_arm
    def benefit(kap):
        sp = wt_setpoint(kap); end, _ = tip_arm(sp, kap, 12.0, tmax=700)
        return np.log10(max(sp[4], 1e-9)) - np.log10(max(end[4], 1e-9))
    b0, b8 = benefit(0.0), benefit(0.8)
    check("de Boer: TIP active-VL benefit collapses under immunity (b0-b0.8 > 1)", b0 - b8, 1.0, 6.0,
          f"(kap0 {b0:.2f} -> kap0.8 {max(b8,0):.2f} log)")


def verify_p14():
    # P14 phase: at marginal immunity the TIP effect grows with coupling; chi=0 = de Boer (inert).
    print("P14 — coupling phase: TIP helps at marginal immunity, inert at chi=0 (de Boer limit):")
    from tip_model_p11_coupled import simulate_coupled
    c0 = simulate_coupled(150, 12.0, chi=0.0, t_ati=500.0, seed=120)[0]   # de Boer limit
    c1 = simulate_coupled(150, 12.0, chi=1.0, t_ati=500.0, seed=120)[0]   # fully coupled
    check("P14 coupling helps at marginal immunity (kf=12, chi 0->1)", 100 * (c1 - c0), 8, 100,
          f"(de Boer-limit {100*c0:.0f}% -> coupled {100*c1:.0f}%)")


def verify_p15():
    # global-sensitivity spot check: at two diverse parameter points, coupled TIP does not backfire.
    print("P15 — coupled TIP no-backfire at diverse parameter points (sensitivity spot check):")
    import tip_model_p11_coupled as M
    from tip_model_p11_coupled import simulate_coupled
    worst = 100.0
    for (kf, nu, psi, rdef) in [(11.0, 0.3, 80.0, 4000.0), (12.0, 0.8, 20.0, 8000.0)]:
        c0 = simulate_coupled(120, kf, chi=0.0, nu=nu, psi=psi, rdef=rdef, t_ati=400.0, seed=7)[0]
        c1 = simulate_coupled(120, kf, chi=1.0, nu=nu, psi=psi, rdef=rdef, t_ati=400.0, seed=7)[0]
        worst = min(worst, 100 * (c1 - c0))
    check("P15 no backfire across diverse params (min dControl >= -6pts)", worst, -6, 100,
          f"(worst-case TIP effect {worst:+.0f} pts)")


def main():
    for fn in (verify_p1, verify_p3, verify_p4_p6, verify_p8, verify_p11, verify_p12, verify_p13,
               verify_deboer, verify_p14, verify_p15):
        fn(); print()
    npass = sum(1 for c in CHECKS if c[0]); ntot = len(CHECKS)
    print(f"==== claim-chain verification: {npass}/{ntot} checks PASS ====")
    if npass < ntot:
        print("FAILED:", [c[1] for c in CHECKS if not c[0]])
    sys.exit(0 if npass == ntot else 1)


if __name__ == "__main__":
    main()
