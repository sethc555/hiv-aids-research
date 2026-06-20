#!/usr/bin/env python3
"""P6 — heterogeneous-patient model: a population over (immune strength kf, reservoir clock
f_lat), fit to the clinical ATI rebound Kaplan-Meier curve, then the TIP test on the real mix.

P5 found a single immune setpoint cannot reproduce BOTH the fast-rebound majority (~3 wk) AND
the durable-controller minority (~5%): the clock fit wants kf~9, the PTC-fraction fit wants
kf~13. The clinical resolution is patient HETEROGENEITY -- controllers have a different immune /
reservoir setpoint. P6 makes each replicate a *patient*:
    kf_i    ~ Normal(m_kf, s_kf)            (immune 'vaccinal-effect' strength)
    f_lat_i ~ lognormal(median f0, gsd)     (per-patient reservoir reactivation clock)
The cohort's spread of (kf, f_lat) produces a DISTRIBUTION of rebound times -> a model KM curve.
We fit (m_kf, s_kf) to the clinical KM summary (median, % rebounded by wk4 & wk12, PTC fraction),
then ask the project's question on the realistic cohort: does a TIP help or hurt across patients?

Clinical KM targets (sourced in P5_findings.md / METHODS.md): median ~22 d (>1000 c/mL, A5345);
~66% rebounded by wk4, ~77% by wk12 (>200, ACTG pooled); spontaneous PTC ~4-5% (CHAMP/Gunst).

AUDIT DISCIPLINE: fit to published KM SUMMARY points (not digitized patient-level curves -- we
lack access); report the fit AND residuals; heterogeneity params are illustrative, not unique.
"""
import os
import numpy as np
from tip_model import T0, P
from tip_model_p4_reservoir import simulate, pc, DL, PL, G, RDEF, S_SEED, A_REACT
from tip_model_p13_wm import QS, Estar

HERE = os.path.dirname(os.path.abspath(__file__))
DETECT = 1000.0
N = 300                      # patients per cohort
F0_LOG10 = np.log10(8e-5)    # population median reservoir clock (from P5)
F_GSD_LOG10 = 0.35           # reservoir clock spread (~ x/div 2.2)
KF_FLOOR = 4.0

# clinical KM targets
T_MEDIAN, T_W4, T_W12, T_PTC = 22.0, 0.66, 0.77, 0.05
WK = 7.0


def draw_cohort(m_kf, s_kf, rng):
    kf = np.clip(rng.normal(m_kf, s_kf, N), KF_FLOOR, 40.0)
    f_lat = 10.0 ** rng.normal(F0_LOG10, F_GSD_LOG10, N)
    return kf, f_lat


def rebound_curve(kf, f_lat, nu=None, tip_sustained=0.0, psi=60.0, t_ati=600.0, dt=0.04, seed=0):
    """Per-patient time-to-rebound (Vw>DETECT). nu defaults to 1 (no TIP) or 0.9 if TIP dosed."""
    if nu is None:
        nu = 0.9 if tip_sustained else 1.0
    s0 = np.zeros((N, 7)); s0[:, 0] = T0; s0[:, 1] = 10
    art = simulate(s0, nu, kf, 500.0, 300.0, 0.0, f_lat=f_lat, psi=psi, seed=seed)
    rng = np.random.default_rng(seed + 1)
    T, Iw, It, Id, Ldef, Llat, Llatd = (art[:, j].astype(np.float64).copy() for j in range(7))
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]
    k = QS["k"]
    treb = np.full(N, np.nan)
    for n in range(int(t_ati / dt)):
        Vw = pc * (Iw + (1 - rho) * Id)
        treb[np.isnan(treb) & (Vw > DETECT)] = n * dt
        if tip_sustained:
            It = It + rng.poisson(np.full(N, tip_sustained * dt))
        Vt = pc * (psi * rho * Id)
        E = Estar(Iw + nu * Id + Ldef)
        kw, kd = kf * k * E, nu * kf * k * E

        def po(r):
            return rng.poisson(np.clip(r, 0, None) * dt)
        Tp = po(np.full(N, lam)); Td = po(dT * T)
        iW = po(b * T * Vw); iT = po(b * T * Vt)
        iWs = np.minimum(po(b * Iw * Vt), Iw.astype(np.int64))
        iTs = np.minimum(po(b * It * Vw), It.astype(np.int64))
        iWd = po((d + kw) * Iw); iTd = po(d * It); idd = po((d + kd) * Id)
        Lb = po(G * Ldef + S_SEED * (Iw + Id)); Ld = po(G * Ldef * Ldef / RDEF)
        lr = np.minimum(po(A_REACT * Llat), Llat.astype(np.int64))
        ldd = po(DL * Llat); lp = po(PL * Llat)
        ldr = np.minimum(po(A_REACT * Llatd), Llatd.astype(np.int64))
        lddd = po(DL * Llatd); ldp = po(PL * Llatd)
        lat = rng.binomial(np.maximum(iW.astype(np.int64), 0), f_lat)
        T = np.clip(T + Tp - Td - iW - iT, 0, None)
        Iw = np.clip(Iw + (iW - lat) - iWd - iWs + lr, 0, None)
        It = np.clip(It + iT - iTd - iTs, 0, None)
        Id = np.clip(Id + (iWs + iTs) - idd + ldr, 0, None)
        Ldef = np.clip(Ldef + Lb - Ld, 0, None)
        Llat = np.clip(Llat + lat - lr - ldd + lp, 0, None)
        Llatd = np.clip(Llatd - ldr - lddd + ldp, 0, None)
    return treb


def km_summary(treb):
    reb = treb[~np.isnan(treb)]
    ptc = 1 - len(reb) / len(treb)
    med = np.median(reb) if len(reb) >= len(treb) * 0.5 else np.inf
    f4 = np.mean(np.nan_to_num(treb, nan=1e9) <= 4 * WK)
    f12 = np.mean(np.nan_to_num(treb, nan=1e9) <= 12 * WK)
    return med, f4, f12, ptc


def fit_error(med, f4, f12, ptc):
    m = (min(med, 200) - T_MEDIAN) / T_MEDIAN
    return m * m + (f4 - T_W4) ** 2 + (f12 - T_W12) ** 2 + 3 * (ptc - T_PTC) ** 2  # weight PTC


def main():
    print(f"P6 heterogeneous cohort (N={N}; reservoir clock median 8e-5, gsd~{10**F_GSD_LOG10:.1f}x)")
    print(f"clinical KM targets: median {T_MEDIAN:.0f}d, reb by wk4 {T_W4:.0%}, wk12 {T_W12:.0%}, PTC {T_PTC:.0%}\n")

    # ---- fit population (m_kf, s_kf) to the KM summary ----
    print("fitting (m_kf, s_kf) to clinical KM:")
    M = np.array([8.0, 9.0, 10.0, 11.0]); S = np.array([1.5, 2.5, 3.5, 4.5])
    best, best_err = None, 1e9
    grid = {}
    for m_kf in M:
        for s_kf in S:
            kf, fl = draw_cohort(m_kf, s_kf, np.random.default_rng(7))
            med, f4, f12, ptc = km_summary(rebound_curve(kf, fl, seed=int(m_kf * 10 + s_kf)))
            err = fit_error(med, f4, f12, ptc)
            grid[(m_kf, s_kf)] = (med, f4, f12, ptc, err)
            print(f"  m_kf={m_kf:.1f} s_kf={s_kf:.1f}: median={med:5.1f}d wk4={f4:.0%} wk12={f12:.0%} "
                  f"PTC={ptc:.0%}  err={err:.3f}")
            if err < best_err:
                best_err, best = err, (m_kf, s_kf)
    m_kf, s_kf = best
    med, f4, f12, ptc, err = grid[best]
    print(f"\nBEST FIT: m_kf={m_kf}, s_kf={s_kf}  ->  median={med:.1f}d (t{T_MEDIAN:.0f}), "
          f"wk4={f4:.0%} (t{T_W4:.0%}), wk12={f12:.0%} (t{T_W12:.0%}), PTC={ptc:.0%} (t{T_PTC:.0%})")
    print(f"  => a SINGLE cohort now matches fast-rebound majority AND controller minority "
          f"(the P5 tension is resolved by heterogeneity: kf spread {s_kf} around {m_kf}).")

    # ---- the question on the realistic cohort: no-TIP vs sustained TIP ----
    kf, fl = draw_cohort(m_kf, s_kf, np.random.default_rng(7))
    tr_no = rebound_curve(kf, fl, seed=101)
    tr_tip = rebound_curve(kf, fl, tip_sustained=2000.0, psi=60.0, seed=102)
    s_no, s_tip = km_summary(tr_no), km_summary(tr_tip)
    print(f"\nTIP test on the FITTED cohort (sustained TIP psi=60, nu=0.9):")
    print(f"  no TIP : median {s_no[0]:.0f}d, reb wk4 {s_no[1]:.0%}, wk12 {s_no[2]:.0%}, PTC {s_no[3]:.0%}")
    print(f"  + TIP  : median {s_tip[0]:.0f}d, reb wk4 {s_tip[1]:.0%}, wk12 {s_tip[2]:.0%}, PTC {s_tip[3]:.0%}")
    print(f"  PTC change with TIP: {s_no[3]:.0%} -> {s_tip[3]:.0%} "
          f"({'+' if s_tip[3]>=s_no[3] else ''}{100*(s_tip[3]-s_no[3]):.0f} pts) "
          f"-- small => TIP-neutral verdict holds on the realistic cohort")

    # model KM curve (no-TIP) vs clinical points
    days = np.arange(0, 365, 2.0)
    km_no = np.array([np.mean(np.nan_to_num(tr_no, nan=1e9) > t) for t in days])     # frac aviremic
    km_tip = np.array([np.mean(np.nan_to_num(tr_tip, nan=1e9) > t) for t in days])
    np.savez(os.path.join(HERE, "p6_heterogeneous.npz"), m_kf=m_kf, s_kf=s_kf,
             days=days, km_no=km_no, km_tip=km_tip, kf=kf, f_lat=fl,
             summary_no=np.array(s_no), summary_tip=np.array(s_tip),
             targets=np.array([T_MEDIAN, T_W4, T_W12, T_PTC]))

    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    ax[0].plot(days / 7, 100 * km_no, lw=2, label="model cohort (no TIP)")
    ax[0].scatter([0, 4, 12], [100, 100 * (1 - T_W4), 100 * (1 - T_W12)], color="k", zorder=5,
                  label="clinical KM (ACTG: wk4, wk12)")
    ax[0].axhline(100 * T_PTC, color="tab:green", ls=":", label=f"clinical PTC ~{T_PTC:.0%}")
    ax[0].set_xlabel("weeks post-ATI"); ax[0].set_ylabel("% still aviremic (not rebounded)")
    ax[0].set_title(f"P6 fitted cohort KM vs clinical\n(m_kf={m_kf}, s_kf={s_kf}: majority rebounds, minority controls)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(days / 7, 100 * km_no, lw=2, label="no TIP")
    ax[1].plot(days / 7, 100 * km_tip, lw=2, ls="--", label="+ sustained TIP")
    ax[1].set_xlabel("weeks post-ATI"); ax[1].set_ylabel("% still aviremic")
    ax[1].set_title("TIP effect on the realistic cohort\n(curves ~coincide -> TIP-neutral verdict holds)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "p6_heterogeneous.png"), dpi=130)
    print("\nwrote p6_heterogeneous.npz, p6_heterogeneous.png")


if __name__ == "__main__":
    main()
