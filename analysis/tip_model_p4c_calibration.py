#!/usr/bin/env python3
"""P4c — calibration sanity check: is the model's ATI rebound TIMING clinically realistic?

P4's conclusions are qualitative, but their credibility depends on the rebound being in the
right ballpark. Clinical/animal anchors (from the corpus): macaque ATI controls rebound at a
median ~7.5 days, bNAb-treated ~64 days (delayed); human ATI rebound is typically ~2-3 weeks;
the RIO trial shows LS-bNAbs delay/prevent rebound. So we want: (i) in the rebound regime
(weak CD8) the model rebounds within ~days-to-weeks, and (ii) stronger CD8 ('vaccinal effect')
DELAYS rebound and then controls it -- the RIO/bNAb pattern -- across the kf threshold.

This records, per replicate, the time from ATI onset until WT viral load crosses a
detection-analog threshold, starting from the model's own ART-suppressed state.
"""
import os
import numpy as np
from tip_model import T0, P
from tip_model_p4_reservoir import (simulate, pc, PSI, F_LAT, A_REACT, DL, PL, G, RDEF, S_SEED)
from tip_model_p13_wm import QS, Estar

HERE = os.path.dirname(os.path.abspath(__file__))
DETECT = 1000.0          # Vw detection-analog threshold (/mL); setpoint ~ 4e5
REPS = 120


def ati_rebound_times(art_state, kf, nu=1.0, t_ati=400.0, dt=0.04, seed=0):
    """Run ATI from an ART-suppressed state; return per-rep time-to-rebound (Vw>DETECT), nan if none."""
    rng = np.random.default_rng(seed)
    T, Iw, It, Id, Ldef, Llat, Llatd = (art_state[:, k].astype(np.float64).copy() for k in range(7))
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]
    k = QS["k"]
    N = T.shape[0]
    treb = np.full(N, np.nan)
    nsteps = int(t_ati / dt)
    for n in range(nsteps):
        Vw = pc * (Iw + (1 - rho) * Id)
        newly = np.isnan(treb) & (Vw > DETECT)
        treb[newly] = n * dt
        Vt = pc * (PSI * rho * Id)
        E = Estar(Iw + nu * Id + Ldef)
        kw, kd = kf * k * E, nu * kf * k * E

        def pois(rate):
            return rng.poisson(np.clip(rate, 0, None) * dt)
        Tprod = pois(np.full_like(T, lam)); Tdeath = pois(dT * T)
        infW = pois(b * T * Vw); infTIP = pois(b * T * Vt)
        IwSup = np.minimum(pois(b * Iw * Vt), Iw.astype(np.int64))
        ItSup = np.minimum(pois(b * It * Vw), It.astype(np.int64))
        IwD = pois((d + kw) * Iw); ItD = pois(d * It); IdD = pois((d + kd) * Id)
        Ldef_b = pois(G * Ldef + S_SEED * (Iw + Id)); Ldef_d = pois(G * Ldef * Ldef / RDEF)
        Llat_r = np.minimum(pois(A_REACT * Llat), Llat.astype(np.int64))
        Llat_d = pois(DL * Llat); Llat_p = pois(PL * Llat)
        Llatd_r = np.minimum(pois(A_REACT * Llatd), Llatd.astype(np.int64))
        Llatd_d = pois(DL * Llatd); Llatd_p = pois(PL * Llatd)
        lat = rng.binomial(np.maximum(infW.astype(np.int64), 0), F_LAT)
        T = np.clip(T + Tprod - Tdeath - infW - infTIP, 0, None)
        Iw = np.clip(Iw + (infW - lat) - IwD - IwSup + Llat_r, 0, None)
        It = np.clip(It + infTIP - ItD - ItSup, 0, None)
        Id = np.clip(Id + (IwSup + ItSup) - IdD + Llatd_r, 0, None)
        Ldef = np.clip(Ldef + Ldef_b - Ldef_d, 0, None)
        Llat = np.clip(Llat + lat - Llat_r - Llat_d + Llat_p, 0, None)
        Llatd = np.clip(Llatd + 0 - Llatd_r - Llatd_d + Llatd_p, 0, None)
    return treb


def main():
    # get the model's own ART-suppressed state (chronic -> ART, no ATI) at a reference kf
    s0 = np.zeros((REPS, 7)); s0[:, 0] = T0; s0[:, 1] = 10
    art_state = simulate(s0, 1.0, np.full(REPS, 12.0), 500.0, 300.0, 0.0, seed=7)
    print(f"ART-suppressed start: Iw med={np.median(art_state[:,1]):.0f}, "
          f"latent reservoir L_lat med={np.median(art_state[:,5]):.0f}\n")
    print(f"time-to-rebound (Vw>{DETECT:.0f}/mL) from ATI onset, {REPS} reps:")
    print(f"clinical anchors: macaque controls ~7.5d, human ATI ~2-3 wk, bNAb-treated ~64d (RIO)\n")
    rows = []
    for kf in [6.0, 8.0, 10.0, 12.0, 14.0]:
        tr = ati_rebound_times(art_state, kf, seed=int(kf))
        reb = tr[~np.isnan(tr)]
        frac = len(reb) / len(tr)
        med = np.median(reb) if len(reb) else float("nan")
        iqr = (np.percentile(reb, 25), np.percentile(reb, 75)) if len(reb) else (np.nan, np.nan)
        print(f"  kf={kf:4.1f}: rebound in {frac:3.0%} of reps; median t_rebound="
              f"{med:6.1f} d  (IQR {iqr[0]:.0f}-{iqr[1]:.0f})")
        rows.append((kf, frac, med, iqr[0], iqr[1]))
    np.savez(os.path.join(HERE, "p4c_calibration.npz"), rows=np.array(rows), detect=DETECT)
    print("\nReading: low kf (weak CD8) = fast rebound (~clinical days-to-weeks); raising kf "
          "DELAYS then abolishes rebound -- the bNAb/RIO 'delayed-or-controlled' pattern. "
          "The reactivation rate sets the rebound clock and is the first knob to calibrate to a "
          "specific ATI cohort for quantitative claims.")
    print("wrote p4c_calibration.npz")


if __name__ == "__main__":
    main()
