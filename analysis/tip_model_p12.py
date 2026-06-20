#!/usr/bin/env python3
"""P1.2 — TIP vs *realistic* CD8 immunity (memory + exhaustion).

P1.1 used the simplest effector ODE and found (a) predator-prey oscillations and
(b) that early TIP suppresses CD8 priming -> TIPs may ANTAGONIZE the vaccinal-
effect control the field is betting on. Two caveats drove P1.2:
  - the oscillations are likely an artifact of an effector model with no memory;
  - real chronic-HIV CD8 are EXHAUSTED, which makes the TIP two-sided: lowering
    antigen could *relieve* exhaustion and strengthen control.

Immune model now has memory M (persistence -> damps oscillation) and functional
exhaustion (killing per effector wanes at high antigen):
  A    = Iw + Id                      # antigen = productive infected cells
  stim = A/(A+K_E)
  kill = k*E / (1 + A/Kx)             # EXHAUSTION: high antigen -> weak killing
  dE/dt = a*stim*(E+omega*M)*(1-(E+M)/Emax) - dE*E - chi*(A/(A+Kx))*E
  dM/dt = rM*stim*E - dM*M            # memory (slow decay -> persistence)
productive cells Iw, Id die at (d + kill).

Tests:
  A. does memory damp the P1.1 oscillation into a stable controlled setpoint?
  B. ANTAGONISM SEARCH: over (administration time t_admin x TIP advantage psi),
     time-averaged WT reduction AND retained CD8 (final E+M vs no-TIP). Is there an
     "immune-compatible" TIP (lowers WT *and* preserves CD8)? Or are they exclusive?

Metric: time-AVERAGED Vw over the last 300 d (handles cycles, not t=end snapshots).
Reuses HIV params from tip_model.py.
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tip_model import P, T0, _san

# tuned for a chronic CONTROLLED setpoint (~0.7 log below the no-immunity 7.7e5)
# with persistent, exhausted CD8. NOTE: antigen-driven CD8 control is intrinsically
# oscillatory in these ODEs even WITH memory+exhaustion (a robust feature, not a
# minimal-model artifact) -> we answer the antagonism question on cycle-TIME-AVERAGES.
# mem=0 (memory-contributed killing drove boom-bust relaxation cycles; disabled).
IM = dict(a=1.2, K_E=3000.0, dE=0.25, k=2e-4, rM=0.05, dM=0.005,
          omega=0.5, Emax=3e4, chi=0.30, Kx=2500.0, mem=0.0)


def rhs(t, y, psi):
    T, Iw, It, Id, Vw, Vt, E, M = np.clip(y, 0.0, 1e13)
    lam, dT, d, c, p, b, bt, rho = (P[k] for k in
        ("lam", "dT", "d", "c", "p", "b", "bt", "rho"))
    a, K_E, dE, k, rM, dM, om, Emax, chi, Kx, mem = (IM[k] for k in
        ("a", "K_E", "dE", "k", "rM", "dM", "omega", "Emax", "chi", "Kx", "mem"))
    A = Iw + Id
    stim = A/(A+K_E)
    kill = k*(E + mem*M)/(1.0 + A/Kx)        # memory contributes (fast-recall proxy) -> damps cycle
    return [
        lam - dT*T - b*T*Vw - bt*T*Vt,
        b*T*Vw - (d+kill)*Iw - bt*Iw*Vt,
        bt*T*Vt - d*It - b*It*Vw,
        bt*Iw*Vt + b*It*Vw - (d+kill)*Id,
        p*Iw + (1-rho)*p*Id - c*Vw,
        psi*rho*p*Id - c*Vt,
        a*stim*(E+om*M)*(1-(E+M)/Emax) - dE*E - chi*(A/(A+Kx))*E,
        rM*stim*E - dM*M,
    ]


def integ(y0, psi, t0, t1, teval=None):
    s = solve_ivp(rhs, (t0, t1), _san(y0), args=(psi,), method="LSODA",
                  rtol=1e-6, atol=1e-2, max_step=2.0, t_eval=teval)
    return s


ACUTE = [T0, 0, 0, 0, 1e-3, 0, 1.0, 0.0]


def tail_avg(y0, psi, t0, run=1500, win=500):   # long window -> average over many cycles
    s = integ(y0, psi, t0, t0+run, teval=np.linspace(t0+run-win, t0+run, 300))
    Vw = np.clip(np.nan_to_num(s.y[4], posinf=1e13), 1e-3, None)
    EM = np.clip(np.nan_to_num(s.y[6]+s.y[7], posinf=1e13), 0, None)
    return np.exp(np.mean(np.log(Vw))), np.mean(EM), s   # geometric-mean Vw, mean CD8


def run_tip(psi, t_admin, dose=1e2):
    if t_admin > 0:
        s1 = integ(ACUTE, psi, 0, t_admin)
        st = _san(s1.y[:, -1]); st[5] = dose
    else:
        st = list(ACUTE); st[5] = dose
    return tail_avg(st, psi, t_admin)


def main():
    quick = "quick" in sys.argv
    # ---- baseline: WT + immunity, no TIP ----
    vw0, em0, sb = tail_avg(ACUTE, 1.0, 0.0)
    # detect oscillation: spread of Vw over the averaging window
    Vwwin = np.clip(sb.y[4], 1e-3, None)
    osc = (np.max(Vwwin)/max(np.min(Vwwin), 1e-3))
    print(f"BASELINE (no TIP): time-avg Vw={vw0:.3e}, CD8 E+M={em0:.3e}, "
          f"window max/min Vw={osc:.2f} ({'STABLE' if osc < 1.3 else 'OSCILLATING'})")
    sfull = integ(ACUTE, 1.0, 0, 1300, teval=np.linspace(0, 1300, 700))
    print(f"  setpoint check: Vw(1300)={sfull.y[4,-1]:.3e}, Iw={sfull.y[1,-1]:.2e}, "
          f"E={sfull.y[6,-1]:.2e}, M={sfull.y[7,-1]:.2e}, "
          f"kill={IM['k']*sfull.y[6,-1]/(1+(sfull.y[1,-1]+sfull.y[3,-1])/IM['Kx']):.3f}/day")

    # ---- Exp A: timecourse, early vs late TIP under realistic immunity ----
    se = integ([T0, 0, 0, 0, 1e-3, 1e2, 1.0, 0.0], 15.0, 0, 1300, teval=np.linspace(0, 1300, 700))
    sb2 = integ(ACUTE, 15.0, 0, 400)             # baseline to chronic, then late TIP
    stl = _san(sb2.y[:, -1]); stl[5] = 1e2
    sl = integ(stl, 15.0, 400, 1300, teval=np.linspace(400, 1300, 500))
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ax[0].semilogy(sfull.t, np.clip(sfull.y[4], 1e-3, None), "k:", lw=2, label="no TIP")
    ax[0].semilogy(se.t, np.clip(se.y[4], 1e-3, None), "tab:green", lw=2, label="early TIP (t=0)")
    ax[0].semilogy(sl.t, np.clip(sl.y[4], 1e-3, None), "tab:orange", lw=2, label="late TIP (t=400)")
    ax[0].set_xlabel("days"); ax[0].set_ylabel("WT load Vw"); ax[0].legend(fontsize=9)
    ax[0].set_title("Exp A: WT under realistic CD8 (memory+exhaustion)\nmemory should damp the P1.1 oscillation")
    ax[1].semilogy(sfull.t, np.clip(sfull.y[6]+sfull.y[7], 1e-2, None), "k:", lw=2, label="no TIP")
    ax[1].semilogy(se.t, np.clip(se.y[6]+se.y[7], 1e-2, None), "tab:green", lw=2, label="early TIP")
    ax[1].semilogy(sl.t, np.clip(sl.y[6]+sl.y[7], 1e-2, None), "tab:orange", lw=2, label="late TIP")
    ax[1].set_xlabel("days"); ax[1].set_ylabel("CD8 E+M"); ax[1].legend(fontsize=9)
    ax[1].set_title("Exp A: does the TIP erase the CD8 response?")
    fig.tight_layout(); fig.savefig("/home/seth/dev/hiv-aids-research/analysis/p12_timecourse.png", dpi=130)
    print(f"\nExp A (psi=15): early-TIP avg CD8={np.mean(se.y[6,-100:]+se.y[7,-100:]):.2e}, "
          f"late-TIP CD8={np.mean(sl.y[6,-100:]+sl.y[7,-100:]):.2e}, no-TIP CD8={em0:.2e}")
    if quick:
        print("quick mode: skipping Exp B sweep"); return

    # ---- Exp B: antagonism search over (t_admin, psi) ----
    TADM = np.linspace(0, 250, 14)
    PSI = np.linspace(1.0, 20.0, 14)
    wtred = np.zeros((len(PSI), len(TADM)))
    cd8ret = np.zeros((len(PSI), len(TADM)))
    for j, ta in enumerate(TADM):
        for i, psi in enumerate(PSI):
            vw, em, _ = run_tip(psi, ta)
            wtred[i, j] = np.log10(vw0) - np.log10(vw)
            cd8ret[i, j] = em / max(em0, 1e-9)
    np.savez("/home/seth/dev/hiv-aids-research/analysis/p12_sweep.npz",
             TADM=TADM, PSI=PSI, wtred=wtred, cd8ret=cd8ret, vw0=vw0, em0=em0)

    # immune-compatible = WT down >=1 log AND CD8 retained >=50%
    compat = (wtred >= 1.0) & (cd8ret >= 0.5)
    fig2, ax2 = plt.subplots(1, 2, figsize=(13, 5))
    im0 = ax2[0].pcolormesh(TADM, PSI, wtred, shading="auto", cmap="viridis")
    ax2[0].contour(TADM, PSI, cd8ret, levels=[0.5], colors="white", linewidths=2)
    if compat.any():
        ax2[0].contourf(TADM, PSI, compat.astype(float), levels=[0.5, 1.5],
                        colors="none", hatches=["xx"])
    ax2[0].set_xlabel("TIP administration time (days post-infection)")
    ax2[0].set_ylabel("TIP advantage psi")
    ax2[0].set_title("WT log-reduction (color); CD8>=50% retained (white line);\nimmune-compatible = hatched")
    fig2.colorbar(im0, ax=ax2[0], label="time-avg WT log10 reduction")
    im1 = ax2[1].pcolormesh(TADM, PSI, cd8ret, shading="auto", cmap="RdYlGn", vmin=0, vmax=1.2)
    ax2[1].set_xlabel("TIP administration time (days post-infection)")
    ax2[1].set_ylabel("TIP advantage psi")
    ax2[1].set_title("CD8 retained vs no-TIP\n(<1 = TIP suppressed the response)")
    fig2.colorbar(im1, ax=ax2[1], label="(E+M with TIP) / (E+M no TIP)")
    fig2.tight_layout(); fig2.savefig("/home/seth/dev/hiv-aids-research/analysis/p12_antagonism.png", dpi=130)

    print("\n--- Exp B: antagonism search ---")
    print(f"cells with WT down >=1 log: {(wtred>=1).sum()}/{wtred.size}")
    print(f"cells with CD8 retained >=50%: {(cd8ret>=0.5).sum()}/{cd8ret.size}")
    print(f"IMMUNE-COMPATIBLE (both): {compat.sum()}/{compat.size} "
          f"({'EXISTS' if compat.any() else 'EMPTY -> antagonism is a hard constraint'})")
    if compat.any():
        ii, jj = np.where(compat)
        best = np.argmax(wtred[ii, jj])
        print(f"  best immune-compatible: t_admin={TADM[jj[best]]:.0f}d, psi={PSI[ii[best]]:.1f}, "
              f"WT down {wtred[ii[best],jj[best]]:.2f} log, CD8 {cd8ret[ii[best],jj[best]]:.0%} retained")
    # the antagonism gradient: do high-WT-reduction cells cost CD8?
    strong = wtred >= 1.5
    if strong.any():
        print(f"  among strong-suppression cells (WT down>=1.5 log): "
              f"median CD8 retained = {np.median(cd8ret[strong]):.0%}")
    print("wrote p12_timecourse.png, p12_antagonism.png, p12_sweep.npz")


if __name__ == "__main__":
    main()
