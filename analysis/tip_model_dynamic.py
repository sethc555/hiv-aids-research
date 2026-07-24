#!/usr/bin/env python3
"""P1.1 — TIP vs *dynamic* immunity (the novel step).

P1 held immunity at fixed killing pressure kap and found TIP efficacy collapses
once kap is moderate. But a working TIP lowers viral load -> lowers antigen ->
contracts the immune response -> less killing -> the TIP's dual-infected
production cells survive. P1 omits exactly this feedback.

Immunity is a dynamic effector E(t):  dE/dt = a*E*A/(A+K_E) - dE_*E,  A = Iw+Id;
productive cells Iw, Id die at (d + k*E). 'a' = immune responsiveness.

NOTE: at high responsiveness (a >~ 0.7) the effector-virus system goes into
predator-prey LIMIT CYCLES, so a fixed "baseline" is undefined there. We restrict
the clean analysis to the STABLE moderate-immunity regime (a in [0.25, 0.6]) —
exactly where de Boer's "moderate immune response" claim lives — and flag the
oscillatory regime as out of scope.

Two tests:
  (1) RESCUE map: dynamic TIP benefit minus the P1 static-kap prediction at the
      matched baseline killing pressure. Does the feedback open dead zones?
  (2) TIMING: TIP co-administered at ACUTE infection (before immunity matures)
      vs added at the CHRONIC setpoint (established immunity). Reconciles
      de Boer (chronic, pessimistic) with Pitchai 2024 (early challenge, worked).

Reuses the P1 static model (tip_model.py).
"""
import os
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tip_model import P, T0, wt_setpoint as wt_static, tip_arm as tip_static, _san

HERE = os.path.dirname(os.path.abspath(__file__))   # path-portable (AUDIT2 #14)

IK = dict(k=2e-4, K_E=5000.0, dE=0.2)
E_INIT = 1.0


def rhs_dyn(t, y, a, psi):
    T, Iw, It, Id, Vw, Vt, E = np.clip(y, 0.0, 1e13)
    lam, dT, d, c, p, b, bt, rho = (P[k] for k in
        ("lam", "dT", "d", "c", "p", "b", "bt", "rho"))
    A = Iw + Id
    kE = IK["k"] * E
    return [
        lam - dT*T - b*T*Vw - bt*T*Vt,
        b*T*Vw - (d+kE)*Iw - bt*Iw*Vt,
        bt*T*Vt - d*It - b*It*Vw,
        bt*Iw*Vt + b*It*Vw - (d+kE)*Id,
        p*Iw + (1-rho)*p*Id - c*Vw,
        psi*rho*p*Id - c*Vt,
        a*E*A/(A+IK["K_E"]) - IK["dE"]*E,
    ]


def integ(y0, a, psi, tmax=900, teval=None):
    y0 = _san(y0)
    s = solve_ivp(rhs_dyn, (0, tmax), y0, args=(a, psi), method="LSODA",
                  rtol=1e-6, atol=1e-2, max_step=2.0, t_eval=teval)
    return _san(s.y[:, -1]), s


def base_dyn(a):                       # WT + immunity, no TIP -> chronic setpoint
    end, _ = integ([T0, 0, 0, 0, 1e-3, 0, E_INIT], a, 1.0)
    return end


def tip_late(state, a, psi, dose=1e2):  # add TIP at established setpoint
    y0 = list(state); y0[5] = dose
    return integ(y0, a, psi)


def tip_early(a, psi, dose=1e2):        # co-inoculate TIP with WT at acute infection
    return integ([T0, 0, 0, 0, 1e-3, dose, E_INIT], a, psi)


def main():
    print("a     Vw0        Iw+Id      kap_eff0=k*E0   (stability check)")
    for a in [0.25, 0.35, 0.45, 0.55, 0.6]:
        sp = base_dyn(a)
        Apred = IK["K_E"]*IK["dE"]/max(a-IK["dE"], 1e-9)
        print(f"{a:<5} {sp[4]:.3e}  {sp[1]+sp[3]:.3e}  {IK['k']*sp[6]:.3f}        "
              f"A*pred={Apred:.0f} vs A={sp[1]+sp[3]:.0f}")

    A = np.linspace(0.25, 0.6, 18)
    PSI = np.linspace(1.0, 20.0, 18)
    benefit = np.zeros((len(PSI), len(A)))
    persist = np.zeros((len(PSI), len(A)))
    rescue = np.zeros((len(PSI), len(A)))
    keff0 = np.zeros(len(A))
    for j, a in enumerate(A):
        sp = base_dyn(a); vw0 = max(sp[4], 1e-9); keff0[j] = IK["k"]*sp[6]
        sps = wt_static(keff0[j]); vw0s = max(sps[4], 1e-9)
        for i, psi in enumerate(PSI):
            end, _ = tip_late(sp, a, psi)
            benefit[i, j] = max(np.log10(vw0) - np.log10(max(end[4], 1e-9)), 0)
            persist[i, j] = 1.0 if end[5] > 1.0 else 0.0
            ends, _ = tip_static(sps, keff0[j], psi)
            bstat = max(np.log10(vw0s) - np.log10(max(ends[4], 1e-9)), 0)
            rescue[i, j] = benefit[i, j] - bstat

    # ---- TIMING experiment at strong TIP (psi=15) ----
    PSI_T = 15.0
    vw_noTIP, vw_late, vw_early = [], [], []
    for a in A:
        sp = base_dyn(a); vw_noTIP.append(max(sp[4], 1e-3))
        el, _ = tip_late(sp, a, PSI_T); vw_late.append(max(el[4], 1e-3))
        ee, _ = tip_early(a, PSI_T);    vw_early.append(max(ee[4], 1e-3))
    vw_noTIP, vw_late, vw_early = map(np.array, (vw_noTIP, vw_late, vw_early))

    # ---- kap_eff(t) timecourse: early-TIP at a where chronic immunity is lethal ----
    a_show = A[np.argmin(np.abs(keff0 - 0.7))]
    teval = np.linspace(0, 900, 600)
    _, s_e = integ([T0, 0, 0, 0, 1e-3, 1e2, E_INIT], a_show, PSI_T, teval=teval)
    spc = base_dyn(a_show); _, s_l = integ(list(spc[:5]) + [1e2, spc[6]], a_show, PSI_T, teval=teval)

    np.savez(os.path.join(HERE, "tip_dyn_sweep.npz"),
             A=A, PSI=PSI, benefit=benefit, persist=persist, rescue=rescue,
             keff0=keff0, vw_noTIP=vw_noTIP, vw_late=vw_late, vw_early=vw_early)

    # ---- figures ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    im0 = ax[0].pcolormesh(A, PSI, benefit, shading="auto", cmap="viridis")
    ax[0].contour(A, PSI, persist, levels=[0.5], colors="white", linewidths=2)
    ax0b = ax[0].twiny(); ax0b.set_xlim(keff0[0], keff0[-1])
    ax0b.set_xlabel("-> realized baseline immune pressure k*E0 (/day)")
    ax[0].set_xlabel("immune responsiveness a (/day)"); ax[0].set_ylabel("TIP advantage psi")
    ax[0].set_title("Dynamic immunity: late-TIP marginal benefit\n(white = persistence edge)")
    fig.colorbar(im0, ax=ax[0], label="WT log10 reduction")
    vmax = max(0.5, np.abs(rescue).max())
    im1 = ax[1].pcolormesh(A, PSI, rescue, shading="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax[1].set_xlabel("immune responsiveness a (/day)"); ax[1].set_ylabel("TIP advantage psi")
    ax[1].set_title("RESCUE = dynamic - static-matched\n(red>0 = feedback helps beyond static)")
    fig.colorbar(im1, ax=ax[1], label="extra log10 vs static")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "dynamic_phase.png"), dpi=130)

    fig2, ax2 = plt.subplots(1, 2, figsize=(13, 4.8))
    ax2[0].semilogy(A, vw_noTIP, "k:", lw=2, label="no TIP (immunity only)")
    ax2[0].semilogy(A, vw_late, "tab:orange", lw=2, marker="o", ms=4, label="TIP late (chronic)")
    ax2[0].semilogy(A, vw_early, "tab:green", lw=2, marker="s", ms=4, label="TIP early (acute co-dose)")
    ax2[0].set_xlabel("immune responsiveness a (/day)"); ax2[0].set_ylabel("final WT load Vw")
    ax2[0].set_title(f"TIMING decides it (psi={PSI_T:.0f}): early TIP controls,\nlate TIP fails as immunity rises")
    ax2[0].legend(fontsize=9)
    ax2[1].semilogy(s_e.t, np.clip(s_e.y[4], 1e-3, None), "tab:green", lw=2, label="Vw, early TIP")
    ax2[1].semilogy(s_l.t, np.clip(s_l.y[4], 1e-3, None), "tab:orange", lw=2, label="Vw, late TIP")
    axr = ax2[1].twinx()
    axr.plot(s_e.t, IK["k"]*s_e.y[6], "tab:green", ls="--", alpha=0.7)
    axr.plot(s_l.t, IK["k"]*s_l.y[6], "tab:orange", ls="--", alpha=0.7)
    axr.set_ylabel("k*E (/day), dashed")
    ax2[1].set_xlabel("days"); ax2[1].set_ylabel("WT load Vw (solid)")
    ax2[1].set_title(f"a={a_show:.2f} (chronic k*E0~{keff0[np.argmin(np.abs(keff0-0.7))]:.2f}): early keeps\nimmunity+virus low; late can't start")
    ax2[1].legend(fontsize=9, loc="lower left")
    fig2.tight_layout(); fig2.savefig(os.path.join(HERE, "timing.png"), dpi=130)

    # ---- summary ----
    print("\n--- (1) rescue map (dynamic vs static-matched) ---")
    pos = rescue > 0.3
    print(f"cells where dynamic >0.3 log better than static: {pos.sum()}/{rescue.size} ({pos.mean():.0%})")
    print(f"max extra benefit from feedback: {rescue.max():.2f} log; "
          f"dead-zone reopened? {'NO' if benefit[:, keff0>0.6].max() < 0.5 else 'YES'} "
          f"(max benefit where k*E0>0.6: {benefit[:, keff0>0.6].max() if (keff0>0.6).any() else 0:.2f} log)")
    print("\n--- (2) timing ---")
    for idx in [0, len(A)//2, len(A)-1]:
        print(f"a={A[idx]:.2f} (k*E0~{keff0[idx]:.2f}): noTIP Vw={vw_noTIP[idx]:.2e} | "
              f"late={vw_late[idx]:.2e} ({np.log10(vw_noTIP[idx]/vw_late[idx]):+.2f} log) | "
              f"early={vw_early[idx]:.2e} ({np.log10(vw_noTIP[idx]/vw_early[idx]):+.2f} log)")
    earlywin = np.log10(vw_late/vw_early)
    print(f"\nearly-minus-late WT reduction: mean {earlywin.mean():.2f} log, max {earlywin.max():.2f} log")
    print("wrote dynamic_phase.png, timing.png, tip_dyn_sweep.npz")


if __name__ == "__main__":
    main()
