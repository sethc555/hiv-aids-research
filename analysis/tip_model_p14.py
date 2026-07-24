#!/usr/bin/env python3
"""P1.4 — does the Simonetti 5'-leader-NSV biology break the antagonism?

Built on the STABLE P1.3 Part-1 model (QSS CD8, no oscillation). Two refinements,
both grounded in Simonetti (JCI 2023 / Nat Commun 2026):

  (A) Graded immune visibility of dual (TIP-producing) cells, nu in [0,1]:
      defective-genome cells make less WT antigen (esp. Env) + carry CTL-escape
      mutations -> killed less AND prime CD8 less. I_d killed at nu*kill, and I_d
      contributes nu to the antigen that sustains CD8.
  (B) Persistent defective-clone antigen reservoir A_def: clonally-expanded
      5'-leader-defective proviruses (~10^6-10^7 cells, "1 in 3 infected") transcribe
      and present antigen INDEPENDENT of active WT replication -> a CD8-sustaining
      antigen the TIP cannot remove.

A_eff = Iw + nu*Id + A_def  (drives CD8: E*(A_eff)); kill_w = k*E*, kill_d = nu*k*E*.
Escape = TIP lowers active WT (Vw) AND CD8 stays primed (E* retained). Sweep (A_def, nu)
at strong psi=22. Baseline active antigen ~5e3, so A_def~5e3 ~ "half the infected pool"
is the Simonetti-realistic regime.
"""
import os
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tip_model import P, T0, _san
from tip_model_p13_wm import QS

HERE = os.path.dirname(os.path.abspath(__file__))   # path-portable (AUDIT2 #14)

PSI = 22.0


def Estar(Aeff):
    return QS["Emax"] * Aeff / (Aeff + QS["K_E"])


def rhs(t, y, nu, Adef):
    T, Iw, It, Id, Vw, Vt = np.clip(y, 0.0, 1e13)
    lam, dT, d, c, p, b, bt, rho = (P[k] for k in
        ("lam", "dT", "d", "c", "p", "b", "bt", "rho"))
    Aeff = Iw + nu*Id + Adef
    kill = QS["k"] * Estar(Aeff)
    return [
        lam - dT*T - b*T*Vw - bt*T*Vt,
        b*T*Vw - (d+kill)*Iw - bt*Iw*Vt,
        bt*T*Vt - d*It - b*It*Vw,
        bt*Iw*Vt + b*It*Vw - (d + nu*kill)*Id,
        p*Iw + (1-rho)*p*Id - c*Vw,
        PSI*rho*p*Id - c*Vt,
    ]


def integ(y0, nu, Adef, tmax=1500):
    s = solve_ivp(rhs, (0, tmax), _san(y0), args=(nu, Adef), method="LSODA",
                  rtol=1e-7, atol=1e-2, max_step=4.0)
    return _san(s.y[:, -1])


def main():
    ADEF = np.linspace(0, 15000, 13)
    NU = np.linspace(0.0, 1.0, 13)
    wtred = np.zeros((len(NU), len(ADEF)))
    cd8ret = np.zeros((len(NU), len(ADEF)))
    Aactive0 = None
    for j, Adef in enumerate(ADEF):
        base = integ([T0, 0, 0, 0, 1e-3, 0], 1.0, Adef)   # nu irrelevant (Id=0 at baseline)
        Vw0 = max(base[4], 1e-9)
        Aeff0 = base[1] + base[3] + Adef; E0 = Estar(Aeff0)
        if j == 0:
            Aactive0 = base[1] + base[3]
        for i, nu in enumerate(NU):
            st = list(base); st[5] = 1e2
            end = integ(st, nu, Adef)
            Aeff_t = end[1] + nu*end[3] + Adef
            wtred[i, j] = np.log10(Vw0) - np.log10(max(end[4], 1e-9))
            cd8ret[i, j] = Estar(Aeff_t) / max(E0, 1e-9)

    escape = (wtred >= 0.5) & (cd8ret >= 0.7)
    np.savez(os.path.join(HERE, "p14_sweep.npz"),
             ADEF=ADEF, NU=NU, wtred=wtred, cd8ret=cd8ret, Aactive0=Aactive0)

    print(f"baseline active antigen A0 (Adef=0) = {Aactive0:.0f}  "
          f"(so Adef~{Aactive0:.0f} ≈ 'half the infected pool', Simonetti-realistic)")
    print(f"\nmax WT reduction anywhere = {wtred.max():.2f} log")
    print(f"IMMUNE-COMPATIBLE (WT down >=0.5 log AND CD8 >=70% kept): "
          f"{escape.sum()}/{escape.size} ({'EXISTS' if escape.any() else 'NONE'})")
    if escape.any():
        ii, jj = np.where(escape)
        bi = np.argmax(wtred[ii, jj])
        a_esc = ADEF[jj[ii.argsort()][0]] if False else ADEF[jj][np.argmin(ADEF[jj])]
        print(f"  best: nu={NU[ii[bi]]:.2f}, A_def={ADEF[jj[bi]]:.0f} "
              f"({ADEF[jj[bi]]/Aactive0:.1f}x active) -> WT down {wtred[ii[bi],jj[bi]]:.2f} log, "
              f"CD8 {cd8ret[ii[bi],jj[bi]]:.0%} kept")
        print(f"  minimum A_def for ANY escape: {ADEF[jj].min():.0f} "
              f"({ADEF[jj].min()/Aactive0:.1f}x active antigen)")
        # is escape reachable at REALISTIC A_def (<= ~2x active) ?
        realistic = escape & (np.broadcast_to(ADEF, wtred.shape) <= 2*Aactive0)
        print(f"  escape at Simonetti-realistic A_def (<=2x active): "
              f"{'YES' if realistic.any() else 'NO (needs implausibly large reservoir)'}")
    # role of each lever alone
    print(f"\nlever isolation:")
    print(f"  graded killing only (A_def=0): max WT red {wtred[:,0].max():.2f} log, "
          f"min CD8ret among WT-down cells = "
          f"{cd8ret[:,0][wtred[:,0]>=0.3].min() if (wtred[:,0]>=0.3).any() else float('nan'):.2f}")
    jmid = np.argmin(np.abs(ADEF - Aactive0))
    print(f"  reservoir at ~1x active (A_def={ADEF[jmid]:.0f}): at nu=0, "
          f"WT red {wtred[0,jmid]:.2f} log, CD8 {cd8ret[0,jmid]:.0%} kept")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    im0 = ax[0].pcolormesh(ADEF/Aactive0, NU, wtred, shading="auto", cmap="viridis")
    ax[0].contour(ADEF/Aactive0, NU, escape.astype(float), levels=[0.5], colors="red", linewidths=2)
    ax[0].set_xlabel("defective-clone reservoir A_def (x active antigen)")
    ax[0].set_ylabel("dual-cell immune visibility nu")
    ax[0].set_title("Active-WT log-reduction by TIP\n(red = immune-compatible: WT down + CD8 kept)")
    fig.colorbar(im0, ax=ax[0])
    im1 = ax[1].pcolormesh(ADEF/Aactive0, NU, cd8ret, shading="auto", cmap="RdYlGn", vmin=0, vmax=1.2)
    ax[1].set_xlabel("defective-clone reservoir A_def (x active antigen)")
    ax[1].set_ylabel("dual-cell immune visibility nu")
    ax[1].set_title("CD8 retained vs no-TIP  (E*(A_eff))")
    fig.colorbar(im1, ax=ax[1])
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "p14_escape.png"), dpi=130)
    print("\nwrote p14_escape.png, p14_sweep.npz")


if __name__ == "__main__":
    main()
