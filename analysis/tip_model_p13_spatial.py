#!/usr/bin/env python3
"""P1.3 Part 2 — the spatial escape hatch.

Well-mixed (Part 1, stabilized) gives clean antagonism: TIP suppression == antigen
removal == CD8 loss, one axis. The escape hatch: in tissue, antigen presentation and
TIP interference may be SPATIALLY separated. Two compartments:
  P (periphery/blood): TIP works here; systemic viral load = Vw_P (what's measured).
  F (follicle/germinal center): TIP penetrates with efficiency phi (phi=0 => sanctuary).
CD8 are systemic, primed by TOTAL antigen A = (Iw+Id)_P + (Iw+Id)_F, kill productive
cells in both. Free virus diffuses P<->F at rate g; TIP enters F at phi*g.

Question: can the TIP suppress SYSTEMIC WT (Vw_P) while the follicle retains antigen
that keeps CD8 up? i.e. does (peripheral WT) DECOUPLE from (total antigen/CD8),
breaking the well-mixed -1 anti-correlation? If yes, the antagonism is a well-mixed
artifact; if no, it's structural even in tissue.
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tip_model import P, T0, _san
from tip_model_p13_wm import QS, Estar, killrate


def rhs(t, y, psi, phi, g):
    Tp, Iwp, Itp, Idp, Vwp, Vtp, Tf, Iwf, Itf, Idf, Vwf, Vtf = np.clip(y, 0.0, 1e13)
    lam, dT, d, c, p, b, bt, rho = (P[k] for k in
        ("lam", "dT", "d", "c", "p", "b", "bt", "rho"))
    A = Iwp + Idp + Iwf + Idf
    kl = killrate(A)

    def loc(T, Iw, It, Id, Vw, Vt):
        return [lam - dT*T - b*T*Vw - bt*T*Vt,
                b*T*Vw - (d+kl)*Iw - bt*Iw*Vt,
                bt*T*Vt - d*It - b*It*Vw,
                bt*Iw*Vt + b*It*Vw - (d+kl)*Id,
                p*Iw + (1-rho)*p*Id - c*Vw,
                psi*rho*p*Id - c*Vt]
    dP = loc(Tp, Iwp, Itp, Idp, Vwp, Vtp)
    dF = loc(Tf, Iwf, Itf, Idf, Vwf, Vtf)
    dP[4] += g*(Vwf - Vwp); dF[4] += g*(Vwp - Vwf)         # WT diffuses freely
    dP[5] += g*Vtf - phi*g*Vtp; dF[5] += phi*g*Vtp - g*Vtf  # TIP enters F at phi*g
    return dP + dF


ACUTE = [T0, 0, 0, 0, 1e-3, 0, T0, 0, 0, 0, 1e-3, 0]


def run(psi, phi, g, dose=0.0, tmax=1600):
    y0 = list(ACUTE)
    if dose > 0:                      # introduce TIP into periphery at setpoint first
        s1 = solve_ivp(rhs, (0, 400), _san(ACUTE), args=(psi, phi, g),
                       method="LSODA", rtol=1e-7, atol=1e-2, max_step=4.0)
        y0 = _san(s1.y[:, -1]); y0[5] = dose
    s = solve_ivp(rhs, (0, tmax), _san(y0), args=(psi, phi, g),
                  method="LSODA", rtol=1e-7, atol=1e-2, max_step=4.0)
    e = _san(s.y[:, -1])
    return dict(VwP=e[4], VwF=e[10], A=e[1]+e[3]+e[7]+e[9],
                AP=e[1]+e[3], AF=e[7]+e[9])


def main():
    PSI = 22.0
    PHI = np.array([0.0, 0.1, 0.25, 0.5, 1.0])
    G = np.array([0.01, 0.03, 0.1, 0.3, 1.0])
    sysred = np.zeros((len(PHI), len(G)))
    cd8ret = np.zeros((len(PHI), len(G)))
    print(f"psi={PSI:.0f}.  systemic WT reduction (Vw_P) / CD8 retention  vs (phi, g)")
    print("phi\\g   " + "  ".join(f"{g:>6.2f}" for g in G))
    for i, phi in enumerate(PHI):
        cells = []
        for j, g in enumerate(G):
            base = run(PSI, phi, g, dose=0.0)
            tip = run(PSI, phi, g, dose=1e2)
            sysred[i, j] = np.log10(max(base["VwP"], 1e-9)) - np.log10(max(tip["VwP"], 1e-9))
            cd8ret[i, j] = Estar(tip["A"]) / max(Estar(base["A"]), 1e-9)
            cells.append(f"{sysred[i,j]:>4.2f}/{cd8ret[i,j]:>3.2f}")
        print(f"{phi:>4.2f}   " + " ".join(cells))

    # escape = systemic WT down AND CD8 retained (decoupled vs the well-mixed -1 line)
    escape = (sysred >= 0.3) & (cd8ret >= 0.7)
    print(f"\nDECOUPLING/escape cells (systemic WT down >=0.3 log AND CD8 >=70% kept): "
          f"{escape.sum()}/{escape.size}")
    if escape.any():
        ii, jj = np.where(escape)
        b = np.argmax(sysred[ii, jj])
        print(f"  best: phi={PHI[ii[b]]:.2f}, g={G[jj[b]]:.2f} -> "
              f"systemic WT down {sysred[ii[b],jj[b]]:.2f} log, CD8 {cd8ret[ii[b],jj[b]]:.0%} kept")
        print("  => antagonism BREAKS in tissue (follicular antigen sanctuary decouples"
              " systemic suppression from CD8 priming)")
    else:
        print("  => antagonism is STRUCTURAL even with spatial sanctuary (no decoupling)")

    # also report: at the escape point, is the follicle a persistent WT reservoir?
    if escape.any():
        phi, g = PHI[ii[b]], G[jj[b]]
        tip = run(PSI, phi, g, dose=1e2)
        print(f"  at that point: follicular WT Vw_F={tip['VwF']:.2e}, "
              f"follicular antigen A_F={tip['AF']:.2e} (reservoir persists => not a cure)")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    im0 = ax[0].pcolormesh(G, PHI, sysred, shading="auto", cmap="viridis")
    ax[0].set_xscale("log"); ax[0].set_xlabel("P<->F coupling g (/day)")
    ax[0].set_ylabel("TIP follicle access phi"); ax[0].set_title("Systemic WT log-reduction (Vw_P)")
    fig.colorbar(im0, ax=ax[0])
    im1 = ax[1].pcolormesh(G, PHI, cd8ret, shading="auto", cmap="RdYlGn", vmin=0, vmax=1.1)
    ax[1].contour(G, PHI, (escape).astype(float), levels=[0.5], colors="blue", linewidths=2)
    ax[1].set_xscale("log"); ax[1].set_xlabel("P<->F coupling g (/day)")
    ax[1].set_ylabel("TIP follicle access phi")
    ax[1].set_title("CD8 retention (blue = escape: systemic down + CD8 kept)")
    fig.colorbar(im1, ax=ax[1])
    fig.tight_layout(); fig.savefig("/home/seth/dev/hiv-aids-research/analysis/p13_spatial.png", dpi=130)
    print("wrote p13_spatial.png")


if __name__ == "__main__":
    main()
