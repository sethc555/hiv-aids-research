#!/usr/bin/env python3
"""P1.5 — the stable-by-construction build (the re-audit's identified next step).

Every realistic refinement oscillated. P1.3 Part 1 showed QSS *immunity* removes the
effector-lag oscillation, but P1.4 (graded visibility nu<1) re-oscillated — because that
instability is the TIP<->WT COMPETITION cycling through the free-virus stage, not immune
lag. Fix: ALSO adiabatically eliminate the fast free-virus variables (Vw, Vt clear at
c=23/day >> cell turnover d=1/day), collapsing to a 4-variable CELL-ONLY competition with
no fast lag. Then re-run the P1.4 escape test (does an immune-compatible TIP exist:
suppress active WT AND keep CD8 primed?) with the two Simonetti levers (nu, A_def).

QSS free virus:  Vw* = (p/c)(Iw + (1-rho)Id) ,  Vt* = (p/c)(psi*rho*Id)
Cells: dT  = lam - dT*T - b*T*Vw* - b*T*Vt*
       dIw = b*T*Vw* - (d+kw)*Iw - b*Iw*Vt*       kw = k*E*(Aeff)
       dIt = b*T*Vt* - d*It      - b*It*Vw*
       dId = b*Iw*Vt* + b*It*Vw* - (d+kd)*Id       kd = nu*k*E*(Aeff)
QSS immunity: E*(Aeff)=Emax*Aeff/(Aeff+K_E);  Aeff = Iw + nu*Id + A_def.

AUDIT DISCIPLINE: report stability (tail max/min) for every cell; never call an "escape"
real unless it is a converged fixed point.
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tip_model import P, T0, _san
from tip_model_p13_wm import QS

PSI = 22.0
pc = P["p"] / P["c"]                       # virion-per-cell QSS factor


def Estar(A):
    return QS["Emax"] * A / (A + QS["K_E"])


def rhs(t, y, nu, Adef):
    T, Iw, It, Id = np.clip(y, 0.0, 1e13)
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]
    Vw = pc * (Iw + (1 - rho) * Id)        # QSS free WT
    Vt = pc * (PSI * rho * Id)             # QSS free TIP
    E = Estar(Iw + nu * Id + Adef)
    kw, kd = QS["k"] * E, nu * QS["k"] * E
    return [
        lam - dT * T - b * T * Vw - b * T * Vt,
        b * T * Vw - (d + kw) * Iw - b * Iw * Vt,
        b * T * Vt - d * It - b * It * Vw,
        b * Iw * Vt + b * It * Vw - (d + kd) * Id,
    ]


def integ(y0, nu, Adef, tmax=2000, teval=None):
    s = solve_ivp(rhs, (0, tmax), _san(y0), args=(nu, Adef), method="LSODA",
                  rtol=1e-8, atol=1e-3, max_step=5.0, t_eval=teval)
    return s


def vload(y):
    return pc * (y[1] + (1 - P["rho"]) * y[3])   # WT viral load Vw*


def tail_osc(s, win=400):
    vw = np.clip(pc * (s.y[1] + (1 - P["rho"]) * s.y[3]), 1e-6, None)
    t = s.t[s.t > s.t[-1] - win]
    v = vw[s.t > s.t[-1] - win]
    return v.max() / max(v.min(), 1e-9)


def main():
    R0 = P["b"] * pc * T0 / P["d"]
    print(f"reduced model: WT R0 = {R0:.2f}  (should match 8.70)")

    # ---- (1) stability check: does the reduced model stay stable where P1.4 oscillated? ----
    print("\n--- stability check (tail max/min; <1.05 = STABLE) ---")
    for nu, Adef, tag in [(1.0, 0, "no refinements"), (0.42, 11250, "P1.4 'best escape' (was 606x)"),
                          (0.5, 5214, "P1.4 'realistic' (was 84x)"), (0.1, 8000, "strong evasion")]:
        base = integ([T0, 1e-3, 0, 0], 1.0, Adef)
        st = list(base.y[:, -1]); st[2] += 50; st[3] += 50          # seed TIP carriers/dual
        s = integ(st, nu, Adef)
        print(f"  nu={nu}, A_def={Adef:<6} [{tag}]: tail max/min={tail_osc(s):.3f} "
              f"{'STABLE' if tail_osc(s) < 1.05 else 'OSCILLATING'}")

    # ---- (2) escape sweep (nu, A_def) at psi=22 ----
    NU = np.linspace(0.0, 1.0, 13)
    ADEF = np.linspace(0, 15000, 13)
    wtred = np.zeros((len(NU), len(ADEF)))
    cd8ret = np.zeros((len(NU), len(ADEF)))
    stable = np.zeros((len(NU), len(ADEF)))
    Aact0 = None
    for j, Adef in enumerate(ADEF):
        base = integ([T0, 1e-3, 0, 0], 1.0, Adef)
        bend = base.y[:, -1]; vw0 = max(vload(bend), 1e-9)
        E0 = Estar(bend[1] + bend[3] + Adef)
        if j == 0:
            Aact0 = bend[1] + bend[3]
        for i, nu in enumerate(NU):
            st = list(bend); st[2] += 50; st[3] += 50
            s = integ(st, nu, Adef)
            end = s.y[:, -1]
            wtred[i, j] = np.log10(vw0) - np.log10(max(vload(end), 1e-9))
            cd8ret[i, j] = Estar(end[1] + nu * end[3] + Adef) / max(E0, 1e-9)
            stable[i, j] = 1.0 if tail_osc(s) < 1.05 else 0.0

    escape = (wtred >= 0.5) & (cd8ret >= 0.7) & (stable > 0.5)
    np.savez("/home/seth/dev/hiv-aids-research/analysis/p15_sweep.npz",
             NU=NU, ADEF=ADEF, wtred=wtred, cd8ret=cd8ret, stable=stable, Aact0=Aact0)

    print(f"\n--- escape sweep (psi=22) ---")
    print(f"baseline active antigen A0 = {Aact0:.0f}")
    print(f"fraction of grid that is STABLE: {stable.mean():.0%}  (P1.4 was ~0%)")
    print(f"max WT reduction (stable cells only): "
          f"{wtred[stable>0.5].max() if (stable>0.5).any() else 0:.2f} log")
    print(f"IMMUNE-COMPATIBLE & STABLE (WT>=0.5 log AND CD8>=70% AND converged): "
          f"{escape.sum()}/{escape.size} ({'EXISTS' if escape.any() else 'NONE'})")
    if escape.any():
        ii, jj = np.where(escape)
        b = np.argmax(wtred[ii, jj])
        print(f"  best: nu={NU[ii[b]]:.2f}, A_def={ADEF[jj[b]]:.0f} ({ADEF[jj[b]]/Aact0:.1f}x active) "
              f"-> WT down {wtred[ii[b],jj[b]]:.2f} log, CD8 {cd8ret[ii[b],jj[b]]:.0%} kept, STABLE")
        amin = ADEF[jj].min()
        print(f"  min reservoir A_def for any stable escape: {amin:.0f} ({amin/Aact0:.1f}x active)")

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    x = ADEF / Aact0
    for a, M, ttl, cm, vlim in [(ax[0], wtred, "active-WT log-reduction", "viridis", None),
                                 (ax[1], cd8ret, "CD8 retained (E*)", "RdYlGn", (0, 1.2)),
                                 (ax[2], stable, "stable (1) vs oscillating (0)", "Greys", (0, 1))]:
        im = a.pcolormesh(x, NU, M, shading="auto", cmap=cm,
                          vmin=None if not vlim else vlim[0], vmax=None if not vlim else vlim[1])
        a.set_xlabel("reservoir A_def (x active antigen)"); a.set_ylabel("dual-cell visibility nu")
        a.set_title(ttl); fig.colorbar(im, ax=a)
    if escape.any():
        ax[0].contour(x, NU, escape.astype(float), levels=[0.5], colors="red", linewidths=2.5)
    fig.suptitle("P1.5 stable build (QSS virus + QSS immunity): escape test (red = immune-compatible & converged)")
    fig.tight_layout(); fig.savefig("/home/seth/dev/hiv-aids-research/analysis/p15_escape.png", dpi=130)
    print("\nwrote p15_escape.png, p15_sweep.npz")


if __name__ == "__main__":
    main()
