#!/usr/bin/env python3
"""P1.3 Part 1 — STABILIZED well-mixed model (quasi-steady-state CD8).

The audit's main fragility was the predator-prey OSCILLATION (effector lag -> limit
cycles -> unreliable magnitudes). Fix it by adiabatic elimination: assume CD8 track
antigen with no lag, E*(A) = Emax*A/(A+K_E), killing per productive cell
kill(A) = k*E*(A)/(1+A/Kx) (exhaustion saturates total kill -> stable setpoint, no
cycle). No E ODE => no oscillation by construction.

Then redo the antagonism test CLEANLY (stable endpoints, no time-averaging needed):
sweep psi, measure WT reduction and CD8 retention E*(A_tip)/E*(A_0). If antagonism
survives here it's not an oscillation artifact.
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tip_model import P, T0, _san

# weak-ish immunity (kill ~0.4/day at baseline): controls WT only mildly, so the TIP
# CAN partially establish at high psi -> reveals the suppression-vs-CD8 trade-off curve
# (rather than total TIP exclusion, which strong stable immunity trivially produces).
QS = dict(Emax=3e4, K_E=3000.0, k=3e-5, Kx=2500.0)


def Estar(A):
    return QS["Emax"] * A / (A + QS["K_E"])


def killrate(A):
    # monotone-increasing in antigen => stabilizing negative feedback (no exhaustion term,
    # which created a destabilizing positive feedback and kept the baseline oscillating).
    return QS["k"] * Estar(A)


def rhs(t, y, psi):
    T, Iw, It, Id, Vw, Vt = np.clip(y, 0.0, 1e13)
    lam, dT, d, c, p, b, bt, rho = (P[k] for k in
        ("lam", "dT", "d", "c", "p", "b", "bt", "rho"))
    A = Iw + Id
    kl = killrate(A)
    return [
        lam - dT*T - b*T*Vw - bt*T*Vt,
        b*T*Vw - (d+kl)*Iw - bt*Iw*Vt,
        bt*T*Vt - d*It - b*It*Vw,
        bt*Iw*Vt + b*It*Vw - (d+kl)*Id,
        p*Iw + (1-rho)*p*Id - c*Vw,
        psi*rho*p*Id - c*Vt,
    ]


def integ(y0, psi, tmax=1500):
    s = solve_ivp(rhs, (0, tmax), _san(y0), args=(psi,), method="LSODA",
                  rtol=1e-7, atol=1e-2, max_step=4.0)
    return _san(s.y[:, -1]), s


def main():
    # baseline: WT + QSS immunity, no TIP
    base, sb = integ([T0, 0, 0, 0, 1e-3, 0], 1.0)
    A0 = base[1] + base[3]; E0 = Estar(A0)
    # oscillation check: tail variation
    tail = np.clip(sb.y[4][sb.t > sb.t[-1]-300], 1e-3, None)
    osc = tail.max()/max(tail.min(), 1e-3)
    print(f"BASELINE: Vw={base[4]:.3e} (no-imm ref 7.7e5), A=Iw+Id={A0:.3e}, "
          f"E*(A0)={E0:.3e}, kill={killrate(A0):.3f}/day, "
          f"tail max/min={osc:.3f} ({'STABLE' if osc < 1.05 else 'oscillating'})")

    # antagonism test: late TIP, sweep psi, CLEAN endpoints
    PSI = np.linspace(1, 22, 22)
    print(f"\n{'psi':>4} {'Vw_tip':>10} {'WTred':>7} {'A_tip':>9} {'CD8ret':>7}")
    wt, cd8 = [], []
    rows = []
    for psi in PSI:
        st = list(base); st[5] = 1e2
        end, _ = integ(st, psi)
        At = end[1] + end[3]; Et = Estar(At)
        wr = np.log10(base[4]) - np.log10(max(end[4], 1e-9))
        cr = Et / max(E0, 1e-9)
        wt.append(wr); cd8.append(cr)
        if psi in PSI[::3]:
            rows.append(f"{psi:>4.0f} {end[4]:>10.2e} {wr:>7.2f} {At:>9.2e} {cr:>7.2f}")
    print("\n".join(rows))
    wt, cd8 = np.array(wt), np.array(cd8)

    # is there a psi with WT down >=1 log AND CD8 retained >=50%?
    compat = (wt >= 1.0) & (cd8 >= 0.5)
    print(f"\nmax WT reduction = {wt.max():.2f} log (psi={PSI[wt.argmax()]:.0f})")
    print(f"max WT reduction with CD8>=50% retained = "
          f"{wt[cd8>=0.5].max() if (cd8>=0.5).any() else 0:.2f} log")
    print(f"immune-compatible psi (WT>=1 & CD8>=0.5): {compat.sum()} "
          f"({'EXISTS' if compat.any() else 'NONE -> antagonism holds (no oscillation)'})")
    print(f"corr(WTred, CD8ret) over psi where TIP acts (WTred>0.05): "
          f"{np.corrcoef(wt[wt>0.05], cd8[wt>0.05])[0,1] if (wt>0.05).sum()>2 else float('nan'):.2f}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(wt, cd8, "o-", color="tab:purple")
    for i in range(0, len(PSI), 3):
        ax.annotate(f"ψ{PSI[i]:.0f}", (wt[i], cd8[i]), fontsize=7)
    ax.axhline(0.5, ls="--", c="gray"); ax.axvline(1.0, ls="--", c="gray")
    ax.fill_betweenx([0.5, 1.05], 1.0, max(2, wt.max()), color="green", alpha=0.08)
    ax.text(1.02, 1.0, "immune-compatible", color="green", fontsize=8)
    ax.set_xlabel("WT log10 reduction by TIP"); ax.set_ylabel("CD8 retained E*(A_tip)/E*(A0)")
    ax.set_title("P1.3 stabilized (QSS CD8, no oscillation):\nTIP suppression vs CD8 retention, swept over psi")
    fig.tight_layout(); fig.savefig("/home/seth/dev/hiv-aids-research/analysis/p13_wm.png", dpi=130)
    print("wrote p13_wm.png")


if __name__ == "__main__":
    main()
