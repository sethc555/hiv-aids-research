#!/usr/bin/env python3
"""P1 — within-host TIP vs immunity model.

Tests the Dodd & de Boer (2026, J Theor Biol) claim that "even a moderate immune
response against virally infected cells drastically decreases the range of
parameter values for which [TIP] therapy is effective."

Mechanism under test: a therapeutic interfering particle (TIP) is a conditionally-
replicating defective HIV genome. It produces particles ONLY in cells that are
ALSO infected by wild-type (dually-infected, I_d) — where it hijacks WT proteins
and diverts packaging to itself (fraction rho), amplified by a mobilization
advantage psi. Immunity (CTL) kills productive infected cells (I_w, I_d), so it
shortens the TIP's only production window. Does that erase the TIP's benefit?

Within-host ODE (population level), per mL per day:
  T'  = lam - dT*T - b*T*Vw - bt*T*Vt
  Iw' = b*T*Vw            - (d+kap)*Iw - bt*Iw*Vt        # WT-only productive
  It' = bt*T*Vt           - d*It        - b*It*Vw         # TIP-only carrier (no output)
  Id' = bt*Iw*Vt + b*It*Vw - (d+kap)*Id                   # dual: makes WT and TIP
  Vw' = p*Iw + (1-rho)*p*Id - c*Vw                        # WT output reduced in I_d
  Vt' = psi*rho*p*Id        - c*Vt                        # TIP output only from I_d
kap = immune killing pressure on PRODUCTIVE cells (the de Boer knob).
psi = TIP mobilization advantage (the design knob).

Experiment: marginal TIP benefit = log10(WT setpoint without TIP) -
log10(WT setpoint with TIP), measured at each immune level kap. If immunity
"drastically shrinks" efficacy, the benefit collapses as kap rises unless psi is
large -> a design spec psi_required(kap).
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- parameters (standard Perelson-class within-host HIV) ----
P = dict(lam=1e4, dT=0.01, d=1.0, c=23.0, p=2000.0,
         b=1e-7,        # WT infection rate -> R0 ~ 8.7
         bt=1e-7,       # TIP infectivity (complemented Env, ~= WT)
         rho=0.90,      # fraction of dual-cell output diverted to TIP
         )
T0 = P["lam"] / P["dT"]          # uninfected target setpoint = 1e6


def rhs(t, y, kap, psi):
    T, Iw, It, Id, Vw, Vt = np.clip(y, 0.0, 1e13)   # populations >=0, finite
    lam, dT, d, c, p, b, bt, rho = (P[k] for k in
        ("lam", "dT", "d", "c", "p", "b", "bt", "rho"))
    return [
        lam - dT*T - b*T*Vw - bt*T*Vt,
        b*T*Vw - (d+kap)*Iw - bt*Iw*Vt,
        bt*T*Vt - d*It - b*It*Vw,
        bt*Iw*Vt + b*It*Vw - (d+kap)*Id,
        p*Iw + (1-rho)*p*Id - c*Vw,
        psi*rho*p*Id - c*Vt,
    ]


def _san(a):
    return np.clip(np.nan_to_num(np.asarray(a, float), nan=0.0, posinf=1e13, neginf=0.0), 0, 1e13)


def integrate(y0, kap, psi, tmax):
    y0 = _san(y0)
    s = solve_ivp(rhs, (0, tmax), y0, args=(kap, psi), method="LSODA",
                  rtol=1e-6, atol=1e-2, dense_output=False, max_step=2.0)
    return _san(s.y[:, -1]), s


def wt_setpoint(kap, tmax=700):
    """WT-only chronic setpoint at immune pressure kap (no TIP)."""
    y0 = [T0, 0, 0, 0, 1e-3, 0]
    end, _ = integrate(y0, kap, 1.0, tmax)
    return end


def tip_arm(state_wt, kap, psi, tip_dose=1e2, tmax=700):
    """From the WT setpoint state, introduce a TIP bolus; return endpoint."""
    y0 = list(state_wt); y0[5] = tip_dose
    end, s = integrate(y0, kap, psi, tmax)
    return end, s


def main():
    R0 = P["b"]*T0*P["p"]/(P["d"]*P["c"])
    print(f"WT basic reproductive ratio R0 = {R0:.2f}")

    # sanity: baseline WT setpoint, no immunity
    base = wt_setpoint(0.0)
    print(f"WT-only setpoint (kap=0): Vw={base[4]:.3e}, Iw={base[1]:.3e}, T={base[0]:.3e}")

    # invasion threshold (analytic, kap=0): psi ~ 7.5. NB: R0_TIP is sub-linear (~sqrt(psi),
    # two-stage TIP cycle); "0.134*psi" is only a local fit at threshold (audit 2026-06-20).
    # sanity: TIP above threshold should suppress WT
    end0, _ = tip_arm(base, 0.0, 12.0)
    print(f"TIP arm (kap=0, psi=12): Vw={end0[4]:.3e}, Vt={end0[5]:.3e}  "
          f"-> {np.log10(base[4]/max(end0[4],1e-9)):.2f} log WT drop")

    # ---- phase sweep over (kap, psi) ----
    KAP = np.linspace(0.0, 1.2, 25)
    PSI = np.linspace(1.0, 20.0, 25)
    logdrop = np.zeros((len(PSI), len(KAP)))
    persist = np.zeros((len(PSI), len(KAP)))
    for j, kap in enumerate(KAP):
        sp = wt_setpoint(kap)
        vw0 = max(sp[4], 1e-9)
        for i, psi in enumerate(PSI):
            end, _ = tip_arm(sp, kap, psi)
            vw_tip = max(end[4], 1e-9)
            logdrop[i, j] = np.log10(vw0) - np.log10(vw_tip)
            persist[i, j] = 1.0 if end[5] > 1.0 else 0.0
    logdrop = np.clip(logdrop, 0, None)

    # design spec: min psi giving >=1 log marginal benefit, per kap
    psi_req = []
    for j in range(len(KAP)):
        ok = np.where(logdrop[:, j] >= 1.0)[0]
        psi_req.append(PSI[ok[0]] if len(ok) else np.nan)
    psi_req = np.array(psi_req)

    np.savez("/home/seth/dev/hiv-aids-research/analysis/tip_sweep.npz",
             KAP=KAP, PSI=PSI, logdrop=logdrop, persist=persist, psi_req=psi_req)

    # ---- figure 1: phase diagram ----
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    im = ax.pcolormesh(KAP, PSI, logdrop, shading="auto", cmap="viridis")
    cs = ax.contour(KAP, PSI, persist, levels=[0.5], colors="white", linewidths=2)
    ax.plot(KAP, psi_req, "r-", lw=2.2, label="psi required for >=1 log WT drop")
    ax.set_xlabel("immune killing pressure on productive cells  kap  (/day)")
    ax.set_ylabel("TIP mobilization advantage  psi")
    ax.set_title("TIP marginal benefit collapses as immunity rises\n"
                 "(color = log10 WT setpoint reduction from TIP; white = TIP persistence edge)")
    fig.colorbar(im, ax=ax, label="marginal WT log10 reduction by TIP")
    ax.legend(loc="upper right", framealpha=0.9)
    fig.tight_layout()
    fig.savefig("/home/seth/dev/hiv-aids-research/analysis/phase_diagram.png", dpi=130)

    # ---- figure 2: representative timecourses ----
    fig2, ax2 = plt.subplots(figsize=(7.5, 4.8))
    for kap, c_ in [(0.0, "tab:blue"), (0.4, "tab:orange"), (0.8, "tab:red")]:
        sp = wt_setpoint(kap)
        _, s = tip_arm(sp, kap, 12.0)
        vw = np.clip(np.nan_to_num(s.y[4], posinf=1e13), 1e-3, None)
        ax2.semilogy(s.t, vw, c_, lw=2, label=f"Vw, TIP+immunity kap={kap}")
        ax2.axhline(max(sp[4], 1e-3), color=c_, ls=":", lw=1, alpha=0.6)
    ax2.set_xlabel("days after TIP dose"); ax2.set_ylabel("WT viral load Vw (/mL)")
    ax2.set_title("WT load after TIP (psi=12); dotted = pre-TIP setpoint at each immune level")
    ax2.legend(fontsize=8); fig2.tight_layout()
    fig2.savefig("/home/seth/dev/hiv-aids-research/analysis/timecourse.png", dpi=130)

    # ---- printed summary ----
    print("\n--- phase summary ---")
    for jj in [0, 8, 16, 24]:
        kap = KAP[jj]
        col = logdrop[:, jj]
        print(f"kap={kap:.2f}/day : max TIP benefit={col.max():.2f} log "
              f"(at psi={PSI[col.argmax()]:.1f}); "
              f"psi req for 1-log = {psi_req[jj] if not np.isnan(psi_req[jj]) else 'NONE in range'}")
    # how fast does the efficacy window shrink?
    frac0 = (logdrop[:, 0] >= 1.0).mean()
    fracmid = (logdrop[:, 12] >= 1.0).mean()
    print(f"\nfraction of psi-range giving >=1 log benefit: "
          f"kap=0 -> {frac0:.0%};  kap={KAP[12]:.2f} -> {fracmid:.0%}")
    print("wrote phase_diagram.png, timecourse.png, tip_sweep.npz")


if __name__ == "__main__":
    main()
