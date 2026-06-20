#!/usr/bin/env python3
"""P2 — the two-pool (decoupled-antigen) test. Attacks the central audit caveat.

The multi-agent audit's load-bearing finding (AUDIT.md point #1): the TIP<->CD8
antagonism is *partly by-construction* because P1-P1.5 prime CD8 from a SINGLE antigen
pool that is proportional to ACTIVE infection (A = Iw + nu*Id [+ static A_def]). So
"TIP lowers active infection -> starves antigen-driven CD8" is semi-definitional.

P1.5 already carried a *static* decoupled offset A_def, but it could not separate two
explanations for why no stable immune-compatible TIP appeared:
  (i)  the single-pool construction (antagonism is encoded), vs
  (ii) the TIP's intrinsic immune-fragility (it lives at R0~1, so any CD8 that its
       dual-infected factories present sinks it) -- which would survive decoupling.

P2 makes the second pool a DYNAMIC, self-maintaining latent/defective reservoir L,
grounded in Simonetti (JCI 2023 / Nat Commun 2026): ~1e6-1e7-cell defective clones that
persist by antigen-independent CLONAL PROLIFERATION and present HLA antigen WITHOUT
producing infectious virus. Crucially L is:
  - TIP-INDEPENDENT  (a TIP interferes with WT packaging in dual cells; it does nothing
    to a defective clone that only transcribes antigen), and
  - decoupled from active WT on a SLOW timescale (it proliferates, it does not track Iw).

So CD8 is now primed by TWO pools:  Aeff = Iw + nu*Id  +  L   (active + reservoir).
A TIP can suppress active WT (Iw) while the reservoir L keeps CD8 primed. If the
antagonism is merely the single-pool construction, decoupling should now ADMIT a stable
immune-compatible TIP. If it is immune-fragility, it will NOT -- because a primed CD8
pool still kills the TIP's dual-infected factories (kd = nu*k*E).

Second lever this isolates: reservoir TIMESCALE. A large, slow antigen floor pins CD8
~constant, which should DAMP the P1.4/P1.5 immune-feedback oscillation (collapsing the
system toward the stable P1-static limit). We test g_slow (months) vs g_fast (~static
A_def, the P1.5 control) to see whether decoupling buys stability, an escape, or neither.

AUDIT DISCIPLINE (inherited from P1.5): report tail max/min for every run; never call an
"escape" real unless it is a converged fixed point (tail osc < 1.05).
"""
import os
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tip_model import P, T0, _san
from tip_model_p13_wm import QS, Estar

HERE = os.path.dirname(os.path.abspath(__file__))     # write next to this file (no abs path)
PSI = 22.0
pc = P["p"] / P["c"]                                   # virion-per-cell QSS factor
S_SEED = 1e-4                                          # reservoir nucleation from active infection


def rhs(t, y, nu, R, g):
    """4 active-compartment cells (QSS virus) + 1 dynamic reservoir antigen pool L."""
    T, Iw, It, Id, L = np.clip(y, 0.0, 1e13)
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]
    Vw = pc * (Iw + (1 - rho) * Id)                    # QSS free WT
    Vt = pc * (PSI * rho * Id)                         # QSS free TIP
    Aeff = Iw + nu * Id + L                            # TWO POOLS: active + reservoir
    E = Estar(Aeff)
    kw, kd = QS["k"] * E, nu * QS["k"] * E
    return [
        lam - dT * T - b * T * Vw - b * T * Vt,
        b * T * Vw - (d + kw) * Iw - b * Iw * Vt,
        b * T * Vt - d * It - b * It * Vw,
        b * Iw * Vt + b * It * Vw - (d + kd) * Id,
        g * L * (1.0 - L / R) + S_SEED * (Iw + Id),    # self-maintaining clone (TIP-independent)
    ]


def integ(y0, nu, R, g, tmax=3000, teval=None):
    s = solve_ivp(rhs, (0, tmax), _san(y0), args=(nu, R, g), method="LSODA",
                  rtol=1e-8, atol=1e-3, max_step=5.0, t_eval=teval)
    return s


def vload(y):
    return pc * (y[1] + (1 - P["rho"]) * y[3])         # WT viral load Vw*


def tail_osc(s, win=400):
    vw = np.clip(pc * (s.y[1] + (1 - P["rho"]) * s.y[3]), 1e-6, None)
    v = vw[s.t > s.t[-1] - win]
    return v.max() / max(v.min(), 1e-9)


def run_grid(g, label):
    """Escape sweep over (reservoir size R, dual-cell visibility nu) at fixed psi."""
    NU = np.linspace(0.0, 1.0, 13)
    RES = np.linspace(0.0, 15000.0, 13)                # reservoir antigen pool size R
    wtred = np.zeros((len(NU), len(RES)))
    cd8ret = np.zeros((len(NU), len(RES)))
    stable = np.zeros((len(NU), len(RES)))
    A0 = None
    for j, R in enumerate(RES):
        Rsafe = max(R, 1.0)
        # baseline: WT establishes; reservoir self-maintains to ~R (no TIP)
        base = integ([T0, 1e-3, 0, 0, 1.0], 1.0, Rsafe, g)
        bend = base.y[:, -1]
        vw0 = max(vload(bend), 1e-9)
        E0 = Estar(bend[1] + bend[3] + bend[4])
        if j == 0:
            A0 = bend[1] + bend[3]                      # baseline ACTIVE antigen (Iw+Id), R=0
        for i, nu in enumerate(NU):
            st = list(bend); st[2] += 50; st[3] += 50   # seed TIP carriers + dual cells
            s = integ(st, nu, Rsafe, g)
            end = s.y[:, -1]
            wtred[i, j] = np.log10(vw0) - np.log10(max(vload(end), 1e-9))
            cd8ret[i, j] = Estar(end[1] + nu * end[3] + end[4]) / max(E0, 1e-9)
            stable[i, j] = 1.0 if tail_osc(s) < 1.05 else 0.0
    escape = (wtred >= 0.5) & (cd8ret >= 0.7) & (stable > 0.5)
    return dict(NU=NU, RES=RES, wtred=wtred, cd8ret=cd8ret, stable=stable,
                escape=escape, A0=A0, label=label, g=g)


def report(G):
    A0 = G["A0"]
    print(f"\n=== {G['label']} (g={G['g']}/day; psi={PSI}) ===")
    print(f"baseline ACTIVE antigen A0 (R=0) = {A0:.0f}")
    print(f"STABLE fraction of grid: {G['stable'].mean():.0%}  (P1.5 static was ~ small)")
    st = G["stable"] > 0.5
    print(f"max WT reduction (stable cells only): {G['wtred'][st].max() if st.any() else 0:.2f} log")
    e = G["escape"]
    print(f"IMMUNE-COMPATIBLE & STABLE (WT>=0.5 log AND CD8>=70% AND converged): "
          f"{e.sum()}/{e.size}  ->  {'EXISTS' if e.any() else 'NONE'}")
    if e.any():
        ii, jj = np.where(e)
        b = np.argmax(G["wtred"][ii, jj])
        print(f"  best escape: nu={G['NU'][ii[b]]:.2f}, R={G['RES'][jj[b]]:.0f} "
              f"({G['RES'][jj[b]]/A0:.1f}x active antigen) -> WT down "
              f"{G['wtred'][ii[b],jj[b]]:.2f} log, CD8 {G['cd8ret'][ii[b],jj[b]]:.0%} kept, STABLE")


def main():
    R0 = P["b"] * pc * T0 / P["d"]
    print(f"reduced two-pool model: WT R0 = {R0:.2f}  (should match 8.70)")

    # the decisive comparison: a SLOW reservoir (decoupled, months) vs a FAST one
    # (g=2/day ~ instantaneously tracks -> reproduces P1.5's static A_def control).
    G_slow = run_grid(g=0.02, label="SLOW reservoir (half-life ~months, decoupled)")
    G_fast = run_grid(g=2.0, label="FAST reservoir (~static A_def; P1.5 control)")
    report(G_slow)
    report(G_fast)

    np.savez(os.path.join(HERE, "p2_sweep.npz"),
             NU=G_slow["NU"], RES=G_slow["RES"], A0=G_slow["A0"],
             wtred_slow=G_slow["wtred"], cd8ret_slow=G_slow["cd8ret"], stable_slow=G_slow["stable"],
             wtred_fast=G_fast["wtred"], cd8ret_fast=G_fast["cd8ret"], stable_fast=G_fast["stable"])

    # ---- figure: slow vs fast, three panels each ----
    fig, axes = plt.subplots(2, 3, figsize=(17, 9.5))
    for row, G in [(0, G_slow), (1, G_fast)]:
        x = G["RES"] / G["A0"]
        panels = [(G["wtred"], "active-WT log-reduction", "viridis", None),
                  (G["cd8ret"], "CD8 retained (E*/E0)", "RdYlGn", (0, 1.2)),
                  (G["stable"], "stable (1) vs oscillating (0)", "Greys", (0, 1))]
        for k, (M, ttl, cm, vlim) in enumerate(panels):
            a = axes[row, k]
            im = a.pcolormesh(x, G["NU"], M, shading="auto", cmap=cm,
                              vmin=None if not vlim else vlim[0], vmax=None if not vlim else vlim[1])
            a.set_xlabel("reservoir antigen R (x active antigen A0)")
            a.set_ylabel("dual-cell visibility nu")
            a.set_title(f"[{'SLOW' if row == 0 else 'FAST'}] {ttl}")
            fig.colorbar(im, ax=a)
            if G["escape"].any():
                a.contour(x, G["NU"], G["escape"].astype(float), levels=[0.5],
                          colors="red", linewidths=2.5)
    fig.suptitle("P2 two-pool decoupled-antigen test: does a reservoir antigen floor admit a "
                 "stable immune-compatible TIP?  (red = immune-compatible & converged)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "p2_escape.png"), dpi=130)
    print("\nwrote p2_escape.png, p2_sweep.npz")


if __name__ == "__main__":
    main()
