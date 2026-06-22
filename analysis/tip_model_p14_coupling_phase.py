#!/usr/bin/env python3
"""P14 — the coupling phase diagram: WHERE (and how sharply) a coupled TIP flips neutral->helpful.

Pre-work for contacting TIP modelers. P11 showed a coupled TIP helps; P14 makes that a phase
boundary a modeler can scrutinize: the TIP's effect on durable post-treatment control as a
CONTINUOUS function of coupling chi (fraction of reservoir reactivations carrying the TIP) and
immune strength kf. Two things fall out:

  1. At chi=0 the TIP is inert -> this is exactly the de Boer limit (no TIP-reservoir interaction;
     immunity dominates). So our model REDUCES to de Boer's "immunity-only" outcome as chi->0,
     and departs from it continuously as chi rises. (Reduction-to-prior-model check.)
  2. There is a coupling threshold chi*(kf): below it the TIP is ~neutral, above it it helps,
     and the help is largest at intermediate immune strength (where control is marginal and a
     small push matters most).

Output: a (chi x kf) grid of the TIP effect on P(control), the de Boer-limit column, chi* per kf,
and a clean phase figure. Run memory-capped.
"""
import os
import numpy as np
from tip_model_p11_coupled import simulate_coupled

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    N = 180
    CHI = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    KF = np.array([9.0, 10.0, 11.0, 12.0, 13.0])
    print(f"P14 coupling phase diagram (N={N}). TIP effect on durable control = P(ctrl|chi) - P(ctrl|chi=0).")
    print("chi=0 column = the de Boer limit (no TIP-reservoir coupling).\n")
    ctrl = np.zeros((len(CHI), len(KF)))
    for j, kf in enumerate(KF):
        for i, chi in enumerate(CHI):
            ctrl[i, j] = simulate_coupled(N, kf, chi=chi, t_ati=500.0, seed=int(kf * 10 + chi * 11))[0]
    effect = ctrl - ctrl[0:1, :]          # TIP effect vs the chi=0 (de Boer) baseline at each kf

    print("P(control) grid (rows=chi, cols=kf):")
    print("        " + "  ".join(f"kf{int(k)}" for k in KF))
    for i, chi in enumerate(CHI):
        print(f"chi={chi:.1f} " + "  ".join(f"{100*ctrl[i,j]:4.0f}%" for j in range(len(KF))))
    print("\nTIP EFFECT (pts above the chi=0 / de Boer baseline):")
    print("        " + "  ".join(f"kf{int(k)}" for k in KF))
    for i, chi in enumerate(CHI):
        print(f"chi={chi:.1f} " + "  ".join(f"{100*effect[i,j]:+4.0f}" for j in range(len(KF))))

    # coupling threshold chi* per kf: smallest chi giving >= +8 pts (a meaningful help)
    print("\ncoupling threshold chi* (min coupling for >=+8 pts help), per immune strength:")
    for j, kf in enumerate(KF):
        idx = np.where(100 * effect[:, j] >= 8)[0]
        cs = f"{CHI[idx[0]]:.1f}" if len(idx) else ">1.0 (no help in range)"
        print(f"  kf={kf:.0f}: chi* = {cs}   (max help {100*effect[:,j].max():+.0f} pts at chi=1)")

    np.savez(os.path.join(HERE, "p14_coupling_phase.npz"), CHI=CHI, KF=KF, ctrl=ctrl, effect=effect)

    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    im0 = ax[0].imshow(100 * ctrl, origin="lower", aspect="auto", cmap="viridis")
    ax[0].set_xticks(range(len(KF))); ax[0].set_xticklabels([f"{int(k)}" for k in KF])
    ax[0].set_yticks(range(len(CHI))); ax[0].set_yticklabels([f"{c:.1f}" for c in CHI])
    ax[0].set_xlabel("immune strength kf"); ax[0].set_ylabel("TIP-reservoir coupling  chi")
    ax[0].set_title("P(durable control)\n(chi=0 row = de Boer limit)")
    fig.colorbar(im0, ax=ax[0], label="% control")
    im1 = ax[1].imshow(100 * effect, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-40, vmax=40)
    ax[1].set_xticks(range(len(KF))); ax[1].set_xticklabels([f"{int(k)}" for k in KF])
    ax[1].set_yticks(range(len(CHI))); ax[1].set_yticklabels([f"{c:.1f}" for c in CHI])
    ax[1].set_xlabel("immune strength kf"); ax[1].set_ylabel("coupling chi")
    ax[1].set_title("TIP effect vs de Boer limit (pts)\nred=helps, blue=harms, white=neutral")
    for i in range(len(CHI)):
        for j in range(len(KF)):
            ax[1].text(j, i, f"{100*effect[i,j]:+.0f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im1, ax=ax[1], label="TIP effect (pts)")
    fig.suptitle("P14: coupling phase diagram — the TIP flips from neutral (chi=0, de Boer) to helpful (chi->1)")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "p14_coupling_phase.png"), dpi=130)
    print("\nwrote p14_coupling_phase.npz, p14_coupling_phase.png")


if __name__ == "__main__":
    main()
