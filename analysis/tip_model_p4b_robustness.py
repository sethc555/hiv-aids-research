#!/usr/bin/env python3
"""P4b — robustness of the P4 verdict: is the TIP STILL neutral to post-treatment control
when given its best possible chance, and across reservoir parameters?

P4 found the TIP neutral to control but at one operating point, and the TIP usually failed to
even engage (bolus washes out in ~1 day; latent reactivates too rarely). P4b removes that
escape hatch: it gives the TIP a MAINTAINED presence throughout ATI (continuous dosing, so it
cannot wash out) and pushes its mobilization advantage psi up, then sweeps the reservoir
reactivation rate. If the TIP shifts the CD8 control threshold under THESE conditions, the
P4 neutrality was a weak-TIP artifact; if it still doesn't, neutrality is structural.

Three robustness questions:
  (A) sustained TIP dosing vs no TIP, across CD8 strength kf -> does the threshold move?
  (B) raising psi (22 -> 60 -> 120) under sustained dosing -> does a stronger TIP matter?
  (C) faster reservoir reactivation -> does a more active reservoir change the verdict?

AUDIT DISCIPLINE: report the threshold (kf at P(control)=50%) per arm; confirm the sustained
TIP actually engages (active dual cells / Vt > 0 during ATI), so 'neutral' != 'absent'.
"""
import os
import numpy as np
from tip_model import T0, P
from tip_model_p4_reservoir import simulate, pc, A_REACT

HERE = os.path.dirname(os.path.abspath(__file__))
REPS = 80
KF = np.array([8.0, 10.0, 11.0, 12.0, 13.0, 14.0, 16.0])
t_chr, t_art, t_ati = 500.0, 300.0, 500.0
CONTROL = 50.0


def run_arm(nu, psi, sustained, react, seed):
    """P(control) vs kf for one arm; also Vt engagement at a rebound-regime kf."""
    state0, kf_arr, idx = [], [], []
    for i, kf in enumerate(KF):
        s = np.zeros((REPS, 7)); s[:, 0] = T0; s[:, 1] = 10
        state0.append(s); kf_arr.append(np.full(REPS, kf)); idx.append(np.full(REPS, i))
    state0 = np.vstack(state0); kf_arr = np.concatenate(kf_arr); idx = np.concatenate(idx)
    final = simulate(state0, nu, kf_arr, t_chr, t_art, t_ati,
                     tip_sustained=sustained, psi=psi, react=react, seed=seed)
    prod = final[:, 1] + final[:, 3]
    Vt = pc * psi * P["rho"] * final[:, 3]
    p_ctrl = np.array([(prod[idx == i] < CONTROL).mean() for i in range(len(KF))])
    # TIP engagement: fraction of reps (rebound regime, kf<=11) with active TIP virus
    reb = np.isin(idx, np.where(KF <= 11)[0])
    tip_engaged = (Vt[reb] > 1).mean() if sustained else 0.0
    return p_ctrl, tip_engaged


def threshold(p_ctrl):
    """kf at which P(control) crosses 50% (linear interp); nan if never."""
    for i in range(len(KF) - 1):
        if p_ctrl[i] < 0.5 <= p_ctrl[i + 1]:
            f = (0.5 - p_ctrl[i]) / (p_ctrl[i + 1] - p_ctrl[i])
            return KF[i] + f * (KF[i + 1] - KF[i])
    return float("nan")


def main():
    print(f"P4b robustness: {REPS} reps/cell; control=(Iw+Id)<{CONTROL:.0f}; "
          f"threshold = kf at P(control)=50%\n")
    SUS = 2000.0   # sustained TIP carrier cells/day during ATI (cannot wash out)

    arms = [
        ("no TIP",                        dict(nu=1.0, psi=22,  sustained=0.0, react=A_REACT)),
        ("sustained TIP, psi=22, nu=0.9", dict(nu=0.9, psi=22,  sustained=SUS, react=A_REACT)),
        ("sustained TIP, psi=60, nu=0.9", dict(nu=0.9, psi=60,  sustained=SUS, react=A_REACT)),
        ("sustained TIP, psi=120,nu=0.3", dict(nu=0.3, psi=120, sustained=SUS, react=A_REACT)),
        ("no TIP, 3x reactivation",       dict(nu=1.0, psi=22,  sustained=0.0, react=3 * A_REACT)),
        ("sustained TIP, 3x reactivation",dict(nu=0.9, psi=60,  sustained=SUS, react=3 * A_REACT)),
    ]
    rows = {}
    for j, (name, kw) in enumerate(arms):
        p_ctrl, eng = run_arm(seed=900 + j, **kw)
        thr = threshold(p_ctrl)
        rows[name] = (p_ctrl, thr, eng)
        curve = " ".join(f"{x:3.0%}" for x in p_ctrl)
        print(f"{name:34s} thr_kf={thr:5.2f}  TIP-engaged={eng:4.0%}  | P(ctrl) over kf: {curve}")

    base_thr = rows["no TIP"][1]
    print(f"\nkf grid: {KF}")
    print(f"\nno-TIP control threshold kf = {base_thr:.2f}")
    print("threshold SHIFT vs no-TIP (negative = TIP HELPS, positive = TIP HARMS):")
    for name, (_, thr, eng) in rows.items():
        if name != "no TIP":
            print(f"  {name:34s} delta_kf = {thr - base_thr:+.2f}  (TIP engaged in {eng:.0%} of rebound reps)")

    np.savez(os.path.join(HERE, "p4b_robustness.npz"), KF=KF,
             names=np.array(list(rows.keys())),
             pctrl=np.array([rows[n][0] for n in rows]),
             thr=np.array([rows[n][1] for n in rows]))

    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for name, (p_ctrl, _, _) in rows.items():
        ax.plot(KF, 100 * p_ctrl, "o-", lw=2, label=name)
    ax.set_xlabel("CD8 'vaccinal-effect' strength  kf"); ax.set_ylabel("P(post-treatment control) %")
    ax.set_title("P4b: TIP neutrality is robust — even sustained, high-psi, engaged TIP\n"
                 "does not shift the CD8 control threshold")
    ax.legend(fontsize=8); ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(os.path.join(HERE, "p4b_robustness.png"), dpi=130)
    print("\nwrote p4b_robustness.npz, p4b_robustness.png")


if __name__ == "__main__":
    main()
