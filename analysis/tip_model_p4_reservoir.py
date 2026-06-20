#!/usr/bin/env python3
"""P4 — the latent replication-competent reservoir: the model<->clinic closure.

Every phase carried the same caveat (AUDIT.md, HIGH): the model's "control" was extinction
of ACTIVE infection, but the clinical cure target is the LATENT replication-competent
reservoir, which the model never contained. So "TIPs antagonize the 2024-26 immune-control
breakthroughs" was generalizing an active-cell result to a reservoir mechanism not modeled.

P4 fixes that. It adds two biologically distinct reservoirs to the P3 stochastic model:
  - L_def : the defective/Simonetti antigen clone (from P2/P3) -- presents antigen, primes
            CD8, produces no virus, TIP-proof. (unchanged)
  - L_lat : a LATENT, replication-competent provirus pool -- transcriptionally SILENT (NOT in
            the CD8 antigen pool), long-lived, reactivates stochastically to a productive WT
            cell (Iw). This is the reservoir that defines cure.
  - L_latd: a latent DUAL pool -- an engineered TIP that also integrates/persists, reactivating
            to a dual (Id) cell. This lets a TIP ride through ART and intercept rebound.

The experiment is the clinical one (post-treatment control / ATI):
  (1) CHRONIC  : normal infection establishes setpoint + both reservoirs (TIP arm also seeds
                 a latent dual pool, as if engineered TIP were given before ART).
  (2) ART      : block new infection (b->~0); active infection decays, reservoirs persist.
  (3) ATI      : release ART; L_lat reactivates and attempts rebound. Does CD8 (the "vaccinal
                 effect", knob kf) control it -> DURABLE REMISSION, or does it rebound?

Headline question: across CD8 strength kf, does a persistent TIP RAISE or LOWER
P(post-treatment control)? -- and does it matter whether the TIP's dual cells evade CD8 (nu)?

AUDIT DISCIPLINE: stochastic, report P(control) over replicates; confirm the reservoir
persists in "controlled" replicates (functional, NOT sterilizing, cure).
"""
import os
import numpy as np
from tip_model import P, T0
from tip_model_p13_wm import QS, Estar

HERE = os.path.dirname(os.path.abspath(__file__))
pc = P["p"] / P["c"]
PSI = 22.0
F_LAT = 5e-4        # fraction of new infections entering latency
A_REACT = 1e-3      # latent reactivation rate (/cell/day)
DL = 1e-3           # latent death rate
PL = 1e-3           # latent homeostatic proliferation (balances death -> reservoir maintained)
G = 0.02            # defective antigen-clone timescale
RDEF = 5000.0       # defective antigen-clone size
S_SEED = 1e-4
ART_FACTOR = 1e-3   # residual infectivity on ART


def _binsplit(rng, n, f):
    """Split integer event counts n into (latent, productive) by per-event prob f."""
    n = n.astype(np.int64)
    lat = rng.binomial(np.maximum(n, 0), f)
    return lat, n - lat


def simulate(state0, nu, kf, t_chronic, t_art, t_ati, tip_bolus_It=0.0, dt=0.04, seed=0):
    """Vectorized tau-leaping with CHRONIC -> ART -> ATI schedule. Arrays are per-replicate.
    tip_bolus_It: TIP-carrier cells injected at ATI onset (TIP given at treatment interruption)."""
    rng = np.random.default_rng(seed)
    T, Iw, It, Id, Ldef, Llat, Llatd = (state0[:, k].astype(np.float64).copy() for k in range(7))
    lam, dT, d, b, rho = P["lam"], P["dT"], P["d"], P["b"], P["rho"]
    k = QS["k"]
    n_chr, n_art, n_ati = (int(x / dt) for x in (t_chronic, t_art, t_ati))
    for n in range(n_chr + n_art + n_ati):
        if n == n_chr + n_art and tip_bolus_It:           # TIP bolus at ATI onset
            It = It + tip_bolus_It
        art = 1.0 if n < n_chr else (ART_FACTOR if n < n_chr + n_art else 1.0)
        bt = b * art
        Vw = pc * (Iw + (1 - rho) * Id)
        Vt = pc * (PSI * rho * Id)
        E = Estar(Iw + nu * Id + Ldef)            # L_lat / L_latd are SILENT (not primed against)
        kw, kd = kf * k * E, nu * kf * k * E

        def pois(rate):
            return rng.poisson(np.clip(rate, 0, None) * dt)

        Tprod = pois(np.full_like(T, lam)); Tdeath = pois(dT * T)
        infW = pois(bt * T * Vw); infTIP = pois(bt * T * Vt)
        IwSup = np.minimum(pois(bt * Iw * Vt), Iw.astype(np.int64))
        ItSup = np.minimum(pois(bt * It * Vw), It.astype(np.int64))
        IwD = pois((d + kw) * Iw); ItD = pois(d * It); IdD = pois((d + kd) * Id)
        Ldef_b = pois(G * Ldef + S_SEED * (Iw + Id)); Ldef_d = pois(G * Ldef * Ldef / RDEF)
        Llat_r = np.minimum(pois(A_REACT * Llat), Llat.astype(np.int64))
        Llat_d = pois(DL * Llat); Llat_p = pois(PL * Llat)
        Llatd_r = np.minimum(pois(A_REACT * Llatd), Llatd.astype(np.int64))
        Llatd_d = pois(DL * Llatd); Llatd_p = pois(PL * Llatd)

        # latency splits on the freshly infected
        infW_lat, infW_prod = _binsplit(rng, infW, F_LAT)
        dual_new = IwSup + ItSup
        dual_lat, dual_prod = _binsplit(rng, dual_new, F_LAT)

        T = np.clip(T + Tprod - Tdeath - infW - infTIP, 0, None)
        Iw = np.clip(Iw + infW_prod - IwD - IwSup + Llat_r, 0, None)
        It = np.clip(It + infTIP - ItD - ItSup, 0, None)
        Id = np.clip(Id + dual_prod - IdD + Llatd_r, 0, None)
        Ldef = np.clip(Ldef + Ldef_b - Ldef_d, 0, None)
        Llat = np.clip(Llat + infW_lat - Llat_r - Llat_d + Llat_p, 0, None)
        Llatd = np.clip(Llatd + dual_lat - Llatd_r - Llatd_d + Llatd_p, 0, None)

    return np.stack([T, Iw, It, Id, Ldef, Llat, Llatd], axis=1)


def main():
    KF = np.array([8.0, 10.0, 12.0, 14.0, 16.0, 18.0])   # CD8 strength, bracketing the control transition
    # arm: (label, latent_TIP_seed, tip_bolus_It_at_ATI, nu)
    ARMS = [("no TIP", False, 0.0, 1.0),
            ("TIP latent/engineered (nu=0.9)", True, 0.0, 0.9),
            ("TIP bolus@ATI visible (nu=0.9)", False, 1000.0, 0.9),
            ("TIP bolus@ATI evasive (nu=0.3)", False, 1000.0, 0.3)]
    REPS = 80
    t_chronic, t_art, t_ati = 500.0, 300.0, 500.0
    CONTROL = 50.0                                    # productive cells (Iw+Id) below this = controlled

    print(f"P4 ATI experiment: chronic {t_chronic}d -> ART {t_art}d -> ATI {t_ati}d; "
          f"{REPS} reps/cell; control = (Iw+Id)<{CONTROL:.0f} at end")
    results = {}
    for arm, latent, bolus, nu in ARMS:
        # flatten all KF into one vectorized call
        state0, kf_arr, idx = [], [], []
        for i, kf in enumerate(KF):
            s = np.zeros((REPS, 7)); s[:, 0] = T0; s[:, 1] = 10        # T, Iw seed
            if latent:
                s[:, 3] = 10; s[:, 6] = 50                            # Id + latent dual (engineered TIP)
            state0.append(s); kf_arr.append(np.full(REPS, kf)); idx.append(np.full(REPS, i))
        state0 = np.vstack(state0); kf_arr = np.concatenate(kf_arr); idx = np.concatenate(idx)
        final = simulate(state0, nu, kf_arr, t_chronic, t_art, t_ati,
                         tip_bolus_It=bolus, seed=4040 + int(nu * 100) + int(bolus))
        prod = final[:, 1] + final[:, 3]
        resv = final[:, 5] + final[:, 6]
        ctrl = prod < CONTROL
        p_ctrl = np.array([ctrl[idx == i].mean() for i in range(len(KF))])
        resv_in_ctrl = np.array([resv[(idx == i) & ctrl].mean() if (ctrl & (idx == i)).any()
                                 else np.nan for i in range(len(KF))])
        results[arm] = (p_ctrl, resv_in_ctrl)
        print(f"\n{arm}:")
        for i, kf in enumerate(KF):
            print(f"  kf={kf:.1f}: P(post-treatment control)={p_ctrl[i]:.0%}  "
                  f"(reservoir L_lat+L_latd retained in controlled reps ~{resv_in_ctrl[i]:.0f} cells)")

    np.savez(os.path.join(HERE, "p4_ati.npz"), KF=KF,
             **{f"pctrl_{j}": results[a][0] for j, (a, _, _, _) in enumerate(ARMS)},
             arms=np.array([a for a, _, _, _ in ARMS]))

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    for arm, _, _, _ in ARMS:
        ax.plot(KF, 100 * results[arm][0], "o-", lw=2, label=arm)
    ax.set_xlabel("CD8 'vaccinal-effect' strength  kf  (x baseline killing)")
    ax.set_ylabel("P(post-treatment control) after ATI  (%)")
    ax.set_title("P4: does a persistent TIP help or hurt post-treatment control?\n"
                 "(latent replication-competent reservoir; rebound at ATI)")
    ax.legend(); ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(os.path.join(HERE, "p4_ati.png"), dpi=130)
    print("\nwrote p4_ati.npz, p4_ati.png")


if __name__ == "__main__":
    main()
