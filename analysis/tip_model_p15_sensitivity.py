#!/usr/bin/env python3
"""P15 — global sensitivity: does the coupled-TIP benefit survive JOINT random variation, and
what drives it? (Pre-work: kills the 'you cherry-picked parameters' objection.)

For each Monte-Carlo parameter draw (all key knobs varied together over plausible ranges), we
compute the coupled-TIP effect on durable post-treatment control:
    effect = P(control | chi=1)  -  P(control | chi=0)
and the de Boer-limit control P(control | chi=0). We then report the DISTRIBUTION of the effect
(what fraction of parameter space helps / is neutral / harms) and a parameter-importance ranking
(rank correlation of each parameter with the effect). Run memory-capped.

Varied (plausible/illustrative ranges): immune strength kf, dual-cell visibility nu, TIP
mobilization psi, latency fraction f_lat (reservoir clock), antigen floor RDEF, reactivation
rate A_REACT. Coupling fixed at chi=1 (a well-coupled TIP) -- the question is whether, GIVEN good
coupling, the benefit is robust across everything else.
"""
import os
import numpy as np
import tip_model_p11_coupled as M
from tip_model_p11_coupled import simulate_coupled

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    N, NS = 110, 30
    rng = np.random.default_rng(20260621)
    # plausible ranges
    R = dict(kf=(9.0, 14.0), nu=(0.1, 1.0), psi=(10.0, 100.0),
             f_lat=(2e-5, 2e-4), rdef=(2000.0, 10000.0), areact=(5e-4, 2e-3))
    names = list(R)
    samp = {k: rng.uniform(*R[k], NS) for k in names}
    effect = np.zeros(NS); base = np.zeros(NS)
    print(f"P15 global sensitivity: {NS} Monte-Carlo draws (all params varied jointly), N={N}.")
    print("effect = coupled-TIP (chi=1) minus de Boer-limit (chi=0) control.\n")
    for s in range(NS):
        M.A_REACT = float(samp["areact"][s])            # reactivation rate (module global)
        kw = dict(rdef=float(samp["rdef"][s]), nu=float(samp["nu"][s]),
                  psi=float(samp["psi"][s]), f_lat=float(samp["f_lat"][s]), t_ati=400.0)
        kf = float(samp["kf"][s])
        c0 = simulate_coupled(N, kf, chi=0.0, seed=s, **kw)[0]
        c1 = simulate_coupled(N, kf, chi=1.0, seed=s, **kw)[0]
        base[s] = c0; effect[s] = c1 - c0
    M.A_REACT = 1e-3                                     # restore default

    pe = 100 * effect
    helps = (pe >= 5).mean(); neut = ((pe > -5) & (pe < 5)).mean(); harms = (pe <= -5).mean()
    print(f"TIP-effect distribution over parameter space (n={NS}):")
    print(f"  HELPS (>=+5 pts):   {helps:.0%}    mean effect {pe.mean():+.1f} pts")
    print(f"  ~neutral (-5..+5):  {neut:.0%}")
    print(f"  HARMS (<=-5 pts):   {harms:.0%}    min effect {pe.min():+.0f} pts (worst case)")
    print(f"  => no-backfire across the joint space: {(pe > -5).mean():.0%} of draws")

    # parameter importance: Spearman rank corr of each param with the effect
    def spearman(x, y):
        rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
        return float(np.corrcoef(rx, ry)[0, 1])
    print("\nparameter importance (Spearman rank corr with the TIP effect):")
    imp = sorted(((spearman(samp[k], effect), k) for k in names), key=lambda t: -abs(t[0]))
    for r, k in imp:
        print(f"  {k:8s}  rho = {r:+.2f}")

    np.savez(os.path.join(HERE, "p15_sensitivity.npz"),
             **{k: samp[k] for k in names}, effect=effect, base=base)
    print("\nwrote p15_sensitivity.npz")
    print("Reading: if HELPS+neutral ~ 100% and HARMS ~ 0%, the coupled-TIP benefit is robust to")
    print("joint parameter variation, not a tuned corner; the importance ranking says what drives it.")


if __name__ == "__main__":
    main()
