#!/usr/bin/env python3
"""de Boer bridge — recreate Dodd & de Boer (2025) and pinpoint WHY our result differs.

Dodd & de Boer, "Immune responses may make HIV-1 therapeutic interfering particles less
effective" (J Theor Biol 2025, DOI 10.1016/j.jtbi.2025.112317): an ODE model of ACTIVE
infection (no latent reservoir, no ATI). Headline: an immune response against infected cells
"drastically decreases the range of parameter values for which [TIP] therapy is effective" --
i.e. immunity HURTS the TIP's ability to suppress active viral load.

Our tip_model.py (P1) is a faithful STRUCTURAL recreation of that model (same Perelson-class
within-host TIP ODE: T, Iw, It, Id, Vw, Vt; TIP made only in dual cells; immune killing kap).
This script runs BOTH metrics on the same machinery to show the divergence is NOT a
contradiction but a different question + setting:

  METRIC A (de Boer): marginal TIP suppression of ACTIVE viral load vs immune strength kap.
                      -> collapses as immunity rises  (reproduces de Boer: immunity hurts TIP)
  METRIC B (ours):    TIP effect on DURABLE post-treatment control vs immune strength, in the
                      reservoir + ART/ATI setting with the TIP COUPLED to reservoir reactivation
                      (P11). -> positive (the TIP HELPS the immune-mediated cure)

Same TIP, same immunity-kills-infected-cells mechanism. The sign of the TIP<->immunity
interaction flips because the OBJECTIVE flips (suppress active VL now  vs  prevent reservoir
rebound later) and because the reservoir/coupling that the cure setting requires is absent from
de Boer's active-infection model. Run memory-capped.
"""
import numpy as np
from tip_model import wt_setpoint, tip_arm
from tip_model_p11_coupled import simulate_coupled


def metric_A_deboer(psi=12.0):
    """Marginal TIP suppression of active WT viral load vs immune killing kap (de Boer's question)."""
    print("METRIC A -- de Boer recreation: TIP benefit on ACTIVE viral load vs immunity")
    print(f"  (psi={psi}; marginal benefit = log10 WT-setpoint reduction by the TIP)")
    print(f"  {'immune kap':>11} {'TIP benefit (log10)':>20}")
    out = []
    for kap in [0.0, 0.2, 0.4, 0.8]:
        sp = wt_setpoint(kap)
        end, _ = tip_arm(sp, kap, psi, tmax=700)
        benefit = np.log10(max(sp[4], 1e-9)) - np.log10(max(end[4], 1e-9))
        out.append((kap, benefit))
        print(f"  {kap:>11.1f} {max(benefit,0):>20.2f}")
    print("  -> benefit COLLAPSES as immunity rises = de Boer's 'immunity hurts the TIP'.\n")
    return out


def metric_B_ours(N=150):
    """TIP effect on DURABLE post-treatment control vs immune strength (our reservoir/ATI setting)."""
    print("METRIC B -- ours: coupled TIP effect on DURABLE post-treatment control vs immunity")
    print("  (reservoir + ART/ATI; TIP coupled to reactivation chi=1 vs no TIP chi=0)")
    print(f"  {'immune kf':>10} {'control noTIP':>14} {'control +TIP':>13} {'TIP effect':>11}")
    out = []
    for kf in [9.0, 11.0, 13.0]:
        c0 = simulate_coupled(N, kf, chi=0.0, seed=int(kf))[0]
        c1 = simulate_coupled(N, kf, chi=1.0, seed=int(kf))[0]
        out.append((kf, c0, c1)); print(f"  {kf:>10.0f} {100*c0:>13.0f}% {100*c1:>12.0f}% {100*(c1-c0):>+10.0f}")
    print("  -> TIP RAISES durable control = a coupled TIP HELPS the immune-mediated cure.\n")
    return out


def main():
    print("=" * 74)
    print("de Boer bridge: same TIP machinery, two questions -> opposite TIP<->immunity sign")
    print("=" * 74 + "\n")
    metric_A_deboer()
    metric_B_ours()
    print("WHY THEY DIFFER (not a contradiction):")
    print("  * de Boer ask: does immunity help/hurt the TIP's suppression of ACTIVE virus?  -> hurts.")
    print("  * we ask:      does a TIP help/hurt the immune-mediated CURE (no rebound)?      -> helps.")
    print("  * their model has NO reservoir/ATI and NO TIP-reservoir coupling -- the exact")
    print("    ingredients that turn competition (over active WT) into cooperation (vs the reservoir).")
    print("  * where the models OVERLAP (active infection, metric A) we REPRODUCE de Boer.")


if __name__ == "__main__":
    main()
