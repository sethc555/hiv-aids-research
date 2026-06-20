# P4b/P4c — robustness of TIP-neutrality, and a clinical-timing calibration

_Cope Labs / Seth · 2026-06-20 · models [tip_model_p4b_robustness.py](tip_model_p4b_robustness.py),
[tip_model_p4c_calibration.py](tip_model_p4c_calibration.py) · figure [p4b_robustness.png](p4b_robustness.png)._

## P4b — is the TIP still neutral when given its best possible chance?

P4 found the TIP neutral to post-treatment control, but at one operating point and with the
TIP usually failing to even engage (bolus washes out in ~1 day; latent reactivates too
rarely). P4b removes that escape hatch: **maintained TIP dosing throughout ATI** (so it cannot
wash out), **raised ψ** (22→60→120), and **3× reservoir reactivation**. Threshold = the CD8
strength `kf` at which P(control) crosses 50%.

| arm | threshold kf | Δ vs no-TIP | TIP engaged (rebound reps) |
|---|---|---|---|
| no TIP | 12.08 | — | 0% |
| sustained TIP, ψ=22, ν=0.9 | 11.71 | **−0.37** | 61% |
| sustained TIP, ψ=60, ν=0.9 | 11.75 | **−0.33** | 62% |
| no TIP, 3× reactivation | 11.76 | −0.32 | 0% |
| sustained TIP, ψ=60, 3× react | 11.43 | **−0.65** | 61% |

The TIP is now **genuinely active** (engaged in ~61% of rebound replicates), so this is
"neutral," not "absent." Yet the threshold shifts are tiny (**−0.3 to −0.65 kf**, a ~3–5%
move) and **comparable to the measurement noise** — the no-TIP baseline itself moved −0.32 kf
just from a parameter/seed change. So **for CD8-mediated post-treatment control the TIP is
robustly neutral** (a slight *help* if anything, never the harm a naive reading of the
antagonism would predict), across sustained dosing, ψ up to 60, and 3× reactivation.

### The one regime where a TIP changes the outcome — and why it isn't a cure

The **ψ=120, ν=0.3** arm (strong + immune-evading + continuously dosed) gives P(control)=82–86%
**even at sub-curative CD8** (kf=8, where CD8 alone rebounds 100%). But this is **not** immune
control: it is the **P1–P3 antagonism regime reappearing as TIP-dependent suppression** — the
evasive TIP (ν=0.3 lets its dual cells dodge CD8) out-competes each reactivating WT burst, held
down by the continuous dosing. It requires (i) the **risky immune-evasion property**, (ii)
**continuous dosing** (a "second ART" — stop it and the persisting reservoir rebounds), and it
is **not durable immune remission**. So it doesn't move the verdict about the *cure strategy*;
it's the TIP acting as a chronic suppressive therapy on a different axis.

## P4c — calibration: the rebound clock is clinically realistic *without tuning*

Median time from ATI onset to detectable WT (Vw>1000/mL), 120 reps, vs corpus anchors
(macaque controls ~7.5 d, human ATI ~2–3 wk, bNAb-treated ~64 d / RIO):

| CD8 strength kf | rebound in | median t_rebound | IQR |
|---|---|---|---|
| 6 | 100% | 3.3 d | 2–7 |
| 8 | 100% | 4.4 d | 2–8 |
| 10 | 100% | 7.4 d | 3–15 |
| 12 | 100% | 9.8 d | 5–20 |
| 14 | 100% | 21.2 d | 8–44 |

This lands **on the clinical scale with no fitting**: kf≈10 reproduces the macaque-control
~7.5 d; kf≈12–14 reproduces human ATI ~2–3 wk; and raising CD8 further **delays then abolishes**
rebound — exactly the bNAb/RIO "delayed-or-controlled" pattern. That the reactivation rate
(A_REACT=1e-3/day, chosen on biological grounds, not fitted) yields a realistic rebound clock
lends real credibility to P4's qualitative conclusions.

## Net

The P4 verdict is **robust**: for the immune-control cure, a TIP is neutral even when sustained,
strong, engaged, and run against a 3×-faster reservoir (Δthreshold within noise). The only way
a TIP suppresses rebound is the strong-evasive-continuously-dosed regime, which is TIP-dependent
chronic suppression — not the durable immune remission the field is pursuing. And the model's
rebound timing is clinically realistic unfitted, so these are not artifacts of an off-scale
clock. **The reactivation rate is the first knob to calibrate to a specific ATI cohort** for any
quantitative (vs qualitative) claim.

_Caveats: illustrative params throughout; threshold estimates carry ~±0.3 kf noise at 80 reps;
"control" still = productive-cell suppression with the reservoir persisting; calibration is a
ballpark match to summary statistics, not a fit to patient-level rebound curves._
