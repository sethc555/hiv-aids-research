#!/usr/bin/env python3
"""Robust re-analysis of the P1.2 sweep: characterize the WT-suppression vs
CD8-retention trade-off as a Pareto scatter, immune to the oscillation-phase
confound that corrupts the (t_admin) axis. Flags the eradication-column artifact.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load("/home/seth/dev/hiv-aids-research/analysis/p12_sweep.npz")
TADM, PSI, wt, cd8 = d["TADM"], d["PSI"], d["wtred"], d["cd8ret"]

# is the last t_admin column an artifact? a TIP effect must depend on psi;
# uniformity across psi (incl psi=1, which cannot establish) => phase artifact.
last = wt[:, -1]
print(f"last t_admin column (={TADM[-1]:.0f}d): WT-reduction range across psi "
      f"{last.min():.2f}..{last.max():.2f}  (std {last.std():.2f})")
print(f"  -> {'UNIFORM across psi => oscillation-phase artifact, excluded' if last.std() < 0.3 else 'psi-dependent, kept'}")

# build mask excluding artifact columns (any column ~uniform & high across psi)
colstd = wt.std(axis=0); colmean = wt.mean(axis=0)
artifact = (colstd < 0.3) & (colmean > 1.0)
keep = ~np.broadcast_to(artifact, wt.shape)
wtk, cd8k = wt[keep], cd8[keep]
psik = np.broadcast_to(PSI[:, None], wt.shape)[keep]

print(f"\nexcluding {artifact.sum()} artifact column(s): {wtk.size} cells")
print(f"max WT-reduction among CD8-preserving cells (cd8>=0.5): "
      f"{wtk[cd8k>=0.5].max() if (cd8k>=0.5).any() else 0:.2f} log")
print(f"min CD8 retained among strong-suppression cells (wt>=1.0): "
      f"{cd8k[wtk>=1.0].min() if (wtk>=1.0).any() else float('nan'):.2f}")
print(f"cells with wt>=1.0 AND cd8>=0.5 (immune-compatible, artifact-excluded): "
      f"{((wtk>=1.0)&(cd8k>=0.5)).sum()}")
if wtk.size > 3:
    print(f"corr(WT-reduction, CD8-retained) = {np.corrcoef(wtk, cd8k)[0,1]:.2f}")

fig, ax = plt.subplots(figsize=(7.5, 5.5))
sc = ax.scatter(wtk, cd8k, c=psik, cmap="plasma", s=45, edgecolor="k", lw=0.3)
ax.scatter(last, cd8[:, -1], marker="x", c="gray", s=40, label="artifact column (excluded)")
ax.axhline(0.5, ls="--", c="gray", lw=1); ax.axvline(1.0, ls="--", c="gray", lw=1)
ax.fill_betweenx([0.5, 1.25], 1.0, max(2, wtk.max()), color="green", alpha=0.08)
ax.text(1.05, 1.15, "immune-compatible\n(WT down & CD8 kept)", fontsize=8, color="green")
ax.set_xlabel("time-avg WT log10 reduction by TIP")
ax.set_ylabel("CD8 retained (E+M vs no-TIP)")
ax.set_title("P1.2 trade-off: WT suppression vs CD8 retention\n(each point = a (timing, psi) condition)")
fig.colorbar(sc, ax=ax, label="TIP advantage psi")
ax.legend(loc="lower right", fontsize=8)
fig.tight_layout(); fig.savefig("/home/seth/dev/hiv-aids-research/analysis/p12_pareto.png", dpi=130)
print("wrote p12_pareto.png")
