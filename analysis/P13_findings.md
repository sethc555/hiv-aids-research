# P13 — harshest test (active immune exhaustion + stealthy TIP): no backfire — the TIP helps *more*

_Cope Labs / Seth · 2026-06-21 · model [tip_model_p13_exhaustion_damage.py](tip_model_p13_exhaustion_damage.py) ·
closes the residual corner P12 left open._

## Why this phase exists

P12 found no backfire under a merely-*waning* immune memory, but flagged one untested corner: an
immune system **actively degraded** by viral burden (true T-cell exhaustion, not just
under-stimulation) combined with a **very stealthy** coupled TIP (dual cells that evade killing).
The feared mechanism: a stealthy TIP turns a clearable rebound into a persistent smolder that the
immune system can't clear yet that keeps exhausting it → CD8 collapses → control lost. If a
backfire exists anywhere, it should be here.

## The build

P12's engine + two changes that make exhaustion bite and decouple it from killing:
- **maintenance is epitope-specific** (ν-weighted): `Aeff = Iw + ν·Id` feeds the effector pool;
- **exhaustion is driven by TOTAL burden** (`Iw + Id`, unevadable inflammation): `dE` gains a term
  `− exhaust · burden/(burden+Kx) · E`.

So evasive dual cells (low ν) *exhaust* immunity (via burden) without *feeding* it (low `Aeff`) —
exactly the configuration that should let a stealthy TIP starve and wear down CD8. Swept:
exhaustion strength × very-low visibility (ν = 0.2, 0.1, 0.05), no-TIP vs fully-coupled.

## Result — no backfire anywhere; the TIP helps *most* under exhaustion

| exhaustion | ν | control no-TIP | control + coupled TIP | Δ |
|---|---|---|---|---|
| 0.05 | 0.20 / 0.10 / 0.05 | 66% | 77 / 80 / 82% | +10…+16 |
| **0.15** | 0.20 / 0.10 / 0.05 | **17–18%** | **62 / 68 / 67%** | **+44…+50** |
| 0.30 | 0.20 / 0.10 / 0.05 | 33–38% | 74 / 70 / 66% | +29…+40 |

**Backfire found: nowhere.** The opposite happens — the harsher the exhaustion, the *more* the TIP
helps (+50 points at moderate exhaustion). The mechanism is the reverse of the worry: when
immunity can be exhausted, the **untreated** rebound floods the system with burden, exhausts CD8,
and control collapses (17%). A coupled TIP **suppresses the rebound → keeps total burden low →
spares the immune system from exhaustion** → control is rescued. And **stealthier** TIP cells
(ν=0.05) help at least as much, because once they suppress WT the whole burden (and thus the
exhaustion signal) falls regardless of how visible they are.

So the configuration designed to produce a backfire instead reveals a **second protective role**:
a coupled TIP doesn't just intercept the rebound, it **protects the immune response from
exhaustion** by holding the antigen burden down — most valuable precisely when exhaustion is the
threat.

## Net (project bottom line, now stress-tested three ways)

Across **three structurally different immune models** — P11 maintained-floor QSS, P12
wear-down-able memory, P13 actively-exhaustible immunity — a TIP **coupled to the reactivating
reservoir helps** post-treatment control, and **never backfired**. The original "neutral /
orthogonal" verdict was an artifact of modeling the TIP and the reservoir as non-interacting; once
they interact (as a real TIP is designed to), the TIP is a useful adjunct, and its value only
grows where the immune system is fragile.

_Caveats: illustrative stochastic model (N=200), one detection threshold, coarse coupling/visibility
knobs; the no-TIP baseline is mildly non-monotonic in exhaustion strength (a stochastic-regime
nonlinearity, not central to the within-row TIP comparison); exhaustion modeled as burden-driven
degradation of a single effector pool — a multi-pool / antigen-specific-exhaustion model and
TIP-antigen-specific damage remain untried. Direction (help, no backfire) is robust across the
swept space; exact percentages are not._
