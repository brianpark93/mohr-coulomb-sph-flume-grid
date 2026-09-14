# Mohr–Coulomb SPH flume deposits: a complete 7×7×7×7 parameter grid

**▶ Interactive explorer:
<https://brianpark93.github.io/mohr-coulomb-sph-flume-grid/>** — move a slider
and the particle deposit for that parameter combination is drawn. Pin a shape
to compare against while you change one parameter.

2401 LS-DYNA SPH simulations of the same dry granular flume collapse, one for
every combination of seven levels of each of the four Mohr–Coulomb parameters.
Published so that anyone can check what the model actually does as the
parameters move, without running the solver.

The point of a **complete** grid is that any three parameters can be held fixed
while the fourth is stepped, and the simulation for that combination exists.
A space-filling design cannot do this: it scatters its points so that no two
share a coordinate, so "hold φ, ψ and c, sweep G" has no answer in it.

## What was varied

Seven levels each, equally spaced on the axis each parameter was sampled on.

| Parameter | Symbol | Levels | Spacing |
|---|---|---|---|
| Elastic shear modulus | `G` | 1.000, 2.076, 4.309, 8.944, 18.566, 38.540, 80.000 MPa | log-equal, ratio 2.0758 |
| Internal friction angle | `phi` | 0.450 → 0.900 rad | linear, step 0.075 |
| Dilation angle | `psi` | 0 → 0.350 rad | linear, step 0.0583 |
| Cohesion | `c` | 0, 500, 1000, 1500, 2000, 2500, 3000 Pa | linear, step 500 |

`G` is spaced logarithmically because the range spans a factor of 80. Linear
spacing would put the first step at 14 MPa and leave the soft end unresolved,
which is where the stiffness–dilation coupling is strongest.

Held fixed: density 1276.3 kg/m³, Poisson's ratio 0.25, and the geometry —
a block of 645 SPH particles at 10 mm spacing released high on a 60° chute,
running out onto a horizontal plate. Plane strain. All 2401 runs reached
normal termination; total 34.3 h on eight concurrent four-thread solves.

## Files

| File | Contents |
|---|---|
| `index.html` | the interactive explorer; no dependencies, no build step |
| `grid7.json` | case metadata and envelopes, 823 KB |
| `particles.bin` | every particle position, 4.8 MB (see below) |
| `grid7_metrics.csv` | one row per case, scalars only, for a spreadsheet or pandas |

### `particles.bin`

All 1,269,420 deposit particles, as little-endian `int16` pairs — x then z, in
millimetres — with the cases laid end to end. Each case in `grid7.json` carries
`po` (offset, in particles) and `pn` (count), so its slice is

```js
const P = new Int16Array(await (await fetch("particles.bin")).arrayBuffer());
for (let i = 0; i < c.pn; i++) {
  const o = (c.po + i) * 2;
  const x = P[o] / 1000, z = P[o + 1] / 1000;   // metres
}
```

Deposit sizes run from 125 to 643 particles, so the offsets are not a fixed
stride. Millimetres in `int16` loses nothing: the particle spacing in the model
is 10 mm.

## Using `grid7.json`

Cases are a **flat array** indexed by the four level numbers, so four slider
positions give the case directly with no lookup table and no search:

```js
const d = await (await fetch("grid7.json")).json();
const n = d.n_levels;                       // 7
const idx = ((gi * n + pj) * n + sk) * n + cl;   // each index 0..6
const c = d.cases[idx];
```

`d.levels` holds the actual parameter value at each index, for slider labels.

Each case stores a **placement** and a **normalised shape**, which is how the
deposit is described throughout the study — keeping *where it sits* separate
from *what form it takes*. Reconstruct the profile with:

```js
const N = d.n_bins;                          // 40
const x = Array.from({length: N}, (_, i) => c.xmin + (i / (N - 1)) * c.s);
const z = c.z.map(v => v / 1000);            // stored in mm
```

and fill down to `d.geometry.floor_z` to draw the deposit body.

Per case:

| Field | Meaning |
|---|---|
| `xmin`, `s` | back edge (m) and footprint length (m) |
| `xmax`, `zmax` | runout (m) and peak deposit thickness (m) |
| `z` | 40-point upper envelope in mm, on ξ ∈ [0, 1] along the deposit |
| `n` | particles in the deposit |
| `retained` | fraction of the cluster still on the slope face, excluded |
| `ranoff` | fraction that passed the end of the modelled plate, excluded |
| `censored` | see below |

`d.geometry` carries the slide face and the floor plate as **two separate
line segments**. They are different structures and they do not meet — the
slide's lowest point is at (−0.0037, +0.0246) and the floor begins at
(−0.0800, −0.0200) — so joining them into one polyline draws a spurious kink
at the toe.

## The `censored` flag — read this before trusting a runout

The modelled floor plate ends at x = 1.99 m. Material that passes it has no
floor beneath it and falls freely, so its position is set by how long it has
been falling rather than by the material parameters. Those particles are
excluded from the deposit, which means that for such a case the **runout is a
lower bound, not a measurement**.

**317 of the 2401 cells (13.2%) are affected.** They are not scattered: they
concentrate at low friction, high dilation, low stiffness and low cohesion —
the most mobile corner of the box. Anything displaying this data should mark
them rather than drawing them as ordinary deposits.

For scale, the laboratory deposit this model was built against spans only
x = 0.007 to 0.708 m, so the plate is already 2.8× longer than the experiment
ever needed. A parameter set that runs off it is one that does not describe
this apparatus.

## What this dataset is not for

**Training a surrogate.** Grid points sit on a lattice and are strongly
correlated, which is the wrong structure for fitting. Checked rather than
assumed: adding all 2401 cases to a 2736-case space-filling training set —
an 88% increase — moved held-out R² by between −0.005 and +0.019 and the
shape reconstruction error from 2.32 to 2.27 cm. Essentially nothing.

That null result is itself informative: the error that remains is not a
shortage of data.

The grid is, however, a clean **independent test set**. A surrogate fitted to
a separate 3420-case Latin hypercube and never shown a grid point predicts all
2401 of them at R² = 0.84 (runout), 0.78 (thickness), 0.89 (footprint), close
to its own cross-validated figures — so that accuracy is a property of the
parameter-to-deposit map and not of one sampling scheme.

## How it was produced

LS-DYNA SPH with `*MAT_MOHR_COULOMB` (MAT_173), elastic–plastic with a
non-associated flow rule. One solve per grid cell, material parameters written
into the keyword deck and nothing else changed. Each deposit is then taken as
the largest connected particle cluster, clipped to 0 ≤ x ≤ 1.99 m, and reduced
to an upper envelope on 40 equal-width streamwise bins.

The solver deck and the analysis code are not included here yet; they belong
with the manuscript this grid was produced alongside, and will be released
with it. Ask if you need them sooner.

## Caveats worth knowing

- The dilation angle is **fixed** for the whole run, with no cap and no decay
  toward a critical state. At low stiffness and high dilation the material
  expands without bound: in the softest, most dilatant corner the particle
  spacing nearly doubles, which no real sand does. Those cells are a property
  of fixed-ψ perfect plasticity, and are best read as a calibration hazard
  rather than as a prediction about sand.
- Deposits are read at t = 1.5 s. At low stiffness with active dilation the
  flow is still creeping slightly at that point; integrating to 5 s changes
  peak thickness by under 4% but runout by up to 8%.
- Peak thickness is the most reliable measure here. It is vertical, so the
  finite plate cannot censor it, and it is the least sensitive to the
  termination time.
