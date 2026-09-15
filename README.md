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
| `surrogate.json` | emulator hyperparameters, PCA basis and held-out skill |
| `surrogate.bin` | emulator weights, `float64` (see *The emulator*) |

## The emulator

The explorer has a second mode. **Solved grid** shows the 2401 simulations;
**Emulator prediction** replaces the seven-step sliders with continuous ones
and draws what a Gaussian-process emulator predicts for any parameter
combination in range, with a band at its held-out surface error. Nothing in
that mode was solved, and it is drawn as a dashed surface rather than as
particles because a placement plus a shape is all the emulator outputs.

**It has never seen this grid.** It is trained on the 3420 space-filling
(Latin hypercube) simulations of the accompanying paper, so the faint solved
deposit drawn underneath a prediction is a genuine out-of-sample comparison,
not a fit being shown against its own training data.

Inference is a dot product against the stored training inputs, which is why
the whole model is ~270 KB rather than megabytes:

```
z*(x) = [ k(x, X) . alpha ] * y_std + y_mean
k(a,b) = const * (1 + sqrt(3) d) exp(-sqrt(3) d),  d^2 = sum_k ((a_k-b_k)/ls_k)^2
```

`surrogate.bin` is `float64`: first `n_train * n_feat` values are the
standardised training inputs `X` row-major, then `alpha` for each output in
the order given by `outputs` in the header. Six outputs — `x_min`, `spread`,
and four shape-PCA coefficients — reconstruct the profile through
`pca_mean` + `pca_components`, on the same 40 bins as `grid7.json`. The
`skill` block holds the measured held-out R² and RMSE per output and the
surface RMSE the uncertainty band is drawn at. (`float64` rather than
`float32` because the fitted noise level sits at its lower bound, so the
kernel matrix is near-singular and `alpha` is large with heavy cancellation;
in `float32` that cost 3.5e-4 of the output spread on the leading shape mode.)

### What it cannot do

**The starting position is not predictable, and the page does not pretend
otherwise.** `x_min` is bimodal: 97.8% of the solved deposits have their back
edge at the slope toe (sd 1.9 mm), while ~2% detach and start anywhere out to
1.55 m — and those few carry 98% of its variance. Which case detaches is not
a function of the four parameters, so the emulator scores **R² below zero** on
it. The page therefore anchors every prediction at the toe and says so; the
predicted *length and shape* are meaningful where the predicted *position*
would not be. This is the second of the two failure classes named in the
accompanying paper.

Two further caveats. The analysis in the paper uses random forests for
placement; here placement is a GP too, because a 500-tree forest exports to
10–20 MB, and the paper's own comparison shows GP and forest agreeing to
within 0.03 R² on every descriptor. And the band is the emulator's held-out
error, **not** a GP posterior interval: the posterior variance needs the
Cholesky factor of a 3420×3420 matrix (~23 MB), and it would measure
interpolation between design points rather than error against the solver.

### `particles.bin`

All 1,320,684 particles, as little-endian `int16` pairs — x then z, in
millimetres — with the cases laid end to end. Each case in `grid7.json` carries
`po` (offset, in particles) and `pn` (count), so its slice is

```js
const P = new Int16Array(await (await fetch("particles.bin")).arrayBuffer());
for (let i = 0; i < c.pn; i++) {
  const o = (c.po + i) * 2;
  const x = P[o] / 1000, z = P[o + 1] / 1000;   // metres
}
```

**The first `c.pin` of them are the deposit** — the ones every measurement is
taken on. The remaining `c.pn - c.pin` are material still on the chute
(43,009 particles across the grid) or past the end of the plate (8,255), kept
so a viewer can draw them for context. Material that fell below the plate is
not in the file at all.

Membership is decided when the file is written, on the unrounded coordinates,
and not recoverable from the stored ones: rounding to the millimetre moves
particles across the x = 1.99 m boundary, so re-deriving it here would
disagree with the measurements in about 2% of cells.

Counts run from 125 to 682 particles, so the offsets are not a fixed stride.
Millimetres in `int16` loses nothing: the particle spacing in the model is
10 mm.

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
| `po`, `pn` | offset and count of this case's particles in `particles.bin` |
| `pin` | how many of those are the deposit, not context |
| `retained` | fraction of the cluster still on the slope face, excluded |
| `ranoff` | fraction that passed the end of the modelled plate, excluded |
| `censored` | see below |

`d.geometry` carries the chute face and the floor plate as **two separate
line segments**, each `[[x0, z0], [x1, z1]]` in metres. They are different
structures and they do not meet — the chute's lowest point is at
(−0.0037, +0.0246) and the floor begins at (−0.0800, −0.0200) — so joining
them into one polyline draws a spurious kink at the toe. The chute is exactly
straight, `z = −1.7320x + 0.0182`, i.e. 60.00°, so two endpoints are an exact
representation rather than a sampling of it.

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
