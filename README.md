# Anti-Helmholtz fields and toroidal vector potential

Calculate the cylindrical field `(Br, Bz)` of an ideal circular current loop
or a pair of opposite-current loops. The example traces a profile perpendicular
to the field in the `(r, z)` plane and revolves it into `ufo.stl`.

The original research goal is an analogous **vector potential A from opposed
toroids**, followed by investigation of a proposed source of `div(A)`. A new
thin-flux-ring model now computes this analogue in Coulomb gauge, with tests
of divergence, curl, and linked-loop circulation. In standard Maxwell theory,
`div(A)` is gauge-dependent; the calculation does not establish an additional
physical source. See the [derivation, gauge control, and next steps](docs/toroidal-potential.md).

```sh
python -m toroidal_potential
# After installation:
antihemholtz-potential --flux 1e-6 --radius 0.1 --offset 0.05
```

The [cavity audit](docs/cavity-audit.md) reconstructs the separate 2016 wxMaxima
profile, fixes its elliptic modulus/parameter mismatch, and computes actual
Maxwell eigenmodes for an explicitly closed comparison geometry. The new
axisymmetric TM solver is checked against a cylinder's analytic Bessel
spectrum and mesh refinement. It does not include a toroidal source, feed
ports, material loading or losses. Its frequencies do not establish a physical
`div(A)` source.

```sh
python -m cavity_modes
antihemholtz-cavity --profile cylinder --outer-radius 0.1
# Rebuild the committed numerical audit; add --plot after installing .[report]:
python scripts/cavity_audit.py
```

The [driven toroid/pickup baseline](docs/toroid-transfer.md) now predicts
loaded pickup voltage from a Thevenin generator, using the finite rectangular
cross-section flux integral and a passive coupled-coil circuit. It includes
winding resistances, source and probe impedance, backreaction and an explicit
capacitive feedthrough. Opposed-pair, known-mismatch and polarity controls
show ordinary signals that can survive a magnetic null. Component values
are assumed examples; actual discrete windings/lead fields and hardware
measurements remain to be supplied.

```sh
python -m toroid_transfer
antihemholtz-transfer --scenario mismatched-pair
python scripts/transfer_audit.py --plot  # requires .[report]
```

The [scalar source/detector gate](docs/scalar-source-gate.md) derives a zero
additional G-voltage prediction for a complete conserved source, ordinary
detector and retarded zero-scalar data in the specified frameworks. Numerical
controls show how omitted endpoint charge or an omitted scalar-potential
coupling can create a false effect. The audit pins the inspected UBT sources;
it does not add a scalar term or claim a finalized UBT action.

```sh
python scripts/scalar_source_audit.py
```

## Finite-core A and UBT wave matching

The [finite-core extension and UBT bridge](docs/ubt-wave-bridge.md) integrates
Coulomb-gauge A over the existing rectangular toroidal core. It checks the
current-to-flux normalization, linked circulation, quadrature refinement and
traceless central A-gradient. This is an external magnetostatic calculation;
internal core points are explicitly rejected.

The same note develops conditional harmonic/polarity tests for an action-derived
physical Theta mode coupled through F². Such a mode and its coupling remain
undetermined. Ordinary magnetic-energy harmonics are computed as controls;
they are not predicted extra voltages. The existing conserved-source G null
result remains unchanged under its stated assumptions.

```sh
python scripts/finite_core_audit.py --check
```

## Install and run

Use Python 3.10 or newer:

```sh
python -m venv .venv
# Linux/macOS; on Windows use .venv\Scripts\activate
source .venv/bin/activate
python -m pip install -e ".[test]"
python main.py
python main.py surface.stl --steps-vertical 80 --steps-horizontal 180
# The installed command provides the same options:
antihemholtz-stl surface.stl
python -m pytest
```

Running without arguments writes `ufo.stl` in the current directory. Importing
`main` does not generate or overwrite a file. NumPy, SciPy, and numpy-stl are
declared dependencies; no plotting backend is required.

## Field coordinates and units

```python
from field_calculator import AntiHemholtz, LoopField

loop = LoopField(B0=1.0, loop_radius=1.0, axial_offset=0.5)
assert loop.Bz(0, 0.5) == 1.0
br, bz = loop.vector(r=0.2, z=0.8)

pair = AntiHemholtz(B0=1.0, loop_radius=1.0, axial_offset=0.5)
assert pair.vector(0, 0) == (0.0, 0.0)
```

- Supply **global**, finite scalar coordinates, with cylindrical `r >= 0`.
  Every loop applies its own `axial_offset` once, including `Baxial(z)` and
  the dimensionless helper methods.
- `loop_radius` must be positive. All lengths use the same unit.
- `B0` is the signed axial field at the center of one loop. For a single-turn
  vacuum loop in SI units, `B0 = mu0 * I / (2 * loop_radius)`; fields then have
  units of tesla. Negative `B0` reverses the current.
- The pair puts `+B0` at `z = +axial_offset` and `-B0` at
  `z = -axial_offset`. This parameter is half the coil separation, not the
  full separation. The field at the pair's center is zero.
- The ideal zero-thickness wire has a singular field. Evaluation on a
  current-carrying wire raises `ValueError`. A normal direction at a field
  null also raises `ValueError` instead of emitting NaNs.
- Historical names (`AntiHemholtz`, `gama`, `K`, `E`) remain available.

The off-axis components use complete elliptic integrals. The first-kind
integral is evaluated using its complementary parameter to retain precision
near the wire; a local axis expansion avoids cancellation in the radial field
near the axis. The tests compare field values against an independent direct
Biot–Savart quadrature, on-axis formulas, and pair symmetries.

## Surface geometry

`generate_surface()` returns `(vertices, faces)` arrays. It traces **one**
meridional profile with fixed-length Euler steps, then reuses that profile at
every azimuth. The angular seam wraps without duplicating a slice. With `V`
vertical segments and `H` horizontal slices, there are `H * (V + 1)` vertices
and `2 * H * V` triangles. Defaults give 2,520 vertices and 4,800 triangles.

The STL is an **open surface** with two boundary rings. It has no wall thickness
or end caps and is not a watertight printable solid. The trace is a numerical
approximation to a curve perpendicular to the field, not a magnetic field
line; increase `--steps-vertical` to reduce Euler integration error. Changing
the horizontal resolution refines only the angular tessellation.

