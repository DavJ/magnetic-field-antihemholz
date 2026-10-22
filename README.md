# Anti-Helmholtz magnetic field and STL surface

Calculate the cylindrical field `(Br, Bz)` of an ideal circular current loop
or a pair of opposite-current loops. The example traces a profile perpendicular
to the field in the `(r, z)` plane and revolves it into `ufo.stl`.

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
