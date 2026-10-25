# Toroidal vector potential and the proposed div(A) source

The owner's original goal is to reproduce, using toroids, a vector-potential
pattern analogous to the magnetic field of the anti-Helmholtz pair, and then
investigate whether a physical source of `div(A)` can be built.

The first calculation below uses standard Maxwell magnetostatics in free
space. It establishes the analogy and a gauge check. It does not establish
an additional longitudinal electromagnetic source or a hardware design.

## 1. Thin-flux-ring approximation

Represent a toroid by a circular magnetic flux filament of major radius R,
at z = z0, carrying signed cross-sectional flux Phi along +phi. Its small
cross-section has been collapsed to zero area at fixed Phi. Real finite-core
fields cannot be inferred at or near this singular filament.

For a localized divergence-free B, the decaying Coulomb-gauge potential can
be obtained by the inverse-curl integral

```math
\mathbf A_C(\mathbf x)=\frac{1}{4\pi}\int
\mathbf B(\mathbf x')\times
\frac{\mathbf x-\mathbf x'}{|\mathbf x-\mathbf x'|^3}\,d^3x'.
```

Collapsing the toroidal cross-section gives our model's line integral:

```math
\mathbf A_C(\mathbf x)=\frac{\Phi}{4\pi}
\oint_C\frac{d\boldsymbol\ell'\times(\mathbf x-\mathbf x')}
{|\mathbf x-\mathbf x'|^3}.
```

This has the same kernel as the Biot–Savart field of a current loop, with
`mu0*I` replaced by `Phi`. Consequently the existing elliptic-integral loop
solver can be reused with the normalization

```math
B_0\ \text{in the numerical kernel}\quad\longrightarrow\quad
A_0=\frac{\Phi}{2R}.
```

Here A0 has units T m, Phi has units Wb, and R has units m. A0 is the
potential at the center of one isolated thin flux ring in this gauge. It is
not the magnetic field or an ampere-turn setting. Phi is flux through the
toroid's small cross-section, not the N*Phi winding flux linkage.

In particular,

```math
A_z(0,z)=\frac{\Phi R^2}{2[R^2+(z-z_0)^2]^{3/2}},\qquad A_r(0,z)=0.
```

The model has `A_phi = 0`. Off the filament both `div(A_C)` and `curl(A_C)`
vanish. Its B is confined to the singular flux core in the distributional
model. Nevertheless a closed path linking the filament has
`integral A.dl = Phi`, with orientation fixed by the right-hand rule. That
circulation cannot be removed by a smooth, single-valued gauge gradient.
Local zero curl outside the core does not make a linked circulation zero:
every spanning surface for that linked loop encounters the flux core.

## 2. Opposed toroids give a gradient, not a divergence source

Place flux +Phi at z=+d and -Phi at z=-d. Superposition gives exactly the
anti-Helmholtz analogue for A in this idealization and gauge.

Near the center, the leading Cartesian terms are

```math
\mathbf A_C(x,y,z)=g\,(-x/2,-y/2,z)+O(|\mathbf x|^3),\qquad
g=\frac{3\Phi R^2d}{(R^2+d^2)^{5/2}}.
```

There is a nonzero axial derivative, but the full divergence is

```math
\nabla\cdot\mathbf A_C=-g/2-g/2+g=0.
```

In cylindrical coordinates the correct expression is

```math
\nabla\cdot\mathbf A
=\frac{\partial A_r}{\partial r}+\frac{A_r}{r}
+\frac{\partial A_z}{\partial z}.
```

At r=0, smoothness replaces `Ar/r` by `dAr/dr`. Omitting this term or
computing only `dAz/dz` would falsely suggest a source in the middle.

## 3. A direct gauge control

Maxwell potentials obey

```math
\mathbf A'=\mathbf A+\nabla\chi,\qquad
V'=V-\partial_t\chi,
```

so that B and E are unchanged, while

```math
\nabla\cdot\mathbf A'=\nabla\cdot\mathbf A+\nabla^2\chi.
```

The code provides a smooth, time-independent control that also decays at
infinity:

```math
\chi(r,z)=\chi_0e^{-(r^2+z^2)/(2w^2)},\qquad
\nabla^2\chi=\left(\frac{r^2+z^2}{w^4}-\frac{3}{w^2}\right)\chi.
```

It changes `div(A)` while preserving `curl(A)` and every tested closed-loop
circulation. Because this chi is static, V and E are unchanged too. The
Gaussian is a mathematical gauge control, not an extra coil or source.

Lorenz gauge provides another useful distinction:
`div(A) = -(1/c^2) dV/dt`. A nonzero three-dimensional divergence in that
gauge is not, by itself, an additional physical degree of freedom.

## 4. Run the reproducible check

```sh
python -m pip install -e ".[test]"
python -m toroidal_potential
# Or use the installed command; values are SI:
antihemholtz-potential --flux 1e-6 --radius 0.1 --offset 0.05
python -m pytest
```

The JSON report contains numerical divergence and curl outside the cores,
the predicted Gaussian change in divergence, and linked circulations in
both gauges. `diagnose` uses five-point finite differences, with a separate
smooth-axis limit. Its step is a length, and stencils must avoid singular
cores. Small numerical residuals are derivative error, not new sources.

## 5. What remains before a build

The [cavity audit](cavity-audit.md) now supplies a separate, gauge-independent
Maxwell eigenmode calculation for one explicitly closed reference geometry,
with an analytic cylinder control and mesh convergence. It also corrects an
elliptic-parameter error in the old wxMaxima profile. The closed cavity is
empty; its source coupling, loading and losses still need a model.

1. Specify an observable: for example an induced voltage on a defined
   pickup loop, a field measurement, or a predicted force on a specified
   probe. A number for div(A) without a gauge and observable is insufficient.
2. Replace flux filaments with finite toroidal cross-sections and actual
   windings, including leads, permeability, and field leakage. The current
   model takes Phi as input and cannot choose a core, turns, or drive current.
3. If using time-varying flux, add the drive, scalar potential, and dynamic
   fields. The present model is magnetostatic and does not predict radiation,
   losses, bandwidth, or a measurable voltage from a static flux.
4. If proposing a new effect tied to div(A) itself, specify the equations or
   action that makes it observable, its coupling to a detector, and a test
   distinguishing it from ordinary induction and leakage. This must be an
   explicit extension of the current Maxwell model, not a gauge relabeling.

## Sources and derivation status

The flux-filament kernel and center-gradient expressions above are derived
here from the inverse-curl representation and the existing loop solution.
Normalization is checked independently by linked-loop flux integrals. The
gauge and circulation identities can be checked against these primary
teaching sources:

- [Feynman Lectures II, chapter 14, sections 14-1 and 14-2](https://www.feynmanlectures.caltech.edu/II_14.html):
  vector potential, gauge freedom, circulation and steady-current sources.
- [Richard Fitzpatrick, Scalar and Vector Potentials](https://farside.ph.utexas.edu/teaching/jk1/Electromagnetism/node4.html):
  full E/B potential relations and the Lorenz condition.
