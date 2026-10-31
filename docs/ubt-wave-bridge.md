# Finite-core A and a conditional UBT wave bridge

**2026-10-31 · AI-assisted working research; exhaustive human review is not claimed.**

**Česky:** Analogie anti-Helmholtzova pole pro A je nyní spočtena i pro konečný
obdélníkový průřez toroidů. Nenulový gradient A má nulovou divergenci; oběhový
integrál kolem jádra odpovídá jeho toku. To je standardní Maxwellův výsledek.
Samostatný fyzikální mód Θ by potřeboval vlastní odvozený zdroj, propagátor a
detektor. Podmíněná vazba přes F² by měla jinou harmonickou a polaritní odezvu
než lineární indukce, ale její existenci ani velikost současná UBT neurčuje.

## 1. From the filament to a rectangular core

The [thin-flux-ring construction](toroidal-potential.md) defines A in the
Coulomb gauge. Replace the flux filament by the already-declared ideal
rectangular toroid in `toroid_transfer.RectangularToroid`. For inner/outer
radii a,b, height h, turns N and static signed current I, assume linear uniform
permeability and continuous ideal windings:

```math
B_\phi(r')=\frac{\mu_0\mu_rNI}{2\pi r'},\quad
\Phi=\frac{\mu_0\mu_rNIh}{2\pi}\log\frac ba.
```

B is confined to the core in this ideal model. For non-unit permeability this
includes the idealized magnetization; discrete turns and real material response
are not resolved. Take the inverse curl of this complete prescribed B, with
A decaying at infinity. If `K_ring(x;r',z')` is the existing thin-ring Coulomb A
per unit flux, linearity gives

```math
\mathbf A_C(\mathbf x)=\int_a^bdr'\int_{z_0-h/2}^{z_0+h/2}dz'\,
B_\phi(r')\,\mathbf K_{\rm ring}(\mathbf x;r',z').
```

There is **no extra 2*pi*r' volume factor** in this flux integral: that factor
is already in the azimuthal line kernel. `dPhi=B_phi*dr'*dz'`, not winding
flux linkage. `finite_toroid.py` evaluates it by tensor Gauss–Legendre
quadrature. The total quadrature flux is checked against the logarithmic formula.

For observation points outside every core the kernel is smooth. The solver
rejects on-core and internal points; a singular-volume quadrature is required
there and is not implemented. Near an exterior boundary, increase quadrature
order and verify convergence; the fixed-order routine is not an adaptive
error estimator. The core boundary carries the ideal source current; a
finite-difference diagnostic must keep its entire stencil outside the core.

On the axis there is the independent continuous check

```math
A_z(0,z)=\frac12\int d\Phi\,
\frac{r'^2}{[r'^2+(z-z')^2]^{3/2}},\qquad A_r(0,z)=0.
```

## 2. Opposed finite toroids retain the A-gradient

Place identical disjoint cores at +/-d with currents +/-I. Require d>h/2.
Reflection/axial symmetry and `div(A)=0`, `curl(A)=0` in the central current-free
region give

```math
\mathbf A_C(x,y,z)=g(-x/2,-y/2,z)+O(|\mathbf x|^3),
\quad g=\left.\partial_z A_z\right|_0,
\quad\nabla\cdot\mathbf A_C=-g/2-g/2+g=0.
```

The coefficient is computed by differentiating the integrated on-axis kernel.
As h and b-a tend to zero at fixed normalized flux, it tends to the earlier
filament formula `g=3*Phi*R^2*d/(R^2+d^2)^(5/2)`. Finite width produces an
actual correction to A and g; it does not turn the traceless gradient into a
new source. The report gives the correction for explicit assumed dimensions,
not an optimal geometry or a measured effect.

Every meridional loop enclosing the **whole** upper cross-section and avoiding
both cores has `integral A.dl=+Phi`; the corresponding lower loop has -Phi.
A smooth single-valued gauge transformation preserves these circulations and
B. Linked holonomy is a physical flux datum despite zero curl on the external
path. Bare values of A and its local gradient remain gauge-dependent.

```python
from finite_toroid import OpposedFiniteToroidsPotential
from toroid_transfer import RectangularToroid
from toroidal_potential import diagnose, circulation

pair = OpposedFiniteToroidsPotential(RectangularToroid(), current_a=1,
                                    axial_offset=.05, order=10)
print(pair.center_gradient)              # T; gradient of Coulomb-gauge A
print(diagnose(pair, .02, .01))           # current-free point: div/curl about zero
print(circulation(pair, .1, .05, .02))    # Wb; encloses full upper cross-section
```

This is a **magnetostatic** solver. Replacing I by I(t) in A alone is not a
complete dynamic model: scalar potential, retardation, drive and losses must be
included consistently. The existing [loaded transfer model](toroid-transfer.md)
remains the ordinary voltage baseline within its declared quasistatic circuit
assumptions. The finite-core A solver does not silently change that circuit.

## 3. A physical Theta mode must survive the constraints

The common derivation is in
[hyperspace_waves/UBT_ANALYSIS.md](https://github.com/DavJ/hyperspace_waves/blob/research/ubt-wave-bridge-20261031/UBT_ANALYSIS.md).
It fixes the complex BQ inverse/exponential and the previously omitted imaginary
mass-shell equation. The former sqrt(2)*c inference fails that equation.
Complex temporal damping is not itself a new propagating degree of freedom.

The UBT snapshot
[`386b6db`](https://github.com/UBT-Institute/unified-biquaternion-theory/blob/386b6dbed7d35ea63c8fd42fd88d57ac2f9cdaa3/canonical/ACTION.en.md)
keeps the fundamental action `DEFINED_FAMILY_NOT_FINALIZED`. Therefore we do
not assign a finalized physical Theta spectrum, add a fundamental G² action,
choose a component count, or replace a constrained pairing by a positive
coefficient norm. These finite-core calculations are conditional Maxwell
geometry; the microscopic UBT source–detector matching remains open.

For a derived reduced physical-mode sector q, a homogeneous effective comparator
would have `K*q_tt-C*Delta(q)+M*q=s`, with positive K and nonnegative C,M.
An explicitly supplied detector d then reads in Fourier space

```math
\delta V=d^T[M+k^2C-\omega^2K-2i\omega\Gamma]^{-1}s,
\quad \Gamma\geq0.
```

The localized experiment requires spatial Green integrals, boundaries and units.
The matrices and overlaps must come from the **same** action and physical
constraint reduction. Unknown is not the same as nonzero. For the already
specified conserved source, ordinary probe and zero incoming scalar data, the
[existing G gate](scalar-source-gate.md) stays closed: additional G voltage is
zero under those premises. That result does not exclude every possible derived
physical mode of Theta.

## 4. Constructive conditional channel test: F², not bare div(A)

One matching target for an existing physical scalar mode q of Theta is an
action-derived coefficient `-Z(q)*F²/(4*mu0)`. This is a standard scalar/photon
effective interaction, **not an established term or new discovery in UBT**.
With `Z=Z0+z_A*q^A+...` and SI `F²=2(B²-E²/c²)`, its source would be

```math
s_A=\frac{z_A}{2}(\epsilon_0E^2-B^2/\mu_0).
```

A reduction might give z_A=0 or no such physical scalar. Until it is derived,
this is a conditional source identity, not a predicted field strength. A
linear electromagnetic detector response also needs a nonzero background/probe
field or a separately derived material coupling. No fitted pickup gain is
inserted. Gauge invariance of an operator does not guarantee a UBT coefficient.

For a linear air-core toroid the ordinary magnetic energy is
`W=L*I(t)^2/2`. Put `I=Idc+Iac*cos(omega*t)`:

```math
W=L(Idc^2/2+Iac^2/4)+L Idc Iac\cos\omega t
+(L Iac^2/4)\cos2\omega t.
```

For material mu_r!=1 the circuit energy `int B*H/2` differs from the vacuum
`int B²/(2*mu0)` in the proposed source. Therefore this exact circuit-energy
harmonic helper is a vacuum-F² amplitude proxy only for **mu_r=1** and
magnetically dominated fields; with materials, derive their effective action.
Time-varying E terms also need their own calculation.

For disjoint equal cores with negligible cross fields, the two B² integrals
add even when signed flux/pickup induction cancels. Their relative current sign
cannot cancel that even scalar source. This does **not** establish an outgoing
Theta wave; its propagation and detection must still be demonstrated.

| Control | Linear induction | Conditional magnetically dominated F² channel |
|---|---|---|
| Simultaneous DC and AC current reversal | AC sign reverses | unchanged |
| Zero DC bias, pure sinusoidal drive | omega component | DC and 2*omega source |
| Reverse DC only with fixed AC phase | AC induction unchanged in linear model | omega source reverses |
| Opposed disjoint equal cores | can cancel in a symmetric pickup | local even sources add |

These distinctions are useful for falsifying a specified model, not for
identifying new physics by parity alone. Heating, magnetostriction, forces,
permeability changes and electronic nonlinearities can produce the same even
harmonics. The energy coefficients in the report are ordinary stored energy,
not watts emitted into a new channel or detector volts.

## 5. Verification and next gate

```sh
python -m pip install -e '.[test]'
python -m pytest
python scripts/finite_core_audit.py --check
```

The [JSON audit](results/finite-core-bridge.json) records source snapshots,
versions, quadrature refinement, finite-width correction, external div/curl,
linked circulation and ordinary magnetic-energy harmonics. Tests also use
independent adaptive double quadrature on axis, a shrinking-core limit and a
smooth gauge control. No test is an experimental validation. Kernel/calculus
arguments remain analytic; **LEAN-PENDING** for a formal proof.

The next smallest physical step is a **constrained Theta Hessian plus its mixed
electromagnetic derivative** on an upstream nondegenerate on-shell background,
with both source and detector overlaps. Completing a geometry calculation or
observing an unexplained residual does not complete that bridge.

Background primary sources: [Feynman II.14, vector potential](https://www.feynmanlectures.caltech.edu/II_14.html),
[Tong, gauge redundancy](https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S6.html),
[Prokopec, scalar/photon effective operators](https://arxiv.org/abs/astro-ph/0106247).
[Research-origin record](research-origin-2026-10-31.json).
