# Does the proposed scalar G couple to gravity?

Status: **a conditional Maxwell obstruction is established analytically;
a physical UBT scalar-to-metric response is not yet derived.** This is a
theory/calculation audit, not a measurement or a universal exclusion of
additional scalar fields.

The calculation is reproducible with:

```sh
python scripts/gravity_audit.py
python -m pytest tests/test_gravity_gate.py
```

The committed [numerical report](results/scalar-gravity-audit.json) records
versions, constants, assumptions, source commits and blob hashes. It contains
dimensionless curved-background controls and a separate SI circuit calculation.

## 1. Definitions and the exact Maxwell result

Use signature `(-,+,+,+)`, coordinates `x^0=ct`, and the covector
`A_mu=(-phi/c,A_x,A_y,A_z)`. Denote the electromagnetic scalar by

\[
\mathcal G:=\nabla_\mu A^\mu
=\frac1{\sqrt{-g}}\partial_\mu(\sqrt{-g}\,g^{\mu\nu}A_\nu).
\]

In flat Cartesian coordinates this is the preceding audit's
`G=div(A)+phi_t/c^2`. Newton's constant is `G_N`; the Einstein tensor is
`G_munu`, a different object. All curved numerical controls below set
`c=mu0=1` and use the torsion-free Levi-Civita connection.

For every smooth metric and smooth single-valued gauge function `chi`,

\[
A_\mu' = A_\mu+\partial_\mu\chi,\qquad
F_{\mu\nu}'=F_{\mu\nu},\qquad
\mathcal G'=\mathcal G+\Box_g\chi,
\quad \Box_g=\nabla_\mu\nabla^\mu.
\]

The classical Maxwell action and its Hilbert stress tensor are

\[
S_M=-\frac14\int\sqrt{-g}\,F_{\alpha\beta}F^{\alpha\beta}\,d^4x,
\qquad
T_{\mu\nu}=F_{\mu\alpha}F_\nu{}^\alpha
-\frac14g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}.
\]

Consequently `T[g,A+d chi]=T[g,A]` for **every** metric, not only for a flat
background. The action identity holds before metric variation, so varying
the metric cannot turn the gauge shift into a new source. With unchanged
gravitational initial/boundary data, the same spacetime solution is still
allowed. For `A=d chi`, `F=T=0` exactly: there is no additional Maxwell
curvature source even when `G` is nonzero. This says nothing about curvature
already supplied by other matter or a cosmological constant.

For coupled charged matter, the matter phase must also be transformed.
Changing `A` while holding all gauge-charged matter phases fixed is not a
comparison of gauge-equivalent complete states. Ordinary Aharonov-Bohm
holonomy is also preserved; the calculation does not remove linked flux.

The Maxwell action/Hilbert formula agree with David Tong's primary
[GR notes, section 4.5](https://www.damtp.cam.ac.uk/user/tong/gr/grhtml/S4.html).
This audit adds the gauge implication and the concrete controls below.

At a background with `F_bar=0`, a perturbation `F=epsilon*f` has
`T=epsilon^2*T[f]`. Thus an ordinary electromagnetic wave first sources
gravity at quadratic order in its amplitude. A nonzero background field
allows the usual interference term linear in `f`; a pure gauge has `f=0`
and changes neither contribution. This perturbative observation does not
establish the spectrum of the full UBT theory.

## 2. Curvature does not rescue a pure gauge

Choose the prescribed de Sitter metric

\[
ds^2=-dt^2+a(t)^2d\mathbf x^2,\qquad a=e^{Ht},\qquad R=12H^2,
\]

and `chi=t^2+(x^2+y^2+z^2)/2`. Its covector is
`A_mu=(2t,x,y,z)`, so

\[
\mathcal G=-2-6Ht+3e^{-2Ht},\qquad F_{\mu\nu}=0,\qquad T_{\mu\nu}=0.
\]

At `H=0.4`, `t=0.7`, `(x,y,z)=(0.3,-0.2,0.4)`, this gives
`R=1.92` and `G=-1.96637280845...`. These are dimensionless example values,
not measured cosmological parameters. The nonzero background curvature is
unchanged. The polynomial is a local control; it may be smoothly cut off
outside a neighborhood without changing this local statement.

The code independently differentiates `A_mu` and the metric-volume-weighted
raised vector. Step refinement `0.02 -> 0.01 -> 0.005` reduces the absolute
divergence error approximately sixteenfold at each step. It also adds this
gauge to a nonzero homogeneous comoving magnetic field `F_xy=1.2`: `G`
changes while all sixteen stress components remain the same to roundoff.
The latter control avoids relying only on a zero-field example.

Separate tests compare `T` against electric/magnetic energy, Poynting flux
and all spatial stresses in Minkowski space, verify tracelessness, and check
the metric derivative of the action density in all ten symmetric directions
of a non-diagonal Lorentzian metric. These finite checks validate the code;
the general result follows from the action/gauge identity above. There is
no Lean proof or numerical Einstein backreaction solve in this audit.

## 3. The extended-electromagnetism cosmology is a different hypothesis

Beltran Jimenez and Maroto's
[2009 paper, revised 2010](https://arxiv.org/abs/0903.4672) treats an additional
term proportional to `(div_g A)^2` as physical, rather than eliminating
the scalar as a Maxwell gauge degree of freedom. It retains only restricted
gauge transformations. Their equations (7)-(9), (23)-(24) give an additional
stress contribution and a cosmological-constant-like limit for constant
scalar amplitude. This is an external model, not a derived UBT sector or a
verified laboratory effect. We do not independently certify its stability,
quantization or phenomenology here; the preceding Maxwell pure-gauge
equivalence cannot simply be imposed on this different action.

There is an important source limitation even in that model. Taking the
covariant divergence of its vector equation, with conserved current and
nonzero scalar coefficient, gives

\[
\Box_g\mathcal G=0.
\]

On a smooth globally hyperbolic patch, zero scalar Cauchy data
`G|Sigma=0` and `n^mu nabla_mu G|Sigma=0`, with no incoming boundary data,
give the identically zero classical solution. An expanding metric alone
does not change this classical uniqueness statement. The paper's proposed
quantum excitation involves a quantized additional degree of freedom and
vacuum fluctuations; that is not zero classical Cauchy data. For a de Sitter
Fourier mode the same equation reads

\[
\ddot{\mathcal G}_{\mathbf k}+3H\dot{\mathcal G}_{\mathbf k}
+|\mathbf k|^2e^{-2Ht}\mathcal G_{\mathbf k}=0.
\]

The Cauchy uniqueness input is the normally hyperbolic wave-operator result
in Bar, Ginoux and Pfaffle,
[Theorem 3.2.11](https://arxiv.org/pdf/0806.1036), applied to zero smooth
data and source; zero data satisfy the compact-support hypotheses.

A complete ordinary coil current therefore does not acquire a new forced
classical scalar source merely by being placed in curved spacetime under
these assumptions. Nonzero initial scalar data, quantum-state assumptions,
nonconserved effective sources or a different derived action are separate
physical questions. A genuine scalar excitation in a completed theory
could affect ordinary detectors through its induced metric or E/B fields;
the present calculation does not exclude that possibility.

## 4. Ordinary toroid energy and the source-completeness issue

For the reference continuous air-core toroid, `a=0.09 m`, `b=0.11 m`,
`h=0.01 m`, `N=100`,

\[
B_\varphi(r)=\frac{\mu_0NI}{2\pi r},\qquad
U_B=\int\frac{B^2}{2\mu_0}\,dV
=\frac{\mu_0N^2I^2h}{4\pi}\log\frac ba
=\frac12LI^2.
\]

Direct radial volume quadrature agrees with the inductance expression.
At DC `I=1 A`, `L=4.0134139087 microhenry` and
`U_B=2.0067069544 microjoule`; the **field-energy mass equivalent** is
`U_B/c^2=2.2327626052e-23 kg`.

For the existing loaded reference circuit at `10 kHz`, `1 V` source peak,
`50 ohm` source resistance, peak-current phasors are used. The shared-core
field depends on `N1*I1+N2*I2`, including pickup backreaction. The period
mean energies are

\[
\langle U_M\rangle=\frac14
\left(L_1|I_1|^2+L_2|I_2|^2+2M\operatorname{Re}(I_1I_2^*)\right),
\]

\[
\langle U_E\rangle=\frac14
\left(C_p|V_1|^2+C_L|V_2|^2+C_f|V_1-V_2|^2\right).
\]

Core plus pickup-leakage energy independently reproduces the inductance
quadratic form. The mean modeled EM energy is `3.92648e-10 J`, corresponding
to `4.36880e-27 kg`. A point mass of that size at `1 m` would give
`2.91587e-37 m/s^2`. **This last number is only an energy-scale comparison,
not a predicted gravitational signal from the circuit.** Many displayed
digits permit reproduction; the component and geometric inputs are assumed.
No experimental uncertainty or detector sensitivity is inferred.

Why not identify the comparison directly with apparatus gravity? The static
weak-field Einstein equations, with energy-density components in a nearly
Cartesian orthonormal frame, give

\[
\nabla^2\Phi_{\rm grav}
=\frac{4\pi G_N}{c^2}\left(T_{00}+\sum_iT_{ii}\right).
\]

For a vacuum Maxwell field, `T00=u` and `sum Tii=u`, so the field-only local
source is `2u/c^2`. But the field stress is not conserved separately at the
coil/support boundary. The wires, mechanical supports, source and other
energy stores must complete `T_total`; omitting them does not define a
self-consistent isolated static gravitational source.

For a localized stationary conserved system with negligible boundary flux,
integrating `partial_k(x^i T^{kj})` gives `integral T^{ij} dV=0`. The total
leading weak-field monopole is then `E_total/c^2`; support stresses resolve
the apparent factor of two. An actual externally supplied increment of
total rest energy changes that mass. Transferring energy from an internal
battery into a field does **not** increase the rest energy of an otherwise
isolated apparatus; it redistributes energy and stresses. Heat, leakage and
time-dependent multipoles would also need a complete model.

The weak-field derivation uses the trace-reversed Einstein equations in
[Tong's sections 4.5 and 5.1](https://www.damtp.cam.ac.uk/user/tong/gr/grhtml/S5.html).
The source-completeness conclusion is the conservation/integration argument
just given. Constants are taken from `scipy.constants` and recorded in the
report; the values agree with the
[NIST CODATA 2022 table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt).

## 5. What the inspected UBT sources actually establish

Snapshot: `UBT-Institute/unified-biquaternion-theory` at
`fc4c6834b9f052a21a4054418024764ac046fa6e`. Links below are immutable.
This is a read-only source audit, not a rerun of the repository's GR proofs.

| Source | Inspected result and its scope |
| --- | --- |
| [ACTION.en.md](https://github.com/UBT-Institute/unified-biquaternion-theory/blob/fc4c6834b9f052a21a4054418024764ac046fa6e/canonical/ACTION.en.md) | The unique full microscopic action is not finalized. Sector functionals cannot silently replace it. |
| [GR status registry](https://github.com/UBT-Institute/unified-biquaternion-theory/blob/fc4c6834b9f052a21a4054418024764ac046fa6e/canonical/gr_closure/gr_recovery_status.yaml) and [gravity postulate](https://github.com/UBT-Institute/unified-biquaternion-theory/blob/fc4c6834b9f052a21a4054418024764ac046fa6e/canonical/gr_closure/gravity_dynamical_postulate.en.md) | Local classical Einstein-Lambda recovery is marked closed after adopting the one-coupling split-jet Palatini/HT postulate. The EN/CZ postulate files are present on this master snapshot. Full action, microscopic source reduction, global continuation and UV stability remain separately open. |
| [Canonical metric definition](https://github.com/UBT-Institute/unified-biquaternion-theory/blob/fc4c6834b9f052a21a4054418024764ac046fa6e/canonical/CANONICAL_DEFINITIONS.md) | The metric is the real central tetrad anticommutator, not a componentwise real projection or a proposed electromagnetic divergence. Only this definition is used here. |
| [Effective QED sector](https://github.com/UBT-Institute/unified-biquaternion-theory/blob/fc4c6834b9f052a21a4054418024764ac046fa6e/canonical/interactions/qed.tex) | Declares local U(1) and an antisymmetric Maxwell photon term. Those effective premises obey the gauge argument above. This is not endorsement of all older representation/expansion claims in that file. |
| [Independent-metric Einstein variation](https://github.com/UBT-Institute/unified-biquaternion-theory/blob/fc4c6834b9f052a21a4054418024764ac046fa6e/canonical/t_munu/step3_einstein_with_matter.tex) | Explicitly assumes an independently varied metric and an Einstein-Hilbert term. It does not supply the missing physical G mode or its microscopic coupling. |

These status statements answer different questions: an adopted local classical
gravity branch can be closed while the unified microscopic action and source
reduction remain unfinished. It would be inaccurate to call all UBT gravity
open, or to infer a new G-gravity interaction from the closed vacuum branch.
The postulate's integration constant `Lambda0` is not thereby identified with
the electromagnetic scalar or its square.

To obtain a testable additional response, the **same** completed action must
identify a physical excitation in the existing variables, remove gauge and
constraint directions around a nondegenerate on-shell background, and derive
its source and induced metric/ordinary-detector response. The composite
variation must be retained; one cannot simply impose independent metric and
Theta equations by termwise separation. None of the inspected sources
provides that chain for this toroid's proposed G experiment. This is a missing
derivation, not a proof that every possible UBT completion fails. No new
fundamental field, interaction term or fitted detector gain is introduced here.

## 6. Consequence for an experiment

The current calculation supplies a Maxwell null prediction for a **pure gauge
change** and ordinary EM energy scales. It supplies no competing quantitative
UBT gravitational signal, so a hardware comparison between the two cannot yet
be specified from these files. An experiment that drives a real toroid also
changes real E/B fields, heating, stress and source energy; it is not an
experimental operation that varies only a Maxwell gauge.

A later experiment needs a complete source/readout geometry, independently
characterized ordinary forces and energy transfer, and a derived extra signal
with frequency, phase and amplitude. A residual alone would not identify G
or establish a gravitational origin.
