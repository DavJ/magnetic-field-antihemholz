# Scalar G: source and detector conditions for the toroid experiment

**Česky:** Pro úplný obvod se zachováním náboje, nulová počáteční skalární
data a běžný detektor vychází dodatečný signál připisovaný samostatnému G
**nulový**. V probíraném rozšířeném AB modelu je zdrojem G porušení lokální
rovnice kontinuity, nikoli gradient A nebo tvar toroidu. Numerický příklad
ukazuje, že vynecháním nábojů na koncích vodiče dostaneme zdánlivé G; po jejich
doplnění se členy vyruší. Totéž platí pro falešný účinek čisté kalibrace na
detektor, pokud vynecháme skalární potenciál. Tento podmíněný výsledek
nevylučuje všechny možné skalární módy UBT. Aktuální základní akce UBT však
není finalizovaná a z její efektivní U(1) části nový G signál neplyne.

## 1. Conditional result and its domain

Consider smooth localized sources, or wire/surface distributions interpreted
weakly, on flat spacetime. Include the **complete** charge and current:
winding, return/feed, generator and surface/endpoint charge. Assume

```math
I_s=\partial_t\rho_s+\nabla\cdot\mathbf J_s=0.
```

Use the retarded solution with zero scalar initial data and no incoming
scalar wave. Use an ordinary locally conserved detector, without an added
anomalous coupling. Then the independent scalar contribution to the pickup
voltage is

```math
\delta V_G=0.
```

In Maxwell theory G is not an independent local observable in the first
place. In the particular extended AB theory below, a conserved source with
these initial/boundary data also produces no forced G. This does not set
the ordinary [loaded pickup voltage](toroid-transfer.md) to zero, remove
leakage/capacitive backgrounds, or predict an experimental sensitivity.

The result does not forbid homogeneous scalar initial data, externally
supplied scalar waves, a genuinely nonconserved source in an extension,
or a different action-derived UBT mode. Those are different premises and
need their own source, energy and detector derivations.

## 2. Keep gauge representation distinct from a physical excitation

Our conventions are

```math
\mathbf E=-\nabla\phi-\partial_t\mathbf A,
\quad\mathbf B=\nabla\times\mathbf A,
\quad G=c^{-2}\partial_t\phi+\nabla\cdot\mathbf A,
\quad\Box=c^{-2}\partial_t^2-\nabla^2.
```

For a smooth single-valued gauge function chi,

```math
\mathbf A'=\mathbf A+\nabla\chi,
\quad\phi'=\phi-\partial_t\chi,
\quad G'=G-\Box\chi.
```

E and B remain unchanged. The time-dependent localized gauge control in
`scalar_gate.py` has nonzero G but zero E,B, checked by differentiating the
potentials separately. This is a Maxwell representation of zero field.
An arbitrary chi with Box(chi) nonzero is **not** a symmetry of the extended
AB action, and cannot be invoked as a free physical AB-field generator with
the same source. See [Tong, section 6.1.1](https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S6.html)
for standard gauge redundancy.

## 3. Retardation does not turn conserved winding current into a G source

In the lambda=1 extended AB model used by Minotti and Modanese, the potential
equation and its divergence are

```math
\Box A^\mu=\mu_0J^\mu,
\qquad \Box G=\mu_0 I_s.
```

For the stated causal data, let D_ret be the retarded inverse of Box. For
localized sources, integration by parts in the convolution gives

```math
G_{\rm ret}=\partial_\mu(\mu_0D_{\rm ret}*J^\mu)
=\mu_0D_{\rm ret}*(\partial_\mu J^\mu)
=\mu_0D_{\rm ret}*I_s=0.
```

Commuting the derivative requires the stated decay/boundary assumptions;
for distributions the same identity is weak. The implication also needs
the initial-data restriction: Box(G)=0 alone does not imply G=0.
This is a statement about this specified extension, not a derivation of it
from UBT. See [the paper, equations (3), (21g) and (22)](https://arxiv.org/html/2310.13968v1).

For a closed prescribed line current,

```math
\mathbf J(\mathbf x,t)=I(t)\oint_C d\boldsymbol\ell'\,
\delta^3(\mathbf x-\mathbf r'),\qquad \rho=0,
```

its distributional divergence is the integral of a full path derivative
and vanishes. Equivalently the retarded vector potential has divergence

```math
\nabla\cdot\mathbf A
=-\frac{\mu_0}{4\pi}\oint_C d\left[\frac{I(t-R/c)}{R}\right]=0.
```

This includes the time dependence of the retarded I; it is not merely a
magnetostatic cancellation. A real winding need not have spatially uniform
current: its associated charge must then be included in the general
continuity argument. The prescribed uniform-current example does not model
signal transport or the full wire/drive dynamics of a fabricated toroid.

### Open path: the missing endpoint charges

For a path from a to b, with uniform I(t)=dQ/dt,

```math
\nabla\cdot\mathbf J=I(t)(\delta_a-\delta_b),
\qquad\rho=-Q(t)\delta_a+Q(t)\delta_b.
```

Charge accumulation gives rho_t=-I(t)(delta_a-delta_b). It cancels the
current divergence exactly. In units c=mu0=epsilon0=1 the retarded
potentials therefore obey

```math
\nabla\cdot\mathbf A=\frac1{4\pi}
\left[\frac{I(t-R_a)}{R_a}-\frac{I(t-R_b)}{R_b}\right],
\qquad \partial_t\phi=-\nabla\cdot\mathbf A.
```

For a wavy half-loop and I(t)=exp(-t^2), the numerical control gives:

| Source model | div(A) | phi_t | G |
|---|---:|---:|---:|
| Closed path | about 0 | 0 | about 0 |
| Open path, endpoint charges included | 0.03542176219 | -0.03542176219 | about 0 |
| Same open current, endpoint charges omitted | 0.03542176219 | 0 | 0.03542176219 |

These are **dimensionless calculus examples**, not measured G or voltages.
The last row deliberately changes the source to a nonconserved one; it is
not gauge equivalent to the complete source. Its electric field also
changes. It cannot stand in for a complete ordinary laboratory circuit.

## 4. Ordinary detector: retain the complete coupling

The standard electromagnetic interaction action is

```math
S_{\rm int}=\int(\mathbf J_p\cdot\mathbf A-\rho_p\phi)\,d^3x\,dt.
```

Its variation under the gauge change above is

```math
\Delta S_{\rm int}=\int(\mathbf J_p\cdot\nabla\chi+\rho_p\partial_t\chi)
=-\int\chi(\nabla\cdot\mathbf J_p+\partial_t\rho_p)
+\text{boundary terms}.
```

For conserved probe current, with vanishing boundary terms or consistent
endpoint state phases, there is no gauge-dependent local signal. Dropping
the rho*phi term can manufacture one. A localized 1+1-dimensional example,
per unit transverse area, uses P=exp(-x^2-t^2), rho=-P_x and J=P_t. Its
separately integrated gauge contributions are

```math
\int\rho\,\chi_t=+0.306716309993\ldots,\qquad
\int J\,\chi_x=-0.306716309993\ldots.
```

Their sum is zero. Independent whole-space Gaussian moments check the
finite-domain quadrature. These are action-integral controls, not SI detector
gains. A reduced detector model must preserve the full current identity,
including any charge/connection terms it eliminated.

This does not make conserved probes blind to **all** effects of a physical
scalar extension. If a genuine G source changes E,B, an ordinary probe can
detect those changes. The cancellation here concerns pure-gauge variations.
The paper's proposed pure gauge waves have **G=0**, and their direct detection
uses a hypothetical nonconserved probe. They must not be conflated with
generation of nonzero G.

Nor does the argument remove the standard Aharonov-Bohm effect: a loop linked
with toroidal flux can have nonzero circulation of A even with G=0. Removing
that holonomy would require a multivalued or singular gauge outside the
smooth single-valued class used above.

## 5. Consequence of the inspected effective UBT U(1) sector

The inspection is pinned to UBT commit
`fc4c6834b9f052a21a4054418024764ac046fa6e`:

- [canonical/interactions/qed.tex](https://github.com/UBT-Institute/unified-biquaternion-theory/blob/fc4c6834b9f052a21a4054418024764ac046fa6e/canonical/interactions/qed.tex)
  declares D_mu=partial_mu+i*g*A_mu, local phase symmetry and current
  conservation. Only these effective-sector premises are used here. This
  is not an endorsement of every representation, expansion or prediction
  in that file.
- [canonical/ACTION.en.md](https://github.com/UBT-Institute/unified-biquaternion-theory/blob/fc4c6834b9f052a21a4054418024764ac046fa6e/canonical/ACTION.en.md)
  still records DEFINED_FAMILY_NOT_FINALIZED. A postulated effective QED
  sector is not the missing derivation from the unique single-Theta action.

To make the gauge consequence explicit, write components of the existing
matter field as Theta^a, not new fields. For a differentiable effective action
invariant under delta(Theta^a)=i*alpha*Theta^a and
delta(A_mu)=-(1/g)*partial_mu(alpha), its Euler derivatives satisfy

```math
\frac1g\partial_\mu E_A^\mu
+i\sum_a(E_a\Theta^a-\bar E_a\Theta^{a*})=0.
```

Here g is nonzero, E_A^mu=delta(S)/delta(A_mu), E_a=delta(S)/delta(Theta^a)
and bar(E_a)=delta(S)/delta(Theta^{a*}); conjugate field coordinates are
varied independently before taking the physical real slice.

This follows by varying the action, integrating the gauge-parameter
derivative by parts and using arbitrary compactly supported alpha. On the
matter equations E_a=bar(E_a)=0. For a Maxwell term plus the full matter
current, antisymmetry of F then yields partial_mu(J^mu)=0. In curved
spacetime the corresponding density/covariant continuity identity is used.
This is a conditional Noether identity of the stated effective symmetry,
not a construction of a new UBT action or its physical mode spectrum.
It also does not establish the electromagnetic gauge-reduction map of the
full composite-metric/connection single-Theta theory; that map must itself
be derived rather than inferred from the postulated effective sector.

While that U(1) freedom remains, bare G changes under arbitrary local gauge
functions and cannot be declared a new observable. A different physical
Theta mode and its invariant source/detector couplings would need to follow
from the finalized action and constrained spectrum. No extra field, G^2 term
or fitted detector gain is inserted here.

The [existing scalar-energy audit](https://github.com/DavJ/research/blob/49e78d7a389ee7d3e502a9f036b02dc7a67b8e6b/audits/scalar_g/README.md)
is cited for its separate restricted energy obstruction; that suite was not
rerun here and its ansatz is not promoted to the full UBT theory.

## 6. Verification and next physical criterion

The derivative/convolution, boundary-term and Noether arguments above are
analytical arguments with explicit hypotheses. They are not proved by a
finite numerical grid. No Lean formalization is claimed.

The new numerical tests cover a time-dependent pure gauge, physical E/B
diagnostic signs, a retarded nonplanar loop, an independently known open-path
boundary term, endpoint-charge cancellation and Gaussian detector-action
integrals. The [JSON report](results/scalar-source-gate.json) records the
source snapshots, versions and scope. Samples are in dimensionless units.

```sh
python -m pip install -e ".[test]"
python -m pytest
python scripts/scalar_source_audit.py
```

For the proposed ordinary toroid and ordinary pickup, the present prediction
is the Maxwell transfer already computed, with **no additional independent G
voltage** under the stated causal data. Reopening a nonzero-G experimental
prediction requires deriving a changed source or physical mode and detector
coupling, including acceptable energy and constraints. A residual voltage
alone does not supply that derivation.
