"""Finite rectangular-core Coulomb A, for external magnetostatic points.

Continuous ideal windings, linear uniform permeability, B_phi proportional to
1/r inside each core and zero outside. No discrete turns, leads, losses,
radiation or extra div(A) source. See docs/ubt-wave-bridge.md.
"""
from math import isfinite, pi
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import mu_0
from field_calculator import LoopField
from toroid_transfer import RectangularToroid


class FiniteToroidPotential:
    """Integrate dPhi=B_phi(r') dr' dz' against the thin flux-ring kernel.

    Do not evaluate on/inside the core: this external solver does not resolve
    the integrable singular kernel there. current_a is a static signed current.
    """
    def __init__(self, geometry=None, current_a=1.0, axial_offset=0.0, order=10):
        self.geometry = RectangularToroid() if geometry is None else geometry
        if not isinstance(self.geometry, RectangularToroid):
            raise ValueError("geometry must be a RectangularToroid")
        self.current_a, self.axial_offset = float(current_a), float(axial_offset)
        if not isfinite(self.current_a) or not isfinite(self.axial_offset):
            raise ValueError("current_a and axial_offset must be finite")
        if isinstance(order, (bool, np.bool_)) or not isinstance(order, (int, np.integer)) or not 2 <= order <= 64:
            raise ValueError("order must be an integer from 2 to 64")
        self.order = int(order)
        g = self.geometry
        nodes, weights = leggauss(order)
        radii = (g.inner_radius_m+g.outer_radius_m)/2 + nodes*(g.outer_radius_m-g.inner_radius_m)/2
        heights = self.axial_offset + nodes*g.axial_height_m/2
        wr = weights*(g.outer_radius_m-g.inner_radius_m)/2
        wz = weights*g.axial_height_m/2
        coefficient = mu_0*g.relative_permeability*g.drive_turns*self.current_a/(2*pi)
        self._rings = []
        self.quadrature_flux_wb = 0.0
        for radius, radial_weight in zip(radii, wr):
            for height, axial_weight in zip(heights, wz):
                flux = coefficient/radius*radial_weight*axial_weight
                self._rings.append((radius, height, flux, LoopField(flux/(2*radius), radius, height)))
                self.quadrature_flux_wb += flux

    @property
    def flux_wb(self):
        """Analytic core flux, not winding flux linkage."""
        return self.geometry.flux_per_ampere_wb_a*self.current_a

    def contains(self, r, z):
        g = self.geometry
        return g.inner_radius_m <= r <= g.outer_radius_m and abs(z-self.axial_offset) <= g.axial_height_m/2

    def vector(self, r, z):
        r, z = float(r), float(z)
        if not isfinite(r) or not isfinite(z) or r < 0:
            raise ValueError("coordinates must be finite, with r >= 0")
        if self.contains(r, z):
            raise ValueError("external potential solver cannot evaluate on or inside a core")
        values = np.array([kernel.vector(r, z) for _, _, _, kernel in self._rings])
        return tuple(np.sum(values, axis=0))

    def axis_gradient(self, z):
        """Independent analytic derivative of the integrated on-axis kernel."""
        z = float(z)
        if not isfinite(z):
            raise ValueError("z must be finite")
        return sum(-3*flux*radius**2*(z-height)/(2*(radius**2+(z-height)**2)**2.5)
                   for radius, height, flux, _ in self._rings)


class OpposedFiniteToroidsPotential:
    """Identical disjoint rectangular cores, +I above and -I below."""
    def __init__(self, geometry=None, current_a=1.0, axial_offset=0.05, order=10):
        g = RectangularToroid() if geometry is None else geometry
        if not isinstance(g, RectangularToroid):
            raise ValueError("geometry must be a RectangularToroid")
        self.axial_offset = float(axial_offset)
        if not isfinite(self.axial_offset) or self.axial_offset <= g.axial_height_m/2:
            raise ValueError("half-separation must exceed half the core height")
        self.upper = FiniteToroidPotential(g, current_a, self.axial_offset, order)
        self.lower = FiniteToroidPotential(g, -current_a, -self.axial_offset, order)

    def vector(self, r, z):
        return tuple(np.asarray(self.upper.vector(r, z)) + np.asarray(self.lower.vector(r, z)))

    @property
    def center_gradient(self):
        """g in A=g*(-x/2,-y/2,z)+O(|x|^3), in tesla."""
        return self.upper.axis_gradient(0) + self.lower.axis_gradient(0)


def magnetic_energy_harmonics(inductance_h, dc_current_a, ac_peak_a):
    """Coefficients of W=L*(Idc+Iac*cos(wt))^2/2, in joules.

    Returns DC, cos(wt), cos(2wt). W is an ordinary magnetic energy baseline,
    NOT a scalar source normalization or an extra detector voltage.
    Valid for a linear lossless inductance; pair/cross terms belong in L.
    """
    inductance, dc, ac = map(float, (inductance_h, dc_current_a, ac_peak_a))
    if not all(isfinite(v) for v in (inductance, dc, ac)) or inductance < 0:
        raise ValueError("parameters must be finite and inductance nonnegative")
    return {"dc_j": inductance*(dc*dc/2+ac*ac/4),
            "omega_cos_j": inductance*dc*ac,
            "two_omega_cos_j": inductance*ac*ac/4}
