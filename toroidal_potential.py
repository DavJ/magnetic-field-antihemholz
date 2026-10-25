"""Coulomb-gauge A of ideal thin toroidal flux rings and gauge diagnostics.

This is a zero-core-size, magnetostatic model, not a finite-winding hardware
simulation. See docs/toroidal-potential.md for the derivation and limitations.
Coordinates are global cylindrical (r, z), A_phi = 0. SI inputs give A in T m.
"""

import argparse
from dataclasses import asdict, dataclass
import json
from math import cos, exp, isfinite, pi, sin

from scipy.integrate import quad

from field_calculator import AntiHemholtz, LoopField


def _finite(value, name):
    value = float(value)
    if not isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _positive(value, name):
    value = _finite(value, name)
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _coordinates(r, z):
    r, z = _finite(r, "r"), _finite(z, "z")
    if r < 0:
        raise ValueError("cylindrical radius r must be non-negative")
    return r, z


class ThinToroidPotential:
    """Ideal flux filament on the circle r=major_radius, z=axial_offset.

    flux is the signed magnetic flux through the toroid's small cross-section
    (webers), along +phi for positive flux; it is not winding flux linkage.
    The vector potential is singular on the filament itself.
    """

    def __init__(self, flux, major_radius, axial_offset=0):
        self.flux = _finite(flux, "flux")
        self.major_radius = _positive(major_radius, "major_radius")
        self.axial_offset = _finite(axial_offset, "axial_offset")
        self._kernel = LoopField(self.flux / (2 * self.major_radius), self.major_radius, self.axial_offset)

    def vector(self, r, z):
        """Return (Ar, Az) in Coulomb gauge; A_phi is zero."""
        return self._kernel.vector(r, z)


class OpposedToroidsPotential:
    """Two identical thin toroids with fluxes +flux and -flux at z=+/-offset."""

    def __init__(self, flux, major_radius, axial_offset):
        self.flux = _finite(flux, "flux")
        self.major_radius = _positive(major_radius, "major_radius")
        self.axial_offset = _positive(axial_offset, "axial_offset")
        self._kernel = AntiHemholtz(self.flux / (2 * self.major_radius), self.major_radius, self.axial_offset)

    def vector(self, r, z):
        """Return (Ar, Az); opposite fluxes produce an A-gradient at the center."""
        return self._kernel.vector(r, z)


class GaussianGauge:
    """Smooth, single-valued, static gauge chi=amplitude*exp(-|x|^2/(2w^2)).

    amplitude has units of webers, width has units of length. This changes the
    potential representation, not the toroid, its flux, or the physical fields.
    chi and its gradient vanish at infinity.
    """

    def __init__(self, amplitude, width):
        self.amplitude = _finite(amplitude, "amplitude")
        self.width = _positive(width, "width")

    def value(self, r, z):
        r, z = _coordinates(r, z)
        return self.amplitude * exp(-0.5 * ((r / self.width) ** 2 + (z / self.width) ** 2))

    def gradient(self, r, z):
        r, z = _coordinates(r, z)
        factor = -self.value(r, z) / self.width ** 2
        return factor * r, factor * z

    def laplacian(self, r, z):
        r, z = _coordinates(r, z)
        # The THREE-dimensional Laplacian includes (1/r) dchi/dr.
        return self.value(r, z) * (((r / self.width) ** 2 + (z / self.width) ** 2 - 3) / self.width ** 2)


class GaugeTransformedPotential:
    """Represent the same static magnetic field with A' = A + grad(chi)."""

    def __init__(self, potential, gauge):
        self.potential = potential
        self.gauge = gauge

    def vector(self, r, z):
        ar, az = self.potential.vector(r, z)
        gr, gz = self.gauge.gradient(r, z)
        return ar + gr, az + gz


@dataclass(frozen=True)
class PotentialDiagnostics:
    """Numerical divergence and (curl A)_phi; both have units of tesla in SI."""

    divergence: float
    curl_phi: float


def _derivative(function, position, step):
    """Centered five-point derivative, for smooth off-filament evaluations."""
    return (
        -function(position + 2 * step) + 8 * function(position + step)
        - 8 * function(position - step) + function(position - 2 * step)
    ) / (12 * step)


def diagnose(potential, r, z, step=1e-4):
    """Differentiate an axisymmetric poloidal A, including the cylindrical term.

    div A = dAr/dr + Ar/r + dAz/dz, and B_phi = dAr/dz - dAz/dr.
    step uses the same length unit as r and z. Stencils must avoid flux cores.
    On the axis, smoothness gives Ar/r -> dAr/dr and dAz/dr -> 0.
    """
    r, z = _coordinates(r, z)
    step = _positive(step, "step")
    daz_dz = _derivative(lambda coordinate: potential.vector(r, coordinate)[1], z, step)
    dar_dz = _derivative(lambda coordinate: potential.vector(r, coordinate)[0], z, step)
    if r == 0:
        dar_dr = (8 * potential.vector(step, z)[0] - potential.vector(2 * step, z)[0]) / (6 * step)
        divergence = 2 * dar_dr + daz_dz
        daz_dr = 0.0
    else:
        radial_step = min(step, r / 3)
        dar_dr = _derivative(lambda coordinate: potential.vector(coordinate, z)[0], r, radial_step)
        daz_dr = _derivative(lambda coordinate: potential.vector(coordinate, z)[1], r, radial_step)
        divergence = dar_dr + potential.vector(r, z)[0] / r + daz_dz
    return PotentialDiagnostics(float(divergence), float(dar_dz - daz_dr))


def circulation(potential, center_r, center_z, path_radius):
    """Integrate A.dl around a meridional circle, oriented with normal +phi.

    A circle linking one positive thin flux core has circulation +flux even
    though curl A vanishes everywhere on the path. The path must avoid cores.
    """
    center_r = _positive(center_r, "center_r")
    center_z = _finite(center_z, "center_z")
    path_radius = _positive(path_radius, "path_radius")
    if path_radius >= center_r:
        raise ValueError("path_radius must be smaller than center_r to stay at positive r")

    def integrand(angle):
        r = center_r + path_radius * cos(angle)
        z = center_z - path_radius * sin(angle)
        ar, az = potential.vector(r, z)
        return -path_radius * (ar * sin(angle) + az * cos(angle))

    return float(quad(integrand, 0, 2 * pi, epsabs=1e-12, epsrel=1e-10)[0])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flux", type=float, default=1e-6, help="signed flux of the upper toroid in Wb")
    parser.add_argument("--radius", type=float, default=0.1, help="major radius in m")
    parser.add_argument("--offset", type=float, default=0.05, help="positive half-separation in m")
    args = parser.parse_args(argv)
    try:
        pair = OpposedToroidsPotential(args.flux, args.radius, args.offset)
        gauge = GaussianGauge(args.flux, args.radius)
        shifted = GaugeTransformedPotential(pair, gauge)
        # Choose a linking circle strictly smaller than the radius and the
        # distance to the other toroid; no finite core radius is implied.
        path_radius = min(args.radius, args.offset) / 3
        points = [(0, 0), (0.5 * args.radius, 0.4 * args.offset), (1.5 * args.radius, 1.6 * args.offset)]
        samples = []
        for r, z in points:
            original = diagnose(pair, r, z, step=1e-4 * args.radius)
            transformed = diagnose(shifted, r, z, step=1e-4 * args.radius)
            samples.append({
                "r_m": r, "z_m": z,
                "A_coulomb_T_m": pair.vector(r, z),
                "coulomb": asdict(original), "gauge_transformed": asdict(transformed),
                "expected_divergence_change_T": gauge.laplacian(r, z),
            })
        report = {
            "model": "ideal thin flux rings; static; Coulomb gauge and a static Gaussian gauge transform",
            "flux_Wb": pair.flux, "major_radius_m": pair.major_radius, "half_separation_m": pair.axial_offset,
            "diagnostic_units": "divergence and curl_phi: T; A: T m; circulation: Wb",
            "linking_circulation_Wb": {
                "coulomb": circulation(pair, args.radius, args.offset, path_radius),
                "gauge_transformed": circulation(shifted, args.radius, args.offset, path_radius),
                "expected": args.flux,
            },
            "samples": samples,
            "interpretation": "Changing div(A) by a smooth gauge does not create a new electromagnetic source.",
        }
    except ValueError as exc:
        parser.error(str(exc))
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
