"""Reconstruct the 2016 open-integral profile, including its elliptic bug.

The historical expression is retained for comparison only. Neither profile
is a derivation of a conducting cavity's eigenmodes. All lengths use metres;
the normalized source has flux 1 Wb, which cancels from the profile equation.
"""

from math import isfinite, pi, sqrt

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ellipe, ellipk

from toroidal_potential import ThinToroidPotential


HISTORICAL_SOURCE = (
    "https://github.com/DavJ/research/blob/"
    "9f4e22c83aabe5a8e54ff47c52e58c46b6986de3/A-field-cavity/num-sol.wxm"
)


def historical_az(r, z, radius=1.0):
    """Evaluate the old formula literally: Maxima's parameter is given k.

    Maxima elliptic_kc/elliptic_ec and SciPy ellipk/ellipe both take m=k^2.
    This intentionally reproduces the old modulus/parameter mismatch.
    Only the interior r < radius is supported, avoiding the singular ring.
    """
    r, z, radius = float(r), float(z), float(radius)
    if not all(isfinite(v) for v in (r, z, radius)) or not 0 <= r < radius:
        raise ValueError("require finite z and 0 <= r < positive radius")
    q = z * z + (radius + r) ** 2
    k = 2 * sqrt(r * radius / q)
    coefficient = (radius * radius - z * z - r * r) / (z * z + (r - radius) ** 2)
    return (coefficient * ellipe(k) + ellipk(k)) / (2 * pi * sqrt(q))


def open_integral(az, r, height):
    """Integral of Az along an open vertical segment, not a voltage."""
    value, _ = quad(lambda z: az(r, z), 0, height, epsabs=1e-11, epsrel=1e-10)
    return float(value)


def integral_profile(radii, axis_half_height, coil_radius, *, historical=False):
    """Solve integral_0^Z Az(r,z) dz = integral_0^Z0 Az(0,z) dz.

    Return positive upper heights; mirroring them about z=0 is an explicit
    geometry choice made by the cavity model. The source ring is outside the
    domain since all radii must be strictly below coil_radius.
    """
    radii = np.asarray(radii, dtype=float)
    axis_half_height, coil_radius = float(axis_half_height), float(coil_radius)
    if (
        radii.ndim != 1 or len(radii) < 2 or not np.all(np.isfinite(radii))
        or not np.isfinite(axis_half_height) or axis_half_height <= 0
        or not np.isfinite(coil_radius) or coil_radius <= 0
        or radii[0] != 0 or np.any(np.diff(radii) <= 0) or radii[-1] >= coil_radius
    ):
        raise ValueError("require increasing radii from 0 to < coil_radius and positive finite lengths")
    # Work in dimensionless coordinates; tolerances must not depend on metres.
    source = ThinToroidPotential(1, 1)
    az = historical_az if historical else lambda r, z: source.vector(r, z)[1]
    z0 = axis_half_height / coil_radius
    target = open_integral(az, 0, z0)
    heights = []
    for r in radii / coil_radius:
        upper = max(z0, 1.0)
        for _ in range(30):
            if open_integral(az, r, upper) >= target:
                break
            upper *= 2
        else:
            raise ValueError("could not bracket a positive integral-profile height")
        height = brentq(lambda z: open_integral(az, r, z) - target, 0, upper, xtol=1e-12)
        heights.append(height * coil_radius)
    return np.asarray(heights)
