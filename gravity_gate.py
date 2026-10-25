"""Maxwell stress/gauge controls and SI energy scales for the toroid audit.

Covariant controls use dimensionless c=mu0=1, signature (-,+,+,+), and
the covector A_mu=(-phi,A_x,A_y,A_z). Thus div_g(A) agrees with the
previous scalar-G definition in flat space. The SI energy conversion is
separate: it is not a gravitational solution for an apparatus or a UBT mode.
"""

from math import isfinite, sqrt

import numpy as np
from scipy.constants import G as NEWTON_G, c


def _metric_data(metric):
    metric = np.asarray(metric, dtype=float)
    if metric.shape != (4, 4) or not np.all(np.isfinite(metric)):
        raise ValueError("metric must be a finite 4x4 matrix")
    if not np.allclose(metric, metric.T, rtol=0, atol=1e-14):
        raise ValueError("metric must be symmetric")
    eigenvalues = np.linalg.eigvalsh(metric)
    if np.count_nonzero(eigenvalues < 0) != 1 or np.any(eigenvalues == 0):
        raise ValueError("metric must be nondegenerate with signature (-,+,+,+)")
    return metric, np.linalg.inv(metric), sqrt(-float(np.linalg.det(metric)))


def maxwell_stress_energy(field_cov, metric_cov):
    """Return covariant Hilbert T_mu,nu for Maxwell theory, in c=mu0=1.

    F is the covariant antisymmetric field tensor. This includes stresses,
    not just energy density. It is a local constitutive calculation and
    does not solve Einstein's equations or verify source completeness.
    """
    field = np.asarray(field_cov, dtype=float)
    if field.shape != (4, 4) or not np.all(np.isfinite(field)):
        raise ValueError("field_cov must be a finite 4x4 matrix")
    if not np.allclose(field, -field.T, rtol=0, atol=1e-14):
        raise ValueError("field_cov must be antisymmetric")
    metric, inverse, _ = _metric_data(metric_cov)
    invariant = float(np.sum(field * (inverse @ field @ inverse)))
    return field @ inverse @ field.T - metric * invariant / 4


def diagnose_covariant_potential(covector, metric, coordinates, step=1e-3):
    """Independently differentiate A_mu and sqrt(-g)*g^mu,nu*A_nu.

    F=dA needs no Christoffel symbols for a Levi-Civita connection.
    G=div_g(A) requires the volume factor. Four-dimensional coordinates
    and finite-difference step are dimensionless; output T is not SI.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.shape != (4,) or not np.all(np.isfinite(coordinates)):
        raise ValueError("coordinates must have four finite components")
    if not isfinite(step) or step <= 0:
        raise ValueError("step must be finite and positive")

    def potential(q):
        vector = np.asarray(covector(q), dtype=float)
        if vector.shape != (4,) or not np.all(np.isfinite(vector)):
            raise ValueError("covector must return four finite components")
        return vector

    def weighted_vector(q):
        _, inverse, volume = _metric_data(metric(q))
        return volume * inverse @ potential(q)

    def derivative(function, axis):
        direction = np.eye(4)[axis] * step
        return (-function(coordinates + 2 * direction)
                + 8 * function(coordinates + direction)
                - 8 * function(coordinates - direction)
                + function(coordinates - 2 * direction)) / (12 * step)

    gradient = np.stack([derivative(potential, axis) for axis in range(4)])
    field = gradient - gradient.T
    local_metric, _, volume = _metric_data(metric(coordinates))
    divergence = sum(derivative(weighted_vector, axis)[axis] for axis in range(4)) / volume
    return {
        "G": float(divergence), "field_cov": field.tolist(),
        "stress_energy_cov": maxwell_stress_energy(field, local_metric).tolist(),
    }


def de_sitter_metric(coordinates, hubble=0.4):
    """Prescribed curved background diag(-1,a^2,a^2,a^2), a=exp(H*t)."""
    scale_squared = np.exp(2 * hubble * coordinates[0])
    return np.diag([-1, scale_squared, scale_squared, scale_squared])


def polynomial_gauge_covector(coordinates):
    """d(chi) for chi=t^2+(x^2+y^2+z^2)/2; a local pure gauge."""
    return np.asarray(coordinates, dtype=float) * np.array([2, 1, 1, 1])


def de_sitter_gauge_divergence(time, hubble=0.4):
    """Analytic div_g(d chi)=-2-6*H*t+3*exp(-2*H*t)."""
    return float(-2 - 6 * hubble * time + 3 * np.exp(-2 * hubble * time))


def energy_equivalent_scales(energy_j, distance_m=1.0):
    """U/c^2 and G_N*U/(c^2*r^2), not an apparatus acceleration prediction.

    The latter is a point-mass comparison for an externally added total
    rest-frame energy U. A field subsystem alone lacks the supports/source
    stress and energy needed for a complete conserved gravitational source.
    """
    energy_j, distance_m = float(energy_j), float(distance_m)
    if not isfinite(energy_j) or energy_j < 0:
        raise ValueError("energy_j must be finite and nonnegative")
    if not isfinite(distance_m) or distance_m <= 0:
        raise ValueError("distance_m must be finite and positive")
    mass = energy_j / c ** 2
    return {
        "energy_j": energy_j, "energy_mass_equivalent_kg": mass,
        "reference_distance_m": distance_m,
        "point_mass_acceleration_scale_m_s2": NEWTON_G * mass / distance_m ** 2,
    }
