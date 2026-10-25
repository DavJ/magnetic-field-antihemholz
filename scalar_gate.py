"""Numerical controls for scalar-G source and ordinary-detector conditions.

All examples use dimensionless natural units c=mu0=epsilon0=1. They check
calculus and current completeness, not SI hardware performance or UBT
dynamics. The retarded path is a prescribed uniform current on a wavy loop
or its open half; endpoint charge is essential for the open path.
See docs/scalar-source-gate.md for the analytical statements and limitations.
"""

from math import pi, sqrt

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erf


def _point(position):
    position = np.asarray(position, dtype=float)
    if position.shape != (3,) or not np.all(np.isfinite(position)):
        raise ValueError("position must have three finite Cartesian coordinates")
    return position


def _derivative(function, value, step):
    return (-function(value + 2 * step) + 8 * function(value + step)
            - 8 * function(value - step) + function(value - 2 * step)) / (12 * step)


def diagnose_potentials(potential, time, position, step=1e-3):
    """Differentiate independent phi,A evaluations: E,B and G=phi_t+div(A)."""
    position, time, step = _point(position), float(time), float(step)
    if not np.isfinite(time) or not np.isfinite(step) or step <= 0:
        raise ValueError("time must be finite and step must be finite and positive")
    phi_gradient = np.empty(3)
    a_gradient = np.empty((3, 3))
    for axis in range(3):
        direction = np.eye(3)[axis]
        phi_gradient[axis] = _derivative(lambda s: potential(time, position + s * direction)[0], 0, step)
        a_gradient[axis] = _derivative(lambda s: potential(time, position + s * direction)[1], 0, step)
    phi_t = _derivative(lambda t: potential(t, position)[0], time, step)
    a_t = _derivative(lambda t: potential(t, position)[1], time, step)
    electric = -phi_gradient - a_t
    magnetic = np.array([
        a_gradient[1, 2] - a_gradient[2, 1],
        a_gradient[2, 0] - a_gradient[0, 2],
        a_gradient[0, 1] - a_gradient[1, 0],
    ])
    spatial_divergence = float(np.trace(a_gradient))
    return {
        "E": electric.tolist(), "B": magnetic.tolist(),
        "div_A": spatial_divergence, "phi_t": float(phi_t),
        "G": float(phi_t + spatial_divergence),
    }


class LocalizedGaugePulse:
    """A=grad(chi), phi=-chi_t for chi=exp(-|x|^2/2-(t-0.4)^2/(2*0.7^2))."""

    def chi(self, time, position):
        position = _point(position)
        return float(np.exp(-np.dot(position, position) / 2 - (time - 0.4) ** 2 / (2 * 0.7 ** 2)))

    def potential(self, time, position):
        position = _point(position)
        chi = self.chi(time, position)
        return (time - 0.4) / 0.7 ** 2 * chi, -position * chi

    def analytic_g(self, time, position):
        position = _point(position)
        return float((np.dot(position, position) - 3 + 1 / 0.7 ** 2
                      - (time - 0.4) ** 2 / 0.7 ** 4) * self.chi(time, position))


def _curve(theta):
    return np.stack((np.cos(theta), np.sin(theta), 0.2 * np.sin(2 * theta)), axis=-1)


def _tangent(theta):
    return np.stack((-np.sin(theta), np.cos(theta), 0.4 * np.cos(2 * theta)), axis=-1)


def retarded_path_potential(time, position, *, closed=True, include_endpoint_charge=True, quadrature_order=128):
    """Retarded potential of a prescribed I(t)=exp(-t^2) on a wavy path.

    Closed: theta in [0,2*pi], rho=0. Open: theta in [0,pi], with
    rho=-Q(t)*delta_start+Q(t)*delta_end and Q'=I. Disabling endpoint
    charge deliberately makes the open source nonconserved. This is not a
    propagation model for a real wire or an anomalous-current experiment.
    """
    position = _point(position)
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if isinstance(quadrature_order, (bool, np.bool_)) or not isinstance(quadrature_order, (int, np.integer)) or quadrature_order < 8:
        raise ValueError("quadrature_order must be an integer >= 8")
    upper = 2 * pi if closed else pi
    nodes, weights = leggauss(quadrature_order)
    theta = upper * (nodes + 1) / 2
    distance = np.linalg.norm(position - _curve(theta), axis=1)
    if np.any(distance < 1e-10):
        raise ValueError("evaluation must avoid the ideal filament")
    factor = upper / (8 * pi) * weights * np.exp(-(time - distance) ** 2) / distance
    vector = np.sum(factor[:, None] * _tangent(theta), axis=0)
    phi = 0.0
    if not closed and include_endpoint_charge:
        endpoint_distance = np.linalg.norm(position - _curve(np.array([0, upper])), axis=1)
        if np.any(endpoint_distance < 1e-10):
            raise ValueError("evaluation must avoid endpoint point charges")
        q = sqrt(pi) / 2 * (1 + erf(time - endpoint_distance))
        phi = float((q[1] / endpoint_distance[1] - q[0] / endpoint_distance[0]) / (4 * pi))
    return phi, vector


def analytic_open_divergence(time, position):
    """Boundary term of div(A): [I(t-Rstart)/Rstart-I(t-Rend)/Rend]/(4*pi)."""
    position = _point(position)
    distance = np.linalg.norm(position - _curve(np.array([0, pi])), axis=1)
    if np.any(distance == 0):
        raise ValueError("evaluation must avoid endpoints")
    values = np.exp(-(time - distance) ** 2) / distance
    return float((values[0] - values[1]) / (4 * pi))


def probe_gauge_integral(order=80):
    """A localized conserved rho,J probe coupled to a time-dependent pure gauge.

    In 1+1 dimensions, per unit transverse area: P=exp(-x^2-t^2),
    rho=-P_x and J=P_t, so rho_t+J_x=0 exactly. Independent quadrature
    integrates rho*chi_t and J*chi_x; their sum is delta(S_int).
    The closed form uses whole-space Gaussian moments. Finite integration
    bounds [-8,8]^2 are used only for the numerical control.
    """
    if isinstance(order, (bool, np.bool_)) or not isinstance(order, (int, np.integer)) or order < 8:
        raise ValueError("order must be an integer >= 8")
    nodes, weights = leggauss(order)
    x, t = np.meshgrid(8 * nodes, 8 * nodes, indexing="ij")
    quadrature_weights = 64 * weights[:, None] * weights[None, :]
    a, b, x0, t0 = 1.3, 0.8, 0.7, 0.4
    p = np.exp(-x ** 2 - t ** 2)
    chi = np.exp(-a * (x - x0) ** 2 - b * (t - t0) ** 2)
    rho, current = 2 * x * p, -2 * t * p
    chi_t, chi_x = -2 * b * (t - t0) * chi, -2 * a * (x - x0) * chi
    scalar_part = float(np.sum(quadrature_weights * rho * chi_t))
    vector_part = float(np.sum(quadrature_weights * current * chi_x))
    normalization = pi / sqrt((1 + a) * (1 + b)) * np.exp(
        -a * x0 ** 2 / (1 + a) - b * t0 ** 2 / (1 + b)
    )
    exact_scalar_part = float(4 * a * b * x0 * t0 / ((1 + a) * (1 + b)) * normalization)
    return {
        "rho_chi_t_integral": scalar_part, "J_chi_x_integral": vector_part,
        "complete_interaction_change": scalar_part + vector_part,
        "whole_space_exact_scalar_part": exact_scalar_part,
        "whole_space_exact_vector_part": -exact_scalar_part,
        "omitting_scalar_potential_gives_false_change": vector_part,
    }
