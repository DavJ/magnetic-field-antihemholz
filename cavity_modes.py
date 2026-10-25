"""Maxwell eigenmodes of an empty, closed, axisymmetric PEC cavity.

Only the azimuth-independent sector H = H_phi e_phi, E = (Er,0,Ez) is
solved (TM to z). The regular variable w = H_phi/r is approximated with
bilinear finite elements on the meridian. PEC walls impose the natural
condition d_n(r H_phi)=0. No gauge potential enters this eigenproblem.

This computes frequencies and relative field shapes, not drive coupling,
losses, a loaded resonator, or a scalar div(A) source. See docs/cavity-audit.md.
"""

import argparse
from dataclasses import dataclass
import json
from math import isfinite, pi

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import c
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from scipy.special import jn_zeros

from cavity_profile import integral_profile


def _positive(value, name):
    value = float(value)
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _count(value, name, minimum=1):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


@dataclass
class MeridianMesh:
    """A structured (r,z) mesh with a piecewise-linear upper/lower wall."""

    nodes: np.ndarray
    cells: np.ndarray
    radial_cells: int
    axial_cells: int


@dataclass
class CavityModes:
    mesh: MeridianMesh
    frequencies_hz: np.ndarray
    eigenvalues_m2: np.ndarray
    w_vectors: np.ndarray
    relative_residuals: np.ndarray
    electric_magnetic_energy_ratios: np.ndarray

    def relative_h_phi(self, mode=0):
        """Nodal H_phi, scaled to unit peak magnitude (no absolute amplitude)."""
        h = self.mesh.nodes[:, 0] * self.w_vectors[:, mode]
        pivot = np.argmax(np.abs(h))
        return h / h[pivot]


def meridian_mesh(radii, half_heights, axial_cells):
    """Close z=+/-Z(r) by a cylindrical wall at the last radius.

    The axis is included. This construction is separate from the older open
    STL surface. Positive heights make every quadrilateral non-inverted.
    """
    radii, half_heights = np.asarray(radii, dtype=float), np.asarray(half_heights, dtype=float)
    axial_cells = _count(axial_cells, "axial_cells", 2)
    if (
        radii.ndim != 1 or len(radii) < 3 or half_heights.shape != radii.shape
        or not np.all(np.isfinite(radii)) or not np.all(np.isfinite(half_heights))
        or radii[0] != 0 or np.any(np.diff(radii) <= 0) or np.any(half_heights <= 0)
    ):
        raise ValueError("require >=3 increasing radii from zero and matching positive finite half-heights")
    nodes = np.stack((
        np.broadcast_to(radii[:, None], (len(radii), axial_cells + 1)),
        half_heights[:, None] * np.linspace(-1, 1, axial_cells + 1),
    ), axis=-1).reshape(-1, 2)
    lower = np.arange(len(radii) - 1)[:, None] * (axial_cells + 1) + np.arange(axial_cells)
    cells = np.stack((lower, lower + axial_cells + 1, lower + axial_cells + 2, lower + 1), axis=-1)
    return MeridianMesh(nodes, cells.reshape(-1, 4), len(radii) - 1, axial_cells)


def _shape(xi, eta):
    sx = np.array([-1, 1, 1, -1])
    sy = np.array([-1, -1, 1, 1])
    values = (1 + sx * xi) * (1 + sy * eta) / 4
    derivatives = np.stack((sx * (1 + sy * eta), sy * (1 + sx * xi)), axis=-1) / 4
    return values, derivatives


def assemble(mesh, quadrature_order=3):
    """Return symmetric stiffness, mass and length scale for w=H_phi/r.

    Coordinates are normalized by the outer radius. The 2*pi azimuth factor
    cancels from the generalized eigenproblem. Three-point Gauss quadrature
    is used in each reference-coordinate direction.
    """
    scale = float(mesh.nodes[:, 0].max())
    coordinates = mesh.nodes[mesh.cells] / scale
    stiffness = np.zeros((len(mesh.cells), 4, 4))
    mass = np.zeros_like(stiffness)
    points, weights = leggauss(_count(quadrature_order, "quadrature_order", 3))
    for xi, wx in zip(points, weights):
        for eta, wy in zip(points, weights):
            n, dn = _shape(xi, eta)
            jacobian = np.einsum("cai,aj->cij", coordinates, dn)
            determinant = np.linalg.det(jacobian)
            if np.any(determinant <= 0):
                raise ValueError("mesh has a non-positive cell Jacobian")
            gradient = np.einsum("ai,cij->caj", dn, np.linalg.inv(jacobian))
            r = coordinates[:, :, 0] @ n
            br = 2 * n[None, :] + r[:, None] * gradient[:, :, 0]
            bz = r[:, None] * gradient[:, :, 1]
            factor = determinant * wx * wy
            stiffness += (factor * r)[:, None, None] * (
                br[:, :, None] * br[:, None, :] + bz[:, :, None] * bz[:, None, :]
            )
            mass += (factor * r ** 3)[:, None, None] * (n[:, None] * n[None, :])
    rows = np.broadcast_to(mesh.cells[:, :, None], stiffness.shape).ravel()
    cols = np.broadcast_to(mesh.cells[:, None, :], stiffness.shape).ravel()
    shape = (len(mesh.nodes), len(mesh.nodes))
    k = coo_matrix((stiffness.ravel(), (rows, cols)), shape=shape).tocsr()
    m = coo_matrix((mass.ravel(), (rows, cols)), shape=shape).tocsr()
    return k, m, scale


def solve_modes(mesh, count=4, *, quadrature_order=3):
    """Solve K w = (k*scale)^2 M w and retain the lowest positive modes."""
    count = _count(count, "count")
    if count >= len(mesh.nodes):
        raise ValueError("mode count must be smaller than the number of mesh nodes")
    k, m, scale = assemble(mesh, quadrature_order)
    values, vectors = eigsh(
        k, k=count, M=m, sigma=0, which="LM", tol=1e-10,
        v0=np.random.default_rng(0).standard_normal(len(mesh.nodes)),
    )
    order = np.argsort(values)
    values, vectors = values[order], vectors[:, order]
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("non-positive or non-finite cavity eigenvalue")
    kw, mw = k @ vectors, m @ vectors
    residuals = np.linalg.norm(kw - mw * values, axis=0) / (
        np.linalg.norm(kw, axis=0) + values * np.linalg.norm(mw, axis=0)
    )
    ratios = np.sum(vectors * kw, axis=0) / (values * np.sum(vectors * mw, axis=0))
    eigenvalues = values / scale ** 2
    return CavityModes(mesh, c * np.sqrt(eigenvalues) / (2 * pi), eigenvalues, vectors, residuals, ratios)


def cylindrical_tm_modes(radius, length, count=4):
    """Analytic axisymmetric TM_0np frequencies of an empty PEC cylinder."""
    radius, length = _positive(radius, "radius"), _positive(length, "length")
    count = _count(count, "count")
    candidates = []
    for n, root in enumerate(jn_zeros(0, count), start=1):
        for p in range(count):
            frequency = c / (2 * pi) * np.hypot(root / radius, p * pi / length)
            candidates.append({"mode": f"TM_0,{n},{p}", "frequency_hz": float(frequency)})
    return sorted(candidates, key=lambda item: item["frequency_hz"])[:count]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("corrected", "historical", "cylinder"), default="corrected")
    parser.add_argument("--coil-radius", type=float, default=0.1, help="reference flux-ring radius, m")
    parser.add_argument("--outer-radius", type=float, default=0.09, help="explicit cylindrical closing wall radius, m")
    parser.add_argument("--axis-half-height", type=float, default=0.1, help="Z0 on the axis, m; cylinder half-length")
    parser.add_argument("--radial-cells", type=int, default=24)
    parser.add_argument("--axial-cells", type=int, default=48)
    parser.add_argument("--modes", type=int, default=4)
    args = parser.parse_args(argv)
    try:
        radius = _positive(args.outer_radius, "outer_radius")
        height = _positive(args.axis_half_height, "axis_half_height")
        nr = _count(args.radial_cells, "radial_cells", 2)
        radii = np.linspace(0, radius, nr + 1)
        heights = np.full_like(radii, height) if args.profile == "cylinder" else integral_profile(
            radii, height, args.coil_radius, historical=args.profile == "historical"
        )
        modes = solve_modes(meridian_mesh(radii, heights, args.axial_cells), args.modes)
    except ValueError as exc:
        parser.error(str(exc))
    report = {
        "model": "empty closed vacuum PEC cavity; axisymmetric H_phi sector only",
        "profile": args.profile,
        "outer_radius_m": radius, "axis_half_height_m": height,
        "radial_cells": nr, "axial_cells": args.axial_cells,
        "frequencies_hz": modes.frequencies_hz.tolist(),
        "algebraic_relative_residuals": modes.relative_residuals.tolist(),
        "electric_magnetic_energy_ratios": modes.electric_magnetic_energy_ratios.tolist(),
        "limitations": "No ports, source, toroidal core, losses, drive amplitude or physical scalar mode; refine the mesh.",
    }
    if args.profile == "cylinder":
        report["analytic_modes"] = cylindrical_tm_modes(radius, 2 * height, args.modes)
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
