"""Trace a meridional profile and revolve it into an open STL surface."""

import argparse
from math import atan2, isfinite, pi
from pathlib import Path

import numpy as np

from field_calculator import AntiHemholtz


def _positive(value, name):
    value = float(value)
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _count(value, name, minimum):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def trace_profile(generator, start_r, start_z, steps, step_length):
    """Take fixed-length Euler steps perpendicular to the magnetic field."""
    steps = _count(steps, "steps", 1)
    step_length = _positive(step_length, "step_length")
    r = _positive(start_r, "start_r")
    z = float(start_z)
    if not isfinite(z):
        raise ValueError("start_z must be finite")
    profile = np.empty((steps + 1, 2), dtype=float)
    profile[0] = r, z
    for index in range(steps):
        nr, nz = generator.normal_vector(r, z)
        r, z = r + nr * step_length, z + nz * step_length
        if not isfinite(r) or not isfinite(z) or r <= 0:
            raise ValueError("profile left the finite, positive-radius domain")
        profile[index + 1] = r, z
    return profile


def surface_of_revolution(profile, steps_horizontal):
    """Return vertices and consistently wound triangles, with a wrapped seam.

    The two end rings remain open. Axis points are rejected because they
    require a separate triangle-fan topology.
    """
    steps_horizontal = _count(steps_horizontal, "steps_horizontal", 3)
    profile = np.asarray(profile, dtype=float)
    if profile.ndim != 2 or profile.shape[1] != 2 or len(profile) < 2:
        raise ValueError("profile must contain at least two (r, z) points")
    if not np.all(np.isfinite(profile)) or np.any(profile[:, 0] <= 0):
        raise ValueError("profile coordinates must be finite with positive radii")
    if np.any(np.all(np.diff(profile, axis=0) == 0, axis=1)):
        raise ValueError("consecutive profile points must be distinct")

    # Generate one profile, then reuse it at every azimuth. Do not duplicate
    # the 2*pi slice: modulo indexing closes the angular seam exactly.
    angles = np.arange(steps_horizontal) * (2 * pi / steps_horizontal)
    radii, heights = profile.T
    vertices = np.stack(
        (
            np.cos(angles)[:, None] * radii,
            np.sin(angles)[:, None] * radii,
            np.broadcast_to(heights, (steps_horizontal, len(profile))),
        ),
        axis=-1,
    ).reshape(-1, 3)

    slices = np.arange(steps_horizontal)[:, None]
    rows = np.arange(len(profile) - 1)
    current = slices * len(profile) + rows
    following = ((slices + 1) % steps_horizontal) * len(profile) + rows
    first = np.stack((current, current + 1, following), axis=-1)
    second = np.stack((current + 1, following + 1, following), axis=-1)
    faces = np.stack((first, second), axis=2).reshape(-1, 3)
    return vertices, faces


def generate_surface(
    generator_radius=1.0,
    offset_radius=0.05,
    generator_offset=0.5,
    offset_offset=0.025,
    steps_vertical=20,
    steps_horizontal=120,
):
    """Generate the original example's surface with corrected mesh topology.

    steps_vertical counts profile segments; each slice has one more vertex.
    The trace length scales with the coil radius, preserving the unit-radius
    example while allowing geometrically similar coils in other length units.
    """
    generator_radius = _positive(generator_radius, "generator_radius")
    offset_radius = _positive(offset_radius, "offset_radius")
    generator_offset = _positive(generator_offset, "generator_offset")
    offset_offset = float(offset_offset)
    if not isfinite(offset_offset) or not 0 <= offset_offset < generator_offset:
        raise ValueError("offset_offset must be finite and satisfy 0 <= offset_offset < generator_offset")
    steps_vertical = _count(steps_vertical, "steps_vertical", 1)
    steps_horizontal = _count(steps_horizontal, "steps_horizontal", 3)
    generator = AntiHemholtz(1, generator_radius, generator_offset)
    start_r = generator_radius + offset_radius
    theta0 = atan2(start_r, generator_offset - offset_offset)
    step_length = 2 * pi * generator_radius * (pi / 2 - theta0) / steps_vertical
    profile = trace_profile(generator, start_r, generator_offset, steps_vertical, step_length)
    return surface_of_revolution(profile, steps_horizontal)


def save_stl(vertices, faces, output):
    """Write indexed triangles as STL; importing this module writes no files."""
    from stl import mesh

    surface = mesh.Mesh(np.zeros(len(faces), dtype=mesh.Mesh.dtype))
    surface.vectors[:] = vertices[faces]
    surface.save(str(output))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path, default=Path("ufo.stl"))
    parser.add_argument("--steps-vertical", type=int, default=20, help="profile segment count (>= 1)")
    parser.add_argument("--steps-horizontal", type=int, default=120, help="azimuthal slice count (>= 3)")
    args = parser.parse_args(argv)
    try:
        vertices, faces = generate_surface(
            steps_vertical=args.steps_vertical, steps_horizontal=args.steps_horizontal
        )
    except ValueError as exc:
        parser.error(str(exc))
    save_stl(vertices, faces, args.output)
    print(f"Saved {args.output}: {len(vertices)} vertices, {len(faces)} triangles")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
