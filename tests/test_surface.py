"""Geometry and actual STL export regressions for the original generator."""

from collections import Counter
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from stl import mesh

from main import generate_surface, main, save_stl, surface_of_revolution


@pytest.mark.parametrize("vertical,horizontal", [(20, 120), (5, 9), (1, 3)])
def test_mesh_is_axisymmetric_nondegenerate_and_has_only_two_boundary_rings(vertical, horizontal):
    vertices, faces = generate_surface(steps_vertical=vertical, steps_horizontal=horizontal)
    assert vertices.shape == (horizontal * (vertical + 1), 3)
    assert faces.shape == (2 * horizontal * vertical, 3)
    assert np.all(np.isfinite(vertices))
    assert faces.min() == 0
    assert faces.max() == len(vertices) - 1
    assert len(np.unique(faces)) == len(vertices)

    slices = vertices.reshape(horizontal, vertical + 1, 3)
    for angular_slice in slices:
        np.testing.assert_allclose(np.hypot(angular_slice[:, 0], angular_slice[:, 1]), slices[0, :, 0], atol=1e-14)
        np.testing.assert_allclose(angular_slice[:, 2], slices[0, :, 2], atol=1e-14)
    assert len(np.unique(vertices, axis=0)) == len(vertices)
    triangles = vertices[faces]
    areas = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]), axis=1)
    assert np.all(areas > 1e-12)
    assert len(np.unique(np.sort(faces, axis=1), axis=0)) == len(faces)

    directed = Counter()
    undirected = Counter()
    for a, b, c in faces:
        for start, end in [(a, b), (b, c), (c, a)]:
            directed[start, end] += 1
            undirected[tuple(sorted((start, end)))] += 1
    assert set(undirected.values()) <= {1, 2}
    boundary = {edge for edge, count in undirected.items() if count == 1}
    assert len(boundary) == 2 * horizontal
    for start, end in boundary:
        assert start % (vertical + 1) == end % (vertical + 1)
        assert start % (vertical + 1) in (0, vertical)
    for (start, end), count in undirected.items():
        if count == 2:
            assert directed[start, end] == directed[end, start] == 1


def test_changing_only_angular_resolution_preserves_the_profile():
    coarse, _ = generate_surface(steps_horizontal=6)
    fine, _ = generate_surface(steps_horizontal=30)
    np.testing.assert_array_equal(coarse[:21], fine[:21])


def test_scaling_all_lengths_scales_the_surface():
    original, original_faces = generate_surface()
    scaled, scaled_faces = generate_surface(
        generator_radius=7, offset_radius=0.35, generator_offset=3.5, offset_offset=0.175
    )
    np.testing.assert_allclose(scaled, 7 * original, rtol=1e-11, atol=1e-12)
    np.testing.assert_array_equal(scaled_faces, original_faces)


def test_stl_round_trip_preserves_finite_nondegenerate_triangles(tmp_path):
    vertices, faces = generate_surface(steps_vertical=12, steps_horizontal=18)
    output = tmp_path / "surface.stl"
    save_stl(vertices, faces, output)
    loaded = mesh.Mesh.from_file(str(output))
    np.testing.assert_allclose(loaded.vectors, vertices[faces], rtol=1e-6, atol=1e-7)
    assert np.all(np.isfinite(loaded.vectors))
    assert np.all(np.linalg.norm(loaded.normals, axis=1) > 0)


def test_import_main_has_no_output_side_effects(tmp_path):
    source = str(Path(__file__).resolve().parents[1])
    result = subprocess.run(
        [sys.executable, "-c", "import main"], cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": source}, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert list(tmp_path.iterdir()) == []


def test_script_without_arguments_generates_valid_default_stl(tmp_path):
    script = Path(__file__).resolve().parents[1] / "main.py"
    result = subprocess.run([sys.executable, str(script)], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    loaded = mesh.Mesh.from_file(str(tmp_path / "ufo.stl"))
    assert len(loaded.vectors) == 4800
    assert np.all(np.isfinite(loaded.vectors))
    assert np.all(np.linalg.norm(loaded.normals, axis=1) > 0)


def test_cli_accepts_output_and_resolution(tmp_path):
    output = tmp_path / "custom.stl"
    assert main([str(output), "--steps-vertical", "4", "--steps-horizontal", "8"]) == 0
    assert len(mesh.Mesh.from_file(str(output)).vectors) == 64


@pytest.mark.parametrize("kwargs", [
    {"steps_vertical": 0}, {"steps_vertical": 1.5}, {"steps_vertical": True},
    {"steps_horizontal": 2}, {"generator_radius": 0}, {"offset_radius": 0},
    {"offset_offset": 0.5}, {"offset_offset": float("nan")},
])
def test_invalid_generator_parameters_are_rejected(kwargs):
    with pytest.raises(ValueError):
        generate_surface(**kwargs)


@pytest.mark.parametrize("profile", [[], [[1, 0]], [[1, 0], [1, 0]], [[0, 0], [1, 0]], [[1, 0], [1, float("nan")]]])
def test_invalid_profiles_are_rejected(profile):
    with pytest.raises(ValueError):
        surface_of_revolution(profile, 8)


def test_cli_reports_invalid_resolution_without_writing(tmp_path):
    output = tmp_path / "invalid.stl"
    with pytest.raises(SystemExit) as exc:
        main([str(output), "--steps-horizontal", "2"])
    assert exc.value.code == 2
    assert not output.exists()
