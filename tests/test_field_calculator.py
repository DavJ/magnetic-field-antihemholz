"""Physical reference checks independent of the elliptic-integral formulas."""

from math import cos, pi

import numpy as np
import pytest
from scipy.integrate import quad

from field_calculator import AntiHemholtz, LoopField


def biot_savart(B0, radius, offset, r, z):
    """Integrate dl cross (point - wire) / distance**3 around one loop."""
    local_z = z - offset

    def integrand(phi, component):
        distance_squared = radius ** 2 + r ** 2 + local_z ** 2 - 2 * radius * r * cos(phi)
        cross = (radius * local_z * cos(phi), radius ** 2 - radius * r * cos(phi))
        return cross[component] / distance_squared ** 1.5

    return np.array([
        B0 * radius / (2 * pi) * quad(
            integrand, 0, 2 * pi, args=(component,), epsabs=1e-11, epsrel=1e-11
        )[0]
        for component in range(2)
    ])


@pytest.mark.parametrize("B0", [1.2, -0.8])
@pytest.mark.parametrize("r,z", [(0, -0.4), (0, 0.8), (0.2, 0.8), (1.2, 0.1), (2.4, -1.1)])
def test_loop_matches_direct_biot_savart(B0, r, z):
    loop = LoopField(B0, 1.7, -0.4)
    np.testing.assert_allclose(loop.vector(r, z), biot_savart(B0, 1.7, -0.4, r, z), rtol=2e-11, atol=2e-13)


def test_loop_offset_applies_once_to_all_field_components():
    centered = LoopField(2.1, 1.3)
    shifted = LoopField(2.1, 1.3, 0.7)
    assert shifted.Bz(0, 0.7) == 2.1
    for r, local_z in [(0, 0.4), (0.3, -0.9), (2.0, 0.2)]:
        np.testing.assert_allclose(shifted.vector(r, local_z + 0.7), centered.vector(r, local_z), rtol=1e-13, atol=1e-14)
        assert shifted.Baxial(local_z + 0.7) == pytest.approx(centered.Baxial(local_z))
        assert shifted.q(r, local_z + 0.7) == pytest.approx(centered.q(r, local_z))
        assert shifted.k(r, local_z + 0.7) == pytest.approx(centered.k(r, local_z))


def test_pair_matches_two_independent_biot_savart_loops():
    pair = AntiHemholtz(0.8, 1.4, 0.6)
    for r, z in [(0, 0.3), (0.5, 0.8), (1.8, -0.2)]:
        reference = biot_savart(0.8, 1.4, 0.6, r, z) + biot_savart(-0.8, 1.4, -0.6, r, z)
        np.testing.assert_allclose(pair.vector(r, z), reference, rtol=2e-11, atol=2e-13)


def test_pair_center_symmetry_and_axial_gradient():
    B0, radius, offset = 0.8, 1.4, 0.6
    pair = AntiHemholtz(B0, radius, offset)
    assert pair.vector(0, 0) == (0.0, 0.0)
    assert pair.Bz(0.4, -0.3) == pytest.approx(-pair.Bz(0.4, 0.3))
    assert pair.Br(0.4, -0.3) == pytest.approx(pair.Br(0.4, 0.3))
    h = 1e-5
    expected = 6 * B0 * radius ** 3 * offset / (radius ** 2 + offset ** 2) ** 2.5
    observed = (pair.Bz(0, h) - pair.Bz(0, -h)) / (2 * h)
    assert observed == pytest.approx(expected, rel=1e-9)


def test_normal_is_unit_and_perpendicular():
    pair = AntiHemholtz(1, 1, 0.5)
    field = pair.vector(1.05, 0.5)
    normal = pair.normal_vector(1.05, 0.5)
    assert np.linalg.norm(normal) == pytest.approx(1)
    assert np.dot(field, normal) == pytest.approx(0, abs=1e-14)


def test_field_null_and_ideal_wire_raise_clear_errors():
    pair = AntiHemholtz(1, 1, 0.5)
    with pytest.raises(ValueError, match="normal is undefined"):
        pair.normal_vector(0, 0)
    for z in [-0.5, 0.5]:
        with pytest.raises(ValueError, match="singular"):
            pair.vector(1, z)
    assert LoopField(0, 1, 0.5).vector(1, 0.5) == (0.0, 0.0)


def test_near_wire_field_remains_finite():
    loop = LoopField(1, 1)
    with np.errstate(all="raise"):
        field = loop.vector(1, 1e-10)
    assert np.all(np.isfinite(field))
    assert field[0] > 0


@pytest.mark.parametrize("r", [1e-5, 1e-8, 1e-10])
def test_near_axis_radial_field_has_correct_linear_limit(r):
    loop = LoopField(1, 1)
    z = 0.37
    # From div B = 0: Br/r -> -0.5 * dBaxis/dz as r -> 0.
    expected_slope = 1.5 * z / (1 + z ** 2) ** 2.5
    assert loop.Br(r, z) / r == pytest.approx(expected_slope, rel=1e-6)
    assert loop.Bz(r, z) == pytest.approx(loop.Baxial(z), rel=1e-9)


@pytest.mark.parametrize("r", [0.0008, 0.0012])
@pytest.mark.parametrize("z", [-1.7, 0, 0.37])
def test_near_axis_expansion_and_elliptic_formula_match_quadrature(r, z):
    np.testing.assert_allclose(
        LoopField(1, 1).vector(r, z), biot_savart(1, 1, 0, r, z),
        rtol=1e-9, atol=1e-13,
    )


@pytest.mark.parametrize("kwargs", [
    {"B0": float("nan")}, {"B0": float("inf")},
    {"loop_radius": 0}, {"loop_radius": -1}, {"loop_radius": float("inf")},
    {"axial_offset": float("nan")},
])
def test_invalid_loop_parameters_are_rejected(kwargs):
    parameters = dict(B0=1, loop_radius=1, axial_offset=0)
    parameters.update(kwargs)
    with pytest.raises(ValueError):
        LoopField(**parameters)


@pytest.mark.parametrize("r,z", [(-1, 0), (float("nan"), 0), (1, float("inf"))])
def test_invalid_coordinates_are_rejected(r, z):
    with pytest.raises(ValueError):
        LoopField(1, 1).vector(r, z)
