"""Independent Bessel-spectrum, field-shape, scaling and gauge controls."""

import json
from math import exp, pi

import numpy as np
import pytest
from scipy.constants import c
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import j1, jn_zeros, spherical_jn

from cavity_modes import assemble, cylindrical_tm_modes, main, meridian_mesh, solve_modes
from cavity_profile import historical_az, integral_profile, open_integral
from toroidal_potential import ThinToroidPotential


def cylinder(nr, nz=None, scale=1):
    r = np.linspace(0, 0.1 * scale, nr + 1)
    return meridian_mesh(r, np.full_like(r, 0.1 * scale), nz or 2 * nr)


def test_first_six_cylinder_frequencies_converge_to_independent_bessel_spectrum():
    # Includes axial p=0,1,2,3 and the second radial family: detects incorrect
    # PEC Robin conditions, scalar-Laplacian substitution and missing parity.
    exact = np.array([row["frequency_hz"] for row in cylindrical_tm_modes(0.1, 0.2, 6)])
    errors = []
    for nr in (8, 16, 32):
        result = solve_modes(cylinder(nr), 6)
        error = result.frequencies_hz / exact - 1
        assert np.all(error > -1e-10)
        errors.append(error)
    assert np.max(errors[-1]) < 0.00085
    # Doubling each mesh dimension reduces these errors by about four.
    assert np.all(errors[1] < 0.28 * errors[0])
    assert np.all(errors[2] < 0.28 * errors[1])
    assert exact[0] == pytest.approx(c * jn_zeros(0, 1)[0] / (2 * pi * 0.1))


def test_cylinder_fundamental_has_the_maxwell_magnetic_field_shape():
    result = solve_modes(cylinder(32), 1)
    r = result.mesh.nodes[:, 0]
    expected = j1(jn_zeros(0, 1)[0] * r / 0.1)
    expected /= expected.max()
    np.testing.assert_allclose(result.relative_h_phi(), expected, atol=3e-4)
    assert np.all(result.relative_h_phi()[r == 0] == 0)


def test_curved_wall_sequence_converges_to_a_spherical_maxwell_mode():
    # The sphere's electric l=1 mode obeys d_x[x*j1(x)]=0, not J0(x)=0.
    # This independently exercises the sloping-wall condition and Jacobian.
    root = brentq(lambda x: spherical_jn(1, x) + x * spherical_jn(1, x, derivative=True), 2, 3)
    exact = c * root / (2 * pi * 0.1)
    errors = []
    for nr in (16, 32, 64):
        # A small equatorial truncation avoids a zero-height, degenerate cell.
        # Its distance also tends to zero: this is a domain + mesh limit.
        r = np.linspace(0, 0.1 * (1 - 1 / (4 * nr ** 2)), nr + 1)
        z = np.sqrt(0.1 ** 2 - r ** 2)
        result = solve_modes(meridian_mesh(r, z, 2 * nr), 1)
        errors.append(abs(result.frequencies_hz[0] / exact - 1))
    assert errors[-1] < 0.001
    assert errors[1] < 0.45 * errors[0]
    assert errors[2] < 0.45 * errors[1]


def test_all_dimensions_scale_frequencies_inversely():
    baseline = solve_modes(cylinder(12), 4).frequencies_hz
    for scale in (0.01, 7):
        scaled = solve_modes(cylinder(12, scale=scale), 4).frequencies_hz
        np.testing.assert_allclose(scaled * scale, baseline, rtol=1e-10)


def test_lossless_modes_have_small_algebraic_residual_and_balanced_energy():
    r = np.linspace(0, 0.09, 25)
    mesh = meridian_mesh(r, integral_profile(r, 0.1, 0.1), 48)
    result = solve_modes(mesh, 4)
    k, m, _ = assemble(mesh)
    assert np.max(result.relative_residuals) < 1e-10
    np.testing.assert_allclose(result.electric_magnetic_energy_ratios, 1, atol=1e-10)
    np.testing.assert_allclose(result.w_vectors.T @ m @ result.w_vectors, np.eye(4), atol=1e-10)
    assert np.max(abs((k - k.T).data), initial=0) < 1e-13
    assert np.max(abs((m - m.T).data), initial=0) < 1e-13
    higher_quadrature = solve_modes(mesh, 4, quadrature_order=5)
    np.testing.assert_allclose(higher_quadrature.frequencies_hz, result.frequencies_hz, rtol=1e-9)


def test_profile_spectra_stabilize_with_refinement():
    values = []
    for nr in (12, 24, 48):
        r = np.linspace(0, 0.09, nr + 1)
        mesh = meridian_mesh(r, integral_profile(r, 0.1, 0.1), 2 * nr)
        values.append(solve_modes(mesh, 4).frequencies_hz)
    coarse_change = abs(values[1] / values[0] - 1)
    fine_change = abs(values[2] / values[1] - 1)
    assert np.all(fine_change < 0.3 * coarse_change)
    assert np.max(fine_change) < 0.0013


def test_reconstructed_profile_satisfies_a_direct_flux_ring_integral():
    # Integrate the inverse-curl ring kernel directly, independently of the
    # elliptic implementation used to obtain the profile.
    r = np.array([0, 0.3, 0.6, 0.9])
    height = integral_profile(r, 1, 1)
    for radius, top in zip(r, height):
        def az(z):
            return quad(
                lambda theta: (1 - radius * np.cos(theta))
                / (1 + radius ** 2 - 2 * radius * np.cos(theta) + z ** 2) ** 1.5,
                0, 2 * pi, epsabs=1e-11,
            )[0] / (4 * pi)
        assert open_integral(lambda unused, z: az(z), radius, top) == pytest.approx(1 / (2 * np.sqrt(2)), rel=1e-9)


def test_legacy_modulus_error_changes_off_axis_profile_and_respects_scaling():
    r = np.array([0, 0.25, 0.5, 0.75, 0.9])
    old = integral_profile(r, 1, 1, historical=True)
    fixed = integral_profile(r, 1, 1)
    assert old[0] == pytest.approx(fixed[0])
    assert np.all(fixed[1:] > old[1:] + 0.07)
    np.testing.assert_allclose(integral_profile(0.1 * r, 0.1, 0.1), 0.1 * fixed, rtol=1e-10)
    assert historical_az(0, 0) == pytest.approx(0.5)


def test_static_decaying_gauge_changes_open_integral_without_changing_axis_target():
    # chi=0.2*r^2*z*exp(-(r^2+z^2)/2) vanishes on the axis and at infinity.
    source = ThinToroidPotential(1, 1)
    az = lambda r, z: source.vector(r, z)[1]
    chi = lambda r, z: 0.2 * r ** 2 * z * exp(-(r ** 2 + z ** 2) / 2)
    shifted = lambda r, z: az(r, z) + 0.2 * r ** 2 * (1 - z ** 2) * exp(-(r ** 2 + z ** 2) / 2)
    r, height = 0.7, 0.6
    change = open_integral(shifted, r, height) - open_integral(az, r, height)
    assert change == pytest.approx(chi(r, height) - chi(r, 0), rel=1e-10)
    assert change > 0.03
    assert open_integral(shifted, 0, 1) == pytest.approx(open_integral(az, 0, 1))


def test_cavity_cli_produces_a_frequency_report_with_an_analytic_control(capsys):
    assert main(["--profile", "cylinder", "--radial-cells", "8", "--axial-cells", "16", "--modes", "3"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert len(report["frequencies_hz"]) == 3
    assert report["frequencies_hz"][0] == pytest.approx(report["analytic_modes"][0]["frequency_hz"], rel=1e-4)


def test_nonphysical_geometries_are_rejected():
    with pytest.raises(ValueError):
        meridian_mesh([0, 0.5, 1], [1, 1, 0], 8)
    with pytest.raises(ValueError):
        integral_profile([0, 0.5, 1], 1, 1)
    with pytest.raises(ValueError):
        solve_modes(cylinder(2, 2), 9)
    with pytest.raises(ValueError):
        meridian_mesh([0, 0.5, 1], [1, 1, 1], True)
    with pytest.raises(SystemExit) as exc:
        main(["--outer-radius", "0"])
    assert exc.value.code == 2
