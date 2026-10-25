"""Retardation, charge completeness and full detector gauge coupling."""

import numpy as np
import pytest

from scalar_gate import (
    LocalizedGaugePulse, analytic_open_divergence, diagnose_potentials,
    probe_gauge_integral, retarded_path_potential,
)


def test_diagnostic_signs_against_uniform_physical_fields():
    def potential(t, xyz):
        x, y, z = xyz
        return -x, np.array([0, 2 * x, 0])
    result = diagnose_potentials(potential, 0.3, [0.4, -0.1, 0.7])
    np.testing.assert_allclose(result["E"], [1, 0, 0], atol=1e-10)
    np.testing.assert_allclose(result["B"], [0, 0, 2], atol=1e-10)
    assert abs(result["G"]) < 1e-10


@pytest.mark.parametrize("time,position", [
    (0.2, [0, 0, 0]), (0.8, [0.3, -0.2, 0.4]), (1.3, [-0.7, 0.2, 0.8]),
])
def test_time_dependent_pure_gauge_has_zero_fields_but_nonzero_four_divergence(time, position):
    gauge = LocalizedGaugePulse()
    result = diagnose_potentials(gauge.potential, time, position)
    np.testing.assert_allclose(result["E"], 0, atol=1e-9)
    np.testing.assert_allclose(result["B"], 0, atol=1e-9)
    assert abs(result["G"]) > 0.01
    assert result["G"] == pytest.approx(gauge.analytic_g(time, position), abs=1e-9)


@pytest.mark.parametrize("time,position", [
    (1.1, [0.3, -0.2, 0.4]), (0.4, [-0.2, 0.35, 0.6]),
])
def test_conserved_closed_loop_retarded_potential_has_zero_G_and_nonzero_fields(time, position):
    result = diagnose_potentials(retarded_path_potential, time, position)
    assert abs(result["G"]) < 1e-9
    assert np.linalg.norm(result["E"]) > 0.01
    assert np.linalg.norm(result["B"]) > 0.01


def test_endpoint_charge_cancels_the_exact_retarded_boundary_term():
    t, xyz = 1.1, [0.3, -0.2, 0.4]
    complete = lambda time, position: retarded_path_potential(time, position, closed=False)
    incomplete = lambda time, position: retarded_path_potential(time, position, closed=False, include_endpoint_charge=False)
    good, bad = diagnose_potentials(complete, t, xyz), diagnose_potentials(incomplete, t, xyz)
    boundary = analytic_open_divergence(t, xyz)
    assert abs(boundary) > 0.03
    assert bad["G"] == pytest.approx(boundary, abs=1e-9)
    assert good["div_A"] == pytest.approx(boundary, abs=1e-9)
    assert good["phi_t"] == pytest.approx(-boundary, abs=1e-9)
    assert abs(good["G"]) < 1e-9
    # Adding the physical endpoint charges is not a gauge transformation:
    # it corrects E as well as the source, while leaving this J-generated B.
    np.testing.assert_allclose(good["B"], bad["B"], atol=1e-10)
    assert np.linalg.norm(np.array(good["E"]) - np.array(bad["E"])) > 0.1


def test_retarded_boundary_control_converges_under_independent_path_quadrature():
    vectors = [retarded_path_potential(1.1, [0.3, -0.2, 0.4], quadrature_order=n)[1] for n in (64, 96, 128)]
    np.testing.assert_allclose(vectors[0], vectors[-1], atol=1e-12)
    np.testing.assert_allclose(vectors[1], vectors[-1], atol=1e-12)


def test_full_probe_interaction_cancels_against_independent_gaussian_moments():
    result = probe_gauge_integral()
    assert result["rho_chi_t_integral"] == pytest.approx(result["whole_space_exact_scalar_part"], abs=1e-12)
    assert result["J_chi_x_integral"] == pytest.approx(result["whole_space_exact_vector_part"], abs=1e-12)
    assert abs(result["complete_interaction_change"]) < 1e-12
    assert abs(result["omitting_scalar_potential_gives_false_change"]) > 0.3


def test_invalid_calculus_inputs_are_rejected():
    with pytest.raises(ValueError):
        diagnose_potentials(LocalizedGaugePulse().potential, 0, [0, 0], step=1e-3)
    with pytest.raises(ValueError):
        diagnose_potentials(LocalizedGaugePulse().potential, 0, [0, 0, 0], step=0)
    with pytest.raises(ValueError):
        retarded_path_potential(0, [0, 0, 0], quadrature_order=2)
