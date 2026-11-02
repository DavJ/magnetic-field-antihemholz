"""Independent stress variation, curved-gauge and toroid energy controls."""

import numpy as np
import pytest
from scipy.constants import c, mu_0
from scipy.integrate import quad

from gravity_gate import (
    de_sitter_gauge_divergence, de_sitter_metric,
    diagnose_covariant_potential, energy_equivalent_scales,
    maxwell_stress_energy, polynomial_gauge_covector,
)
from toroid_transfer import RectangularToroid, reference_circuit


def test_stress_tensor_matches_independent_electric_magnetic_formulas():
    electric, magnetic = np.array([0.3, -0.5, 0.7]), np.array([0.2, 0.6, -0.4])
    field = np.zeros((4, 4))
    field[0, 1:], field[1:, 0] = -electric, electric
    field[1, 2], field[2, 3], field[3, 1] = magnetic[2], magnetic[0], magnetic[1]
    field[2, 1], field[3, 2], field[1, 3] = -magnetic[2], -magnetic[0], -magnetic[1]
    eta = np.diag([-1, 1, 1, 1])
    stress = maxwell_stress_energy(field, eta)
    u = (electric @ electric + magnetic @ magnetic) / 2
    assert stress[0, 0] == pytest.approx(u)
    np.testing.assert_allclose(stress[0, 1:], -np.cross(electric, magnetic), atol=1e-15)
    np.testing.assert_allclose(stress[1:, 1:],
                               u * np.eye(3) - np.outer(electric, electric) - np.outer(magnetic, magnetic))
    assert np.sum(eta * stress) == pytest.approx(0, abs=1e-15)
    # Static weak-field source is (T00+sum_i Tii)/c^2, not just u/c^2.
    assert stress[0, 0] + np.trace(stress[1:, 1:]) == pytest.approx(2 * u)


def test_stress_is_the_metric_variation_even_with_nondiagonal_metric():
    basis = np.array([[1, 0.2, 0, 0], [0, 1.1, 0.1, 0],
                      [0, 0, 0.9, 0.15], [0.1, 0, 0, 1.2]])
    metric = basis.T @ np.diag([-1, 1, 1, 1]) @ basis
    inverse = np.linalg.inv(metric)
    matrix = np.arange(16).reshape(4, 4) / 10
    field = matrix - matrix.T
    stress = maxwell_stress_energy(field, metric)

    def density(contravariant_metric):
        volume = np.sqrt(-1 / np.linalg.det(contravariant_metric))
        return -volume / 4 * np.sum(field * (contravariant_metric @ field @ contravariant_metric))

    for i in range(4):
        for j in range(i, 4):
            direction = np.zeros((4, 4))
            direction[i, j] = direction[j, i] = 1
            h = 2e-5
            derivative = (density(inverse + h * direction) - density(inverse - h * direction)) / (2 * h)
            expected = -0.5 * np.sqrt(-np.linalg.det(metric)) * np.sum(stress * direction)
            assert derivative == pytest.approx(expected, rel=2e-8, abs=1e-10)
    assert np.sum(inverse * stress) == pytest.approx(0, abs=1e-14)


@pytest.mark.parametrize("time,hubble", [(0.7, 0.4), (-0.3, 0.4), (0.2, 0.0)])
def test_pure_gauge_has_nonzero_G_but_no_maxwell_gravitational_source(time, hubble):
    q = [time, 0.3, -0.2, 0.4]
    metric = lambda point: de_sitter_metric(point, hubble)
    result = diagnose_covariant_potential(polynomial_gauge_covector, metric, q)
    assert result["G"] == pytest.approx(de_sitter_gauge_divergence(time, hubble), abs=1e-9)
    assert abs(result["G"]) > 0.1
    np.testing.assert_allclose(result["field_cov"], 0, atol=1e-12)
    np.testing.assert_allclose(result["stress_energy_cov"], 0, atol=1e-23)


def test_adding_pure_gauge_changes_G_but_preserves_nonzero_stress_on_curved_background():
    def magnetic_potential(q):
        return np.array([0, -0.6 * q[2], 0.6 * q[1], 0])
    transformed = lambda q: magnetic_potential(q) + polynomial_gauge_covector(q)
    q = [0.7, 0.3, -0.2, 0.4]
    original = diagnose_covariant_potential(magnetic_potential, de_sitter_metric, q)
    changed = diagnose_covariant_potential(transformed, de_sitter_metric, q)
    assert abs(changed["G"] - original["G"]) > 1
    np.testing.assert_allclose(changed["field_cov"], original["field_cov"], atol=1e-12)
    np.testing.assert_allclose(changed["stress_energy_cov"], original["stress_energy_cov"], atol=1e-12)
    # Constant comoving F_xy=1.2 gives a locally measured B=1.2/a^2.
    assert original["stress_energy_cov"][0][0] == pytest.approx(1.2 ** 2 / 2 * np.exp(-4 * 0.4 * q[0]), rel=1e-11)


def test_curved_G_check_converges_under_step_refinement():
    q, metric = [0.7, 0.3, -0.2, 0.4], de_sitter_metric
    exact = de_sitter_gauge_divergence(q[0])
    errors = [abs(diagnose_covariant_potential(polynomial_gauge_covector, metric, q, step=h)["G"] - exact)
              for h in (0.02, 0.01, 0.005)]
    assert errors[0] > 12 * errors[1] > 12 ** 2 * errors[2]


def test_toroid_field_energy_agrees_with_volume_integral_and_closed_form():
    toroid, current = RectangularToroid(), 1.0
    radial_integral = quad(lambda r: (mu_0 * toroid.drive_turns * current / (2 * np.pi * r)) ** 2
                           / (2 * mu_0) * 2 * np.pi * r * toroid.axial_height_m,
                           toroid.inner_radius_m, toroid.outer_radius_m, epsabs=1e-20)[0]
    assert radial_integral == pytest.approx(toroid.primary_inductance_h * current ** 2 / 2, rel=1e-13)
    assert energy_equivalent_scales(radial_integral)["energy_mass_equivalent_kg"] == pytest.approx(radial_integral / c ** 2)


def test_loaded_core_and_leakage_energy_match_full_inductance_quadratic_form():
    toroid, circuit = RectangularToroid(), reference_circuit("linked")
    response = circuit.response(10000)
    i1, i2 = response.primary_current_a[0], response.pickup_current_a[0]
    total = (circuit.primary_inductance_h * abs(i1) ** 2 + circuit.pickup_inductance_h * abs(i2) ** 2
             + 2 * circuit.mutual_inductance_h * (i1 * i2.conjugate()).real) / 4
    net_ampere_turns = toroid.drive_turns * i1 + 10 * i2
    core = toroid.permeance_h * abs(net_ampere_turns) ** 2 / 4
    leakage = (circuit.pickup_inductance_h - toroid.permeance_h * 10 ** 2) * abs(i2) ** 2 / 4
    assert core > 0 and leakage > 0
    assert total == pytest.approx(core + leakage, rel=1e-13)


def test_energy_scales_obey_energy_and_inverse_square_scaling():
    base, changed = energy_equivalent_scales(1, 1), energy_equivalent_scales(4, 2)
    assert base["energy_mass_equivalent_kg"] == pytest.approx(1.1126500560536185e-17, rel=1e-14)
    assert changed["energy_mass_equivalent_kg"] == 4 * base["energy_mass_equivalent_kg"]
    assert changed["point_mass_acceleration_scale_m_s2"] == base["point_mass_acceleration_scale_m_s2"]
    assert energy_equivalent_scales(0)["point_mass_acceleration_scale_m_s2"] == 0


@pytest.mark.parametrize("energy,distance", [(-1, 1), (np.nan, 1), (1, 0), (1, np.inf)])
def test_invalid_energy_scale_inputs(energy, distance):
    with pytest.raises(ValueError):
        energy_equivalent_scales(energy, distance)


def test_covariant_control_rejects_invalid_tensors_and_coordinates():
    with pytest.raises(ValueError):
        maxwell_stress_energy(np.eye(4), np.diag([-1, 1, 1, 1]))
    with pytest.raises(ValueError):
        maxwell_stress_energy(np.zeros((4, 4)), np.eye(4))
    with pytest.raises(ValueError):
        diagnose_covariant_potential(polynomial_gauge_covector, de_sitter_metric, [0, 0, 0])
