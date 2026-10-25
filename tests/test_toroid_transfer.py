"""Energy, circuit limits, independent closed form and magnetic-null controls."""

from dataclasses import replace
import json
from math import pi

import numpy as np
import pytest
from scipy.constants import mu_0
from scipy.integrate import quad

from toroid_transfer import CoupledCircuit, RectangularToroid, main, reference_circuit


def test_finite_cross_section_inductance_matches_integrated_magnetic_energy():
    geometry = RectangularToroid(0.03, 0.11, 0.017, 73)
    current = 0.31
    energy = quad(
        lambda r: (mu_0 * geometry.drive_turns * current / (2 * pi * r)) ** 2
        / (2 * mu_0) * 2 * pi * r * geometry.axial_height_m,
        geometry.inner_radius_m, geometry.outer_radius_m,
    )[0]
    assert 2 * energy / current ** 2 == pytest.approx(geometry.primary_inductance_h, rel=1e-11)
    assert geometry.linked_mutual_h(-13) == pytest.approx(-13 * geometry.flux_per_ampere_wb_a)


@pytest.mark.parametrize("load", [50, 1e6])
def test_full_loaded_transfer_matches_independent_reflected_impedance_formula(load):
    circuit = replace(reference_circuit(load_resistance_ohm=load), feedthrough_capacitance_f=0)
    f = np.geomspace(1, 1e6, 61)
    jw = 2j * pi * f
    z1 = circuit.primary_resistance_ohm + jw * circuit.primary_inductance_h
    z2 = circuit.pickup_resistance_ohm + jw * circuit.pickup_inductance_h
    zl = 1 / (1 / load + jw * circuit.load_capacitance_f)
    effective = z1 - (jw * circuit.mutual_inductance_h) ** 2 / (z2 + zl)
    expected_i1 = 1 / (circuit.source_resistance_ohm + (1 + circuit.source_resistance_ohm * jw * circuit.primary_capacitance_f) * effective)
    expected_v2 = jw * circuit.mutual_inductance_h * expected_i1 * zl / (z2 + zl)
    result = circuit.response(f)
    np.testing.assert_allclose(result.primary_current_a, expected_i1, rtol=2e-12)
    np.testing.assert_allclose(result.output_voltage_v, expected_v2, rtol=2e-12)


def test_open_pickup_has_faraday_voltage_and_mutual_reciprocity():
    circuit = replace(reference_circuit(), load_resistance_ohm=None, load_capacitance_f=0, feedthrough_capacitance_f=0)
    f = np.geomspace(10, 1e6, 31)
    for system in (circuit, replace(
        circuit, primary_inductance_h=circuit.pickup_inductance_h,
        pickup_inductance_h=circuit.primary_inductance_h,
        primary_resistance_ohm=circuit.pickup_resistance_ohm,
        pickup_resistance_ohm=circuit.primary_resistance_ohm,
    )):
        result = system.response(f)
        np.testing.assert_allclose(result.pickup_current_a, 0, atol=1e-14)
        np.testing.assert_allclose(result.output_voltage_v / result.primary_current_a, 2j * pi * f * system.mutual_inductance_h, rtol=1e-12)


def test_dc_has_no_induced_or_capacitive_output_and_uses_dc_power():
    circuit = reference_circuit()
    result = circuit.response(0, 2)
    current = 2 / (circuit.source_resistance_ohm + circuit.primary_resistance_ohm)
    assert result.primary_current_a[0] == pytest.approx(current)
    assert result.output_voltage_v[0] == 0
    assert result.source_power_w[0] == pytest.approx(2 * current)
    assert result.resistor_power_w[0] == pytest.approx(2 * current)
    with pytest.raises(ValueError, match="DC"):
        circuit.response(0, 1j)


@pytest.mark.parametrize("polarity", [-1, 0, 1])
def test_passive_coupled_network_conserves_average_power_with_capacitive_feedthrough(polarity):
    l1, l2 = 0.003, 0.0007
    circuit = CoupledCircuit(
        l1, l2, polarity * 0.8 * np.sqrt(l1 * l2), primary_resistance_ohm=2,
        pickup_resistance_ohm=4, load_resistance_ohm=1000,
        primary_capacitance_f=200e-12, load_capacitance_f=1e-9, feedthrough_capacitance_f=40e-12,
    )
    result = circuit.response(np.geomspace(1, 1e6, 83), 1.3 * np.exp(0.7j))
    assert np.max(result.relative_power_errors) < 2e-12
    assert np.all(result.source_power_w >= 0)
    assert np.all(result.resistor_power_w >= 0)
    np.testing.assert_allclose(result.source_power_w, result.resistor_power_w, rtol=2e-12)


def test_reversing_magnetic_polarity_flips_output_when_capacitive_feedthrough_is_absent():
    circuit = reference_circuit(feedthrough_capacitance_f=0)
    reversed_circuit = replace(circuit, mutual_inductance_h=-circuit.mutual_inductance_h)
    f = np.geomspace(100, 1e6, 51)
    forward, backward = circuit.response(f), reversed_circuit.response(f)
    np.testing.assert_allclose(forward.output_voltage_v, -backward.output_voltage_v, rtol=1e-12)
    np.testing.assert_allclose(forward.primary_current_a, backward.primary_current_a, rtol=1e-12)
    np.testing.assert_allclose(forward.source_power_w, backward.source_power_w, rtol=1e-12)


def test_capacitive_feedthrough_survives_a_magnetic_null():
    control = reference_circuit("control")
    ideal_null = replace(control, feedthrough_capacitance_f=0).response([10000, 1e6])
    capacitive_null = control.response([10000, 1e6])
    np.testing.assert_array_equal(ideal_null.output_voltage_v, 0)
    assert np.all(abs(capacitive_null.output_voltage_v) > 1e-10)
    assert np.max(capacitive_null.relative_power_errors) < 1e-12


def test_series_opposed_toroids_cancel_and_known_mismatch_restores_ordinary_induction():
    single = reference_circuit()
    balanced = reference_circuit("balanced-pair", feedthrough_capacitance_f=0)
    mismatch = reference_circuit("mismatched-pair", feedthrough_capacitance_f=0)
    summed = reference_circuit("sum-pair", feedthrough_capacitance_f=0)
    assert balanced.mutual_inductance_h == 0
    assert balanced.primary_inductance_h == pytest.approx(2 * single.primary_inductance_h)
    assert mismatch.mutual_inductance_h == pytest.approx(-0.01 * single.mutual_inductance_h)
    assert summed.mutual_inductance_h == pytest.approx(2 * single.mutual_inductance_h)
    assert balanced.response(10000).output_voltage_v[0] == 0
    assert abs(mismatch.response(10000).output_voltage_v[0]) > 4e-6


def test_zero_source_and_linear_amplitude_scaling():
    circuit = reference_circuit()
    np.testing.assert_array_equal(circuit.response([1000, 10000], 0).output_voltage_v, 0)
    unit = circuit.response([1000, 10000])
    scaled = circuit.response([1000, 10000], 3)
    np.testing.assert_allclose(scaled.output_voltage_v, 3 * unit.output_voltage_v, rtol=1e-12)
    np.testing.assert_allclose(scaled.source_power_w, 9 * unit.source_power_w, rtol=1e-12)


@pytest.mark.parametrize("change", [
    {"mutual_inductance_h": 1}, {"primary_resistance_ohm": -1},
    {"feedthrough_capacitance_f": -1}, {"source_resistance_ohm": 0},
    {"load_resistance_ohm": float("inf")},
])
def test_unphysical_or_nonpassive_component_parameters_are_rejected(change):
    with pytest.raises(ValueError):
        replace(reference_circuit(), **change)


def test_invalid_geometry_and_frequency_inputs_are_rejected():
    with pytest.raises(ValueError):
        RectangularToroid(0.1, 0.09)
    with pytest.raises(ValueError):
        RectangularToroid(drive_turns=True)
    for frequency in ([], [-1], [float("nan")], [[1, 2]]):
        with pytest.raises(ValueError):
            reference_circuit().response(frequency)


def test_transfer_cli_reports_loaded_output_and_power_balance(capsys):
    assert main([]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["primary_current_a"]["peak_magnitude"] == pytest.approx(0.0197822, rel=1e-5)
    assert report["output_voltage_v"]["peak_magnitude"] == pytest.approx(0.000498849, rel=1e-5)
    assert report["relative_power_error"] < 1e-12
    with pytest.raises(SystemExit) as exc:
        main(["--feed-capacitance-pf", "-1"])
    assert exc.value.code == 2
