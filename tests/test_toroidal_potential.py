"""Check A normalization, linked flux, gauge freedom, and cylindrical calculus."""

import json

import numpy as np
import pytest

from toroidal_potential import (
    GaussianGauge, GaugeTransformedPotential, OpposedToroidsPotential,
    ThinToroidPotential, circulation, diagnose, main,
)


@pytest.mark.parametrize("flux", [1.7, -0.8, 0])
def test_single_toroid_axis_potential_has_flux_normalization(flux):
    radius, offset, z = 1.3, 0.4, -0.2
    potential = ThinToroidPotential(flux, radius, offset)
    expected = flux * radius ** 2 / (2 * (radius ** 2 + (z - offset) ** 2) ** 1.5)
    ar, az = potential.vector(0, z)
    assert ar == 0
    assert az == pytest.approx(expected)


@pytest.mark.parametrize("model", [ThinToroidPotential, OpposedToroidsPotential])
@pytest.mark.parametrize("flux", [1.2, -0.8])
@pytest.mark.parametrize("path_radius", [0.1, 0.25])
def test_linking_circulation_equals_core_flux(model, flux, path_radius):
    potential = model(flux, 1, 0.5)
    assert circulation(potential, 1, 0.5, path_radius) == pytest.approx(flux, rel=2e-10, abs=1e-12)


def test_unlinked_circulation_is_zero_and_opposite_core_has_negative_flux():
    pair = OpposedToroidsPotential(1.2, 1, 0.5)
    assert circulation(pair, 1, 0, 0.1) == pytest.approx(0, abs=1e-12)
    assert circulation(pair, 1, -0.5, 0.2) == pytest.approx(-1.2, rel=2e-10)


@pytest.mark.parametrize("model", [ThinToroidPotential, OpposedToroidsPotential])
@pytest.mark.parametrize("r,z", [(0, 0), (0, 0.2), (0.2, 0.1), (0.6, 0.25), (1.4, -0.2)])
def test_coulomb_gauge_is_divergence_free_and_external_B_vanishes(model, r, z):
    diagnostics = diagnose(model(1, 1, 0.5), r, z)
    assert diagnostics.divergence == pytest.approx(0, abs=2e-9)
    assert diagnostics.curl_phi == pytest.approx(0, abs=2e-9)


def test_nonzero_center_gradient_has_zero_three_dimensional_trace():
    flux, radius, offset = 1.7, 1.3, 0.45
    potential = OpposedToroidsPotential(flux, radius, offset)
    assert potential.vector(0, 0) == (0.0, 0.0)
    # These linear-limit quotients have O(h^2) truncation error.
    h = 1e-5 * radius
    gradient = 3 * flux * radius ** 2 * offset / (radius ** 2 + offset ** 2) ** 2.5
    axial = (potential.vector(0, h)[1] - potential.vector(0, -h)[1]) / (2 * h)
    radial = potential.vector(h, 0)[0] / h
    assert gradient > 0
    assert axial == pytest.approx(gradient, rel=1e-7)
    assert radial == pytest.approx(-gradient / 2, rel=1e-7)
    assert 2 * radial + axial == pytest.approx(0, abs=1e-8)


@pytest.mark.parametrize("r,z", [(0, 0), (0, 0.2), (0.2, 0.1), (0.6, -0.3)])
def test_static_gauge_changes_divergence_by_laplacian_but_preserves_curl(r, z):
    potential = OpposedToroidsPotential(0.8, 1, 0.5)
    gauge = GaussianGauge(0.4, 0.7)
    shifted = GaugeTransformedPotential(potential, gauge)
    before = diagnose(potential, r, z)
    after = diagnose(shifted, r, z)
    change = gauge.laplacian(r, z)
    assert abs(change) > 0.1
    assert after.divergence - before.divergence == pytest.approx(change, rel=1e-8, abs=1e-9)
    assert after.curl_phi == pytest.approx(before.curl_phi, abs=2e-9)


@pytest.mark.parametrize("center_z,expected", [(0.5, 0.7), (0, 0)])
def test_static_gauge_preserves_closed_loop_circulation(center_z, expected):
    potential = OpposedToroidsPotential(0.7, 1, 0.5)
    shifted = GaugeTransformedPotential(potential, GaussianGauge(2.4, 0.8))
    assert circulation(potential, 1, center_z, 0.15) == pytest.approx(expected, abs=1e-11)
    assert circulation(shifted, 1, center_z, 0.15) == pytest.approx(expected, abs=1e-11)


@pytest.mark.parametrize("r", [0, 1e-5, 0.5])
def test_diagnostics_include_radial_geometry_instead_of_planar_divergence(r):
    class RadialGradient:
        def vector(self, radius, z):
            return radius, z

    # A = (x, y, z): divergence is 3, including at the cylindrical axis.
    diagnostics = diagnose(RadialGradient(), r, 0.2)
    assert diagnostics.divergence == pytest.approx(3, abs=1e-10)
    assert diagnostics.curl_phi == pytest.approx(0, abs=1e-10)


def test_cli_reports_si_flux_and_a_gauge_change_without_a_new_B(capsys):
    assert main([]) == 0
    report = json.loads(capsys.readouterr().out)
    expected = report["flux_Wb"]
    assert report["linking_circulation_Wb"]["coulomb"] == pytest.approx(expected, rel=1e-9)
    assert report["linking_circulation_Wb"]["gauge_transformed"] == pytest.approx(expected, rel=1e-9)
    for sample in report["samples"]:
        assert np.isfinite(sample["A_coulomb_T_m"]).all()
        assert abs(sample["coulomb"]["divergence"]) < 1e-11
        assert abs(sample["gauge_transformed"]["curl_phi"]) < 1e-11
        assert sample["gauge_transformed"]["divergence"] == pytest.approx(
            sample["expected_divergence_change_T"], rel=1e-8, abs=1e-11,
        )


@pytest.mark.parametrize("model,args", [
    (ThinToroidPotential, (1, 0)), (ThinToroidPotential, (float("nan"), 1)),
    (OpposedToroidsPotential, (1, 1, 0)), (OpposedToroidsPotential, (1, 1, -1)),
    (GaussianGauge, (1, 0)), (GaussianGauge, (float("inf"), 1)),
])
def test_invalid_physical_or_gauge_parameters_are_rejected(model, args):
    with pytest.raises(ValueError):
        model(*args)


def test_invalid_diagnostic_step_and_circulation_path_are_rejected():
    potential = ThinToroidPotential(1, 1)
    with pytest.raises(ValueError, match="step"):
        diagnose(potential, 0.2, 0, step=0)
    with pytest.raises(ValueError, match="path_radius"):
        circulation(potential, 1, 0, 1)


def test_cli_reports_invalid_geometry():
    with pytest.raises(SystemExit) as exc:
        main(["--radius", "0"])
    assert exc.value.code == 2
