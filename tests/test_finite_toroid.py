"""Finite-core A checked with axis quadrature, holonomy and a thin-core limit."""
import numpy as np
import pytest
from scipy.integrate import dblquad
from scipy.constants import mu_0
from toroid_transfer import RectangularToroid
from toroidal_potential import ThinToroidPotential, diagnose, circulation, GaussianGauge, GaugeTransformedPotential
from finite_toroid import FiniteToroidPotential, OpposedFiniteToroidsPotential, magnetic_energy_harmonics


def test_flux_and_axis_potential_match_independent_continuous_integrals():
    geometry=RectangularToroid(.7,1.3,.2,10)
    field=FiniteToroidPotential(geometry,2,.4,order=12)
    assert field.quadrature_flux_wb == pytest.approx(field.flux_wb,rel=1e-12)
    z=.15
    expected=dblquad(lambda h,r: mu_0*10*2/(2*np.pi*r)*r*r/(2*(r*r+(z-h)**2)**1.5),
                     .7,1.3,lambda r:.3,lambda r:.5,epsabs=1e-16,epsrel=1e-11)[0]
    assert field.vector(0,z)[1] == pytest.approx(expected,rel=1e-11)


def test_finite_pair_is_divergence_free_with_nonzero_traceless_gradient():
    pair=OpposedFiniteToroidsPotential(RectangularToroid(.9,1.1,.1,20),3,.5,8)
    h=1e-4
    gradient=pair.center_gradient
    assert gradient>0
    axial=(pair.vector(0,h)[1]-pair.vector(0,-h)[1])/(2*h)
    radial=pair.vector(h,0)[0]/h
    assert axial==pytest.approx(gradient,rel=1e-7)
    assert radial==pytest.approx(-gradient/2,rel=1e-7)
    for r,z in [(0,0),(.2,.1),(.5,-.2)]:
        d=diagnose(pair,r,z,step=2e-4)
        assert abs(d.divergence/gradient)<2e-8
        assert abs(d.curl_phi/gradient)<2e-8


def test_finite_core_holonomy_and_smooth_gauge_control():
    pair=OpposedFiniteToroidsPotential(RectangularToroid(.9,1.1,.1,20),3,.5,6)
    flux=pair.upper.flux_wb
    shifted=GaugeTransformedPotential(pair,GaussianGauge(flux,.8))
    assert circulation(pair,1,.5,.2)==pytest.approx(flux,rel=1e-10)
    assert circulation(shifted,1,.5,.2)==pytest.approx(flux,rel=1e-10)
    assert circulation(pair,1,-.5,.2)==pytest.approx(-flux,rel=1e-10)
    assert circulation(pair,1,0,.1)==pytest.approx(0,abs=abs(flux)*1e-10)


def test_thin_core_limit_at_fixed_current_and_normalized_flux():
    errors=[]
    for width in (.2,.1,.05):
        field=FiniteToroidPotential(RectangularToroid(1-width/2,1+width/2,width,20),1,.5,8)
        thin=ThinToroidPotential(field.flux_wb,1,.5)
        errors.append(np.linalg.norm(np.array(field.vector(.3,.1))-thin.vector(.3,.1))/abs(field.flux_wb))
    assert errors[1]<errors[0]/3
    assert errors[2]<errors[1]/3


def test_energy_harmonics_equal_independent_period_integrals_and_even_reversal():
    L,dc,ac=2e-3,.2,.5
    h=magnetic_energy_harmonics(L,dc,ac)
    phases=np.arange(4096)*2*np.pi/4096
    energy=L*(dc+ac*np.cos(phases))**2/2
    assert h['dc_j']==pytest.approx(np.mean(energy))
    assert h['omega_cos_j']==pytest.approx(2*np.mean(energy*np.cos(phases)))
    assert h['two_omega_cos_j']==pytest.approx(2*np.mean(energy*np.cos(2*phases)))
    assert magnetic_energy_harmonics(L,-dc,-ac)==h
    assert magnetic_energy_harmonics(L,0,ac)['omega_cos_j']==0


@pytest.mark.parametrize('kwargs',[{'order':1},{'order':True},{'current_a':float('nan')}])
def test_invalid_quadrature_or_current(kwargs):
    with pytest.raises(ValueError):
        FiniteToroidPotential(**kwargs)


def test_external_domain_and_pair_overlap_are_explicit():
    field=FiniteToroidPotential()
    with pytest.raises(ValueError,match='inside a core'):
        field.vector(.1,0)
    with pytest.raises(ValueError,match='half-separation'):
        OpposedFiniteToroidsPotential(axial_offset=.001)
