"""Finite-core A audit; no prediction of an additional UBT detector voltage."""
import argparse
import json
import platform
from pathlib import Path
import sys
import numpy as np
import scipy
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from finite_toroid import FiniteToroidPotential, OpposedFiniteToroidsPotential, magnetic_energy_harmonics
from toroid_transfer import RectangularToroid
from toroidal_potential import ThinToroidPotential, diagnose, circulation


def report():
    geometry=RectangularToroid()
    pair=OpposedFiniteToroidsPotential(geometry,1,.05,10)
    records=[]
    for order in (4,8,12):
        potential=FiniteToroidPotential(geometry,1,.05,order)
        records.append({'order_per_dimension':order,'A_axis_T_m':potential.vector(0,.01)[1],
                        'relative_flux_error':abs(potential.quadrature_flux_wb/potential.flux_wb-1)})
    flux=pair.upper.flux_wb
    thin=ThinToroidPotential(flux,.1,.05)
    finite_axis=pair.upper.vector(0,.01)[1]
    diagnostics=diagnose(pair,.02,.01,step=2e-5)
    return {
        'schema':'finite-toroid-ubt-bridge/v1',
        'status':'MAXWELL_BASELINE_PLUS_CONDITIONAL_MATCHING_TARGET',
        'source_snapshots':{
            'antihemholz':'9b35e3facf93a763b5dc921f908e706cc7c2cb0e',
            'hyperspace_waves':'413e34412b4cae2d6481f8da564ef74a4ea1949d',
            'ubt':'386b6dbed7d35ea63c8fd42fd88d57ac2f9cdaa3'},
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'assumptions':['static 1 A current','linear continuous ideal toroid',
                       'air core mu_r=1','no winding/lead leakage','external points only',
                       'no cross-core inductance in energy example'],
        'flux_Wb':flux,'center_gradient_T':pair.center_gradient,
        'external_div_A_T':diagnostics.divergence,'external_curl_A_phi_T':diagnostics.curl_phi,
        'linked_circulation_Wb':circulation(pair,.1,.05,.02),
        'relative_finite_vs_filament_axis_correction':finite_axis/thin.vector(0,.01)[1]-1,
        'quadrature_convergence':records,
        'single_toroid_energy_harmonics_J':magnetic_energy_harmonics(geometry.primary_inductance_h,0,1),
        'opposed_disjoint_pair_energy_harmonics_J':magnetic_energy_harmonics(2*geometry.primary_inductance_h,0,1),
        'additional_G_voltage_V':0,
        'G_zero_scope':'complete conserved source, ordinary detector, retarded zero scalar data in existing gate',
        'physical_Theta_mode_source_detector':'UNDETERMINED, not assigned a fitted gain',
        'lean_status':'LEAN-PENDING: no formal PDE/kernel proof in this repository',
        'not_tested':['final UBT action or reduced Theta Hessian','nonzero extra physical scalar',
                      'time-dependent finite-core A','real lead fields','losses','hardware data']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=report()
    target=Path(__file__).resolve().parents[1]/'docs/results/finite-core-bridge.json'
    if args.check:
        old=json.loads(target.read_text())
        for key in ('status','source_snapshots','assumptions','G_zero_scope','not_tested'):
            assert old[key]==result[key],key
        for key in ('flux_Wb','center_gradient_T','linked_circulation_Wb',
                    'relative_finite_vs_filament_axis_correction'):
            np.testing.assert_allclose(result[key],old[key],rtol=2e-9,atol=1e-18)
        assert abs(result['external_div_A_T'])<1e-12
        assert abs(result['external_curl_A_phi_T'])<1e-12
        assert abs(result['linked_circulation_Wb']/result['flux_Wb']-1)<1e-9
        print('Finite-core Maxwell results agree; physical Theta coupling remains undetermined.')
    else:
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(json.dumps(result,indent=2)+'\n')
        print(target)


if __name__=='__main__':
    main()
