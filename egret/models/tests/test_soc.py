#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
acopf tester
'''
import os
import math
import unittest
from pyomo.opt import SolverFactory, TerminationCondition
import pyomo.environ as pe
from egret.models.acopf import create_psv_acopf_model
from egret.models.ac_relaxations import create_soc_relaxation
from egret.data.model_data import ModelData
from parameterized import parameterized
from egret.parsers.matpower_parser import create_ModelData

current_dir = os.path.dirname(os.path.abspath(__file__))
case_names = ['pglib_opf_case3_lmbd',
              'pglib_opf_case5_pjm',
              'pglib_opf_case30_fsr',
              'pglib_opf_case30_ieee',
              'pglib_opf_case39_epri',
              'pglib_opf_case162_ieee_dtc',
              'pglib_opf_case179_goc',
              'pglib_opf_case240_pserc',
              'pglib_opf_case300_ieee']
test_cases = [os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', '{}.m'.format(i)) for i in case_names]
upper_bounds = {'pglib_opf_case3_lmbd': 5812.6,
                'pglib_opf_case5_pjm': 1.7552e+04,
                'pglib_opf_case30_fsr': 5.7577e+02,
                'pglib_opf_case30_ieee': 8.2085e+03,
                'pglib_opf_case39_epri': 1.3842e+05,
                'pglib_opf_case162_ieee_dtc': 1.0808e+05,
                'pglib_opf_case179_goc': 7.5427e+05,
                'pglib_opf_case240_pserc': 3.3297e+06,
                'pglib_opf_case300_ieee': 5.6522e+05}
gaps = {'pglib_opf_case3_lmbd': 1.32,
        'pglib_opf_case5_pjm': 14.55,
        'pglib_opf_case30_fsr': 0.39,
        'pglib_opf_case30_ieee': 18.84,
        'pglib_opf_case39_epri': 0.56,
        'pglib_opf_case162_ieee_dtc': 5.95,
        'pglib_opf_case179_goc': 0.16,
        'pglib_opf_case240_pserc': 2.78,
        'pglib_opf_case300_ieee': 2.63}


class TestSOC(unittest.TestCase):
    show_output = True

    @classmethod
    def setUpClass(self):
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib_opf import get_pglib_opf
            get_pglib_opf(download_dir)

    @parameterized.expand(case_names)
    def test_soc_model(self, case_name):
        test_case = os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', '{}.m'.format(case_name))
        upper_bound_soln = upper_bounds[case_name]
        gap_soln = gaps[case_name]
        md = create_ModelData(test_case)
        nlp, scaled_md = create_psv_acopf_model(md)
        rel, scaled_md = create_soc_relaxation(md)

        opt = SolverFactory('ipopt')
        opt.options['linear_solver'] = 'mumps'

        res = opt.solve(nlp, tee=False)
        self.assertTrue(res.solver.termination_condition == TerminationCondition.optimal)
        ub = pe.value(nlp.obj)

        res = opt.solve(rel, tee=False)
        self.assertTrue(res.solver.termination_condition == TerminationCondition.optimal)
        lb = pe.value(rel.obj)

        gap = (ub - lb) / ub * 100
        print(ub, upper_bound_soln, gap, gap_soln)
        comparison = math.isclose(ub, upper_bound_soln, rel_tol=1e-4)
        self.assertTrue(comparison)
        comparison = math.isclose(gap, gap_soln, abs_tol=1e-2)
        self.assertTrue(comparison)


if __name__ == '__main__':
    unittest.main()
