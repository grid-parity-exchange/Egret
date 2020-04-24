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
from egret.models.ac_relaxations import create_soc_relaxation, create_polar_acopf_relaxation, \
    create_rectangular_acopf_relaxation
from egret.data.model_data import ModelData
from parameterized import parameterized
from egret.parsers.matpower_parser import create_ModelData
try:
    import coramin
    coramin_available = True
except ImportError:
    coramin_available = False

current_dir = os.path.dirname(os.path.abspath(__file__))
case_names = ['pglib_opf_case3_lmbd',
              'pglib_opf_case5_pjm',
              'pglib_opf_case14_ieee',
              'pglib_opf_case30_fsr',
              'pglib_opf_case30_ieee',
              'pglib_opf_case39_epri',
              'pglib_opf_case57_ieee']
test_cases = [os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', '{}.m'.format(i)) for i in case_names]
upper_bounds = {'pglib_opf_case3_lmbd': 5812.6,
                'pglib_opf_case5_pjm': 1.7552e+04,
                'pglib_opf_case14_ieee': 2.1781e03,
                'pglib_opf_case30_fsr': 5.7577e+02,
                'pglib_opf_case30_ieee': 8.2085e+03,
                'pglib_opf_case39_epri': 1.3842e+05,
                'pglib_opf_case57_ieee': 3.7589e04}
gaps = {'pglib_opf_case3_lmbd': 1.32,
        'pglib_opf_case5_pjm': 14.55,
        'pglib_opf_case14_ieee': 0.11,
        'pglib_opf_case30_fsr': 0.39,
        'pglib_opf_case30_ieee': 18.84,
        'pglib_opf_case39_epri': 0.56,
        'pglib_opf_case57_ieee': 0.16}


def get_obj(m):
    obj = None
    for _obj in m.component_data_objects(pe.Objective, descend_into=True, active=True, sort=True):
        assert obj is None
        obj = _obj
    return obj


@unittest.skipIf(not coramin_available, "coramin is not available")
class TestRelaxations(unittest.TestCase):
    show_output = True

    def setUp(self):
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
        rel, scaled_md = create_soc_relaxation(md, use_linear_relaxation=False)

        opt = SolverFactory('ipopt')
        opt.options['linear_solver'] = 'mumps'

        res = opt.solve(nlp, tee=False)
        self.assertTrue(res.solver.termination_condition == TerminationCondition.optimal)
        ub = pe.value(get_obj(nlp))

        res = opt.solve(rel, tee=False)
        self.assertTrue(res.solver.termination_condition == TerminationCondition.optimal)
        lb = pe.value(get_obj(rel))

        gap = (ub - lb) / ub * 100
        comparison = math.isclose(ub, upper_bound_soln, rel_tol=1e-4)
        self.assertTrue(comparison)
        comparison = math.isclose(gap, gap_soln, abs_tol=1e-2)
        self.assertTrue(comparison)

    @parameterized.expand(case_names)
    def test_polar_relaxation(self, case_name):
        test_case = os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', '{}.m'.format(case_name))
        upper_bound_soln = upper_bounds[case_name]
        md = create_ModelData(test_case)
        nlp, scaled_md = create_polar_acopf_relaxation(md, include_soc=False, use_linear_relaxation=True)
        for b in coramin.relaxations.relaxation_data_objects(nlp, descend_into=True, active=True, sort=True):
            if isinstance(b, (coramin.relaxations.PWCosRelaxationData, coramin.relaxations.PWSinRelaxationData)):
                v = b.get_rhs_vars()[0]
                v.setlb(max(-math.pi / 2, v.lb))
                v.setub(min(math.pi / 2, v.ub))
            b.rebuild(build_nonlinear_constraint=True)

        opt = SolverFactory('ipopt')
        opt.options['linear_solver'] = 'mumps'

        res = opt.solve(nlp, tee=False)
        self.assertTrue(res.solver.termination_condition == TerminationCondition.optimal)
        ub = pe.value(get_obj(nlp))
        comparison = math.isclose(ub, upper_bound_soln, rel_tol=1e-4)
        self.assertTrue(comparison)

    @parameterized.expand(case_names)
    def test_rectangular_relaxation(self, case_name):
        test_case = os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', '{}.m'.format(case_name))
        upper_bound_soln = upper_bounds[case_name]
        md = create_ModelData(test_case)
        nlp, scaled_md = create_rectangular_acopf_relaxation(md, include_soc=False, use_linear_relaxation=True)
        for b in coramin.relaxations.relaxation_data_objects(nlp, descend_into=True, active=True, sort=True):
            b.rebuild(build_nonlinear_constraint=True)

        opt = SolverFactory('ipopt')
        opt.options['linear_solver'] = 'mumps'

        res = opt.solve(nlp, tee=False)
        self.assertTrue(res.solver.termination_condition == TerminationCondition.optimal)
        ub = pe.value(get_obj(nlp))
        comparison = math.isclose(ub, upper_bound_soln, rel_tol=1e-4)
        self.assertTrue(comparison)


if __name__ == '__main__':
    unittest.main()
