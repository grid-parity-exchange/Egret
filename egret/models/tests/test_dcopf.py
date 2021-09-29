#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
dcopf tester
'''
import os
import math
import unittest
from pyomo.opt import SolverFactory, TerminationCondition
from egret.models.dcopf import *
from egret.data.model_data import ModelData
from parameterized import parameterized
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
case_names = ['pglib_opf_case3_lmbd','pglib_opf_case30_ieee','pglib_opf_case300_ieee']
test_cases = [os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', '{}.m'.format(i)) for i in case_names]
soln_cases = [os.path.join(current_dir, 'transmission_test_instances', 'dcopf_solution_files', '{}_dcopf_solution.json'.format(i)) for i in case_names]

builtin_case_names = ['hvdc_test_case3']
test_cases.extend(os.path.join(current_dir, 'transmission_test_instances', 'test_instances', '{}.json'.format(i)) for i in builtin_case_names)
soln_cases.extend(os.path.join(current_dir, 'transmission_test_instances', 'dcopf_solution_files', '{}_dcopf_solution.json'.format(i)) for i in builtin_case_names)

class TestDCOPF(unittest.TestCase):
    show_output = True

    @classmethod
    def setUpClass(self):
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib_opf import get_pglib_opf
            get_pglib_opf(download_dir)

    @parameterized.expand(zip(test_cases, soln_cases))
    def test_btheta_dcopf_model(self, test_case, soln_case, include_kwargs=False):
        dcopf_model = create_btheta_dcopf_model

        md_soln = ModelData.read(soln_case)

        md_dict = ModelData.read(test_case)

        kwargs = {}
        if include_kwargs:
            kwargs = {'include_feasibility_slack':True}
        md, results = solve_dcopf(md_dict, "ipopt", dcopf_model_generator=dcopf_model, solver_tee=False, return_results=True, **kwargs)

        self.assertTrue(results.solver.termination_condition == TerminationCondition.optimal)
        comparison = math.isclose(md.data['system']['total_cost'], md_soln.data['system']['total_cost'], rel_tol=1e-6)
        self.assertTrue(comparison)


    @parameterized.expand(zip(test_cases, soln_cases))
    def test_ptdf_dcopf_model(self, test_case, soln_case, include_kwargs=False):
        dcopf_model = create_ptdf_dcopf_model

        md_soln = ModelData.read(soln_case)

        md_dict = ModelData.read(test_case)

        kwargs = {}
        if include_kwargs:
            kwargs = {'include_feasibility_slack':True}
        md, results = solve_dcopf(md_dict, "ipopt", dcopf_model_generator=dcopf_model, solver_tee=False, return_results=True, **kwargs)

        self.assertTrue(results.solver.termination_condition == TerminationCondition.optimal)
        comparison = math.isclose(md.data['system']['total_cost'], md_soln.data['system']['total_cost'], rel_tol=1e-6)
        self.assertTrue(comparison)


def poly_cost_to_pw_cost(md: ModelData, num_points=10):
    gen_attrs = md.attributes(element_type='generator')
    p_cost = gen_attrs['p_cost']
    for gen_name, cost_dict in p_cost.items():
        assert cost_dict['cost_curve_type'] == 'polynomial'
        cost_dict['cost_curve_type'] = 'piecewise'
        poly_coefs = cost_dict['values']
        pw_values = list()
        p_min = gen_attrs['p_min'][gen_name]
        p_max = gen_attrs['p_max'][gen_name]
        for pt in np.linspace(p_min, p_max, num_points):
            cost = sum(coef*pt**exponent for exponent, coef in poly_coefs.items())
            pw_values.append((pt, cost))
        cost_dict['values'] = pw_values


class TestPWCostBTheta(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib_opf import get_pglib_opf
            get_pglib_opf(download_dir)

    def test_case30_pw_cost_b_theta(self):
        original_md = ModelData.read(os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', 'pglib_opf_case30_as.m'))

        md = original_md.clone()
        m_poly, scaled_md = create_btheta_dcopf_model(md)
        opt = pe.SolverFactory('ipopt')
        res = opt.solve(m_poly)
        poly_obj = pe.value(m_poly.obj.expr)
        self.assertAlmostEqual(poly_obj, 767.60209944437236, places=2)

        poly_cost_to_pw_cost(md, num_points=2)
        m_pw, scaled_md = create_btheta_dcopf_model(md)
        res = opt.solve(m_pw)
        pw_obj = pe.value(m_pw.obj.expr)
        self.assertAlmostEqual(pw_obj, 788.18275424543515, places=2)

        md = original_md.clone()
        poly_cost_to_pw_cost(md, num_points=10)
        m_pw, scaled_md = create_btheta_dcopf_model(md)
        res = opt.solve(m_pw)
        pw_obj = pe.value(m_pw.obj.expr)
        self.assertAlmostEqual(pw_obj, 767.74029197558707, places=2)


class TestPWCostPTDF(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib_opf import get_pglib_opf
            get_pglib_opf(download_dir)

    def test_case30_pw_cost_b_theta(self):
        original_md = ModelData.read(os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', 'pglib_opf_case30_as.m'))

        md = original_md.clone()
        m_poly, scaled_md = create_ptdf_dcopf_model(md)
        opt = pe.SolverFactory('ipopt')
        res = opt.solve(m_poly, tee=False)
        poly_obj = pe.value(m_poly.obj.expr)
        self.assertAlmostEqual(poly_obj, 767.60209944437236, places=2)

        poly_cost_to_pw_cost(md, num_points=2)
        m_pw, scaled_md = create_ptdf_dcopf_model(md)
        res = opt.solve(m_pw, tee=False)
        pw_obj = pe.value(m_pw.obj.expr)
        self.assertAlmostEqual(pw_obj, 783.89649629510222, places=2)

        md = original_md.clone()
        poly_cost_to_pw_cost(md, num_points=10)
        m_pw, scaled_md = create_ptdf_dcopf_model(md)
        res = opt.solve(m_pw, tee=False)
        pw_obj = pe.value(m_pw.obj.expr)
        self.assertAlmostEqual(pw_obj, 767.74029197558707, places=2)


if __name__ == '__main__':
     unittest.main()
