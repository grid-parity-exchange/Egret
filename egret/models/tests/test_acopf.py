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
from egret.models.acopf import *
from egret.data.model_data import ModelData
from parameterized import parameterized
from egret.parsers.matpower_parser import create_ModelData
import numpy.testing as npt
import json
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
case_names = ['pglib_opf_case3_lmbd','pglib_opf_case30_ieee','pglib_opf_case300_ieee','pglib_opf_case3012wp_k']
test_cases = [os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', '{}.m'.format(i)) for i in case_names]
soln_cases = [os.path.join(current_dir, 'transmission_test_instances', 'acopf_solution_files', '{}_acopf_solution.json'.format(i)) for i in case_names]
p_and_v_soln_cases = [os.path.join(current_dir, 'transmission_test_instances', 'acopf_solution_files', '{}_acopf_solution_p_and_v.json'.format(i)) for i in case_names]

def _test_p_and_v(tst, json_fname, md):
    with open(json_fname, 'r') as fd:
        soln = json.load(fd)

    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    gen_attrs = md.attributes(element_type='generator')
    bus_attrs = md.attributes(element_type='bus')
    gens_by_bus = tx_utils.gens_by_bus(buses, gens)

    for b in gens_by_bus:
        genlist = gens_by_bus[b]
        pg1 = 0
        pg2 = 0
        qg1 = 0
        qg2 = 0
        for g in genlist:
            pg1 += gen_attrs['pg'][g]
            pg2 += soln['pg'][g]
            qg1 += gen_attrs['qg'][g]
            qg2 += soln['qg'][g]
        npt.assert_allclose(pg1, pg2, rtol=1e-6, atol=1e-2)
        npt.assert_allclose(qg1, qg2, rtol=5e-2, atol=1e-2)

    bus_attrs = md.attributes(element_type='bus')
    for k in bus_attrs['vm'].keys():
        npt.assert_allclose(bus_attrs['vm'][k], soln['vm'][k], rtol=1e-3, atol=1e-3)
        npt.assert_allclose(bus_attrs['va'][k], soln['va'][k], rtol=1e-3, atol=1e-1)

class TestPSVACOPF(unittest.TestCase):
    show_output = True

    @classmethod
    def setUpClass(self):
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib_opf import get_pglib_opf
            get_pglib_opf(download_dir)

    @parameterized.expand(zip(test_cases, soln_cases, p_and_v_soln_cases))
    def test_acopf_model(self, test_case, soln_case, p_and_v_soln_case, include_kwargs=False):
        acopf_model = create_psv_acopf_model

        md_soln = ModelData.read(soln_case)
        md_dict = create_ModelData(test_case)

        kwargs = {}
        if include_kwargs:
            kwargs = {'include_feasibility_slack':True}
        md, results = solve_acopf(md_dict, "ipopt", acopf_model_generator=acopf_model, solver_tee=False, return_results=True, **kwargs)

        self.assertTrue(results.solver.termination_condition == TerminationCondition.optimal)
        comparison = math.isclose(md.data['system']['total_cost'], md_soln.data['system']['total_cost'], rel_tol=1e-6)
        self.assertTrue(comparison)
        _test_p_and_v(self, p_and_v_soln_case, md)

    def test_keep_vars(self):
        fname = os.path.join(current_dir, 'transmission_test_instances/pglib-opf-master/pglib_opf_case5_pjm.m')
        md = ModelData.read(fname)
        md.data["elements"]["generator"]["1"]["in_service"] = False
        md.data["elements"]["branch"]["2"]["in_service"] = False

        m1, _ = create_psv_acopf_model(md, keep_vars_for_out_of_service_elements=False)
        m2, _ = create_psv_acopf_model(md, keep_vars_for_out_of_service_elements=True)

        opt = SolverFactory('ipopt')
        res1 = opt.solve(m1)
        res2 = opt.solve(m2)

        self.assertEqual(res1.solver.termination_condition, TerminationCondition.optimal)
        self.assertEqual(res2.solver.termination_condition, TerminationCondition.optimal)

        obj1 = pe.value(m1.obj)
        obj2 = pe.value(m2.obj)

        self.assertAlmostEqual(obj1, obj2)
        self.assertTrue(m2.pg["1"].fixed)

        
class TestArctanACOPF(unittest.TestCase):
    show_output = True

    @classmethod
    def setUpClass(self):
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib_opf import get_pglib_opf
            get_pglib_opf(download_dir)

    @parameterized.expand(zip(test_cases, soln_cases))
    def test_acopf_model(self, test_case, soln_case):
        md_soln = ModelData.read(soln_case)
        md_dict = create_ModelData(test_case)

        model, scaled_md = create_atan_acopf_model(md_dict)
        opt = pe.SolverFactory('ipopt')
        res = opt.solve(model)

        self.assertTrue(res.solver.termination_condition == TerminationCondition.optimal)
        self.assertAlmostEqual(pe.value(model.obj)/md_soln.data['system']['total_cost'], 1, 4)

    def test_keep_vars(self):
        fname = os.path.join(current_dir, 'transmission_test_instances/pglib-opf-master/pglib_opf_case5_pjm.m')
        md = ModelData.read(fname)
        md.data["elements"]["generator"]["1"]["in_service"] = False
        md.data["elements"]["branch"]["2"]["in_service"] = False

        m1, _ = create_atan_acopf_model(md, keep_vars_for_out_of_service_elements=False)
        m2, _ = create_atan_acopf_model(md, keep_vars_for_out_of_service_elements=True)

        opt = SolverFactory('ipopt')
        res1 = opt.solve(m1)
        res2 = opt.solve(m2)

        self.assertEqual(res1.solver.termination_condition, TerminationCondition.optimal)
        self.assertEqual(res2.solver.termination_condition, TerminationCondition.optimal)

        obj1 = pe.value(m1.obj)
        obj2 = pe.value(m2.obj)

        self.assertAlmostEqual(obj1, obj2)
        self.assertTrue(m2.pg["1"].fixed)

class TestRSVACOPF(unittest.TestCase):
    show_output = True

    @classmethod
    def setUpClass(self):
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib_opf import get_pglib_opf
            get_pglib_opf(download_dir)

    @parameterized.expand(zip(test_cases, soln_cases, p_and_v_soln_cases))
    def test_acopf_model(self, test_case, soln_case, p_and_v_soln_case, include_kwargs=False):
        acopf_model = create_rsv_acopf_model

        md_soln = ModelData.read(soln_case)

        md_dict = create_ModelData(test_case)

        kwargs = {}
        if include_kwargs:
            kwargs = {'include_feasibility_slack':True}
        md, results = solve_acopf(md_dict, "ipopt", acopf_model_generator=acopf_model, solver_tee=False, return_results=True, **kwargs)

        self.assertTrue(results.solver.termination_condition == TerminationCondition.optimal)
        comparison = math.isclose(md.data['system']['total_cost'], md_soln.data['system']['total_cost'], rel_tol=1e-6)
        self.assertTrue(comparison)
        _test_p_and_v(self, p_and_v_soln_case, md)

    def test_keep_vars(self):
        fname = os.path.join(current_dir, 'transmission_test_instances/pglib-opf-master/pglib_opf_case5_pjm.m')
        md = ModelData.read(fname)
        md.data["elements"]["generator"]["1"]["in_service"] = False
        md.data["elements"]["branch"]["2"]["in_service"] = False

        m1, _ = create_rsv_acopf_model(md, keep_vars_for_out_of_service_elements=False)
        m2, _ = create_rsv_acopf_model(md, keep_vars_for_out_of_service_elements=True)

        opt = SolverFactory('ipopt')
        res1 = opt.solve(m1)
        res2 = opt.solve(m2)

        self.assertEqual(res1.solver.termination_condition, TerminationCondition.optimal)
        self.assertEqual(res2.solver.termination_condition, TerminationCondition.optimal)

        obj1 = pe.value(m1.obj)
        obj2 = pe.value(m2.obj)

        self.assertAlmostEqual(obj1, obj2)
        self.assertTrue(m2.pg["1"].fixed)

        
class TestRIVACOPF(unittest.TestCase):
    show_output = True

    @classmethod
    def setUpClass(self):
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib_opf import get_pglib_opf
            get_pglib_opf(download_dir)

    @parameterized.expand(zip(test_cases, soln_cases, p_and_v_soln_cases))
    def test_acopf_model(self, test_case, soln_case, p_and_v_soln_case, include_kwargs=False):
        acopf_model = create_riv_acopf_model

        md_soln = ModelData.read(soln_case)

        md_dict = create_ModelData(test_case)

        kwargs = {}
        if include_kwargs:
            kwargs = {'include_feasibility_slack':True}
        md, results = solve_acopf(md_dict, "ipopt", acopf_model_generator=acopf_model, solver_tee=False, return_results=True, **kwargs)

        self.assertTrue(results.solver.termination_condition == TerminationCondition.optimal)
        comparison = math.isclose(md.data['system']['total_cost'], md_soln.data['system']['total_cost'], rel_tol=1e-6)
        self.assertTrue(comparison)
        _test_p_and_v(self, p_and_v_soln_case, md)


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


class TestPWCost(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib_opf import get_pglib_opf
            get_pglib_opf(download_dir)

    def test_case30_pw_cost(self):
        original_md = ModelData.read(os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', 'pglib_opf_case30_as.m'))

        md = original_md.clone()
        m_poly, scaled_md = create_atan_acopf_model(md)
        opt = pe.SolverFactory('ipopt')
        res = opt.solve(m_poly)
        poly_obj = pe.value(m_poly.obj.expr)
        self.assertAlmostEqual(poly_obj, 803.12692652061719, places=2)

        poly_cost_to_pw_cost(md, num_points=2)
        m_pw, scaled_md = create_atan_acopf_model(md)
        res = opt.solve(m_pw)
        pw_obj = pe.value(m_pw.obj.expr)
        self.assertAlmostEqual(pw_obj, 827.62708294193681, places=2)

        md = original_md.clone()
        poly_cost_to_pw_cost(md, num_points=10)
        m_pw, scaled_md = create_atan_acopf_model(md)
        res = opt.solve(m_pw)
        pw_obj = pe.value(m_pw.obj.expr)
        self.assertAlmostEqual(pw_obj, 803.56080829371604, places=2)

    def test_pw_cost_with_out_of_service_gens(self):
        md = ModelData.read(os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', 'pglib_opf_case30_as.m'))
        poly_cost_to_pw_cost(md, num_points=3)
        md.data["elements"]["generator"]["2"]["in_service"] = False
        m1, _ = create_atan_acopf_model(md, keep_vars_for_out_of_service_elements=False)
        m2, _ = create_atan_acopf_model(md, keep_vars_for_out_of_service_elements=True)

        opt = pe.SolverFactory('ipopt')
        res1 = opt.solve(m1)
        res2 = opt.solve(m2)

        self.assertEqual(res1.solver.termination_condition, TerminationCondition.optimal)
        self.assertEqual(res2.solver.termination_condition, TerminationCondition.optimal)

        obj1 = pe.value(m1.obj)
        obj2 = pe.value(m2.obj)

        self.assertAlmostEqual(obj1, obj2)

    def test_pw_cost_with_out_of_service_gens_epi(self):
        md = ModelData.read(os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', 'pglib_opf_case30_as.m'))
        poly_cost_to_pw_cost(md, num_points=3)
        md.data["elements"]["generator"]["2"]["in_service"] = False
        m1, _ = create_atan_acopf_model(md, keep_vars_for_out_of_service_elements=False, pw_cost_model='epi')
        m2, _ = create_atan_acopf_model(md, keep_vars_for_out_of_service_elements=True, pw_cost_model='epi')

        opt = pe.SolverFactory('ipopt')
        res1 = opt.solve(m1)
        res2 = opt.solve(m2)

        self.assertEqual(res1.solver.termination_condition, TerminationCondition.optimal)
        self.assertEqual(res2.solver.termination_condition, TerminationCondition.optimal)

        obj1 = pe.value(m1.obj)
        obj2 = pe.value(m2.obj)

        self.assertAlmostEqual(obj1, obj2)


class TestDeltaThetaBounds(unittest.TestCase):
    def test_0_delta_theta_bounds(self):
        md = ModelData.read(os.path.join(current_dir, 'transmission_test_instances', 'test_instances', 'delta_theta_bounds_0.m'))
        m, _ = create_psv_acopf_model(md)
        self.assertAlmostEqual(m.s['1','2'].lb, -1.21)
        self.assertAlmostEqual(m.s['1','2'].ub, 1.21)
        self.assertAlmostEqual(m.c['1','2'].lb, -1.21)
        self.assertAlmostEqual(m.c['1','2'].ub, 1.21)

    def test_minus_10_to_10_delta_theta_bounds(self):
        md = ModelData.read(os.path.join(current_dir, 'transmission_test_instances', 'test_instances', 'delta_theta_bounds_minus_10_to_10.m'))
        m, _ = create_psv_acopf_model(md)
        self.assertAlmostEqual(m.s['1','2'].lb, 1.21 * math.sin(math.radians(-10)))
        self.assertAlmostEqual(m.s['1','2'].ub, 1.21 * math.sin(math.radians(10)))
        self.assertAlmostEqual(m.c['1','2'].lb, 0.81 * math.cos(math.radians(-10)))
        self.assertAlmostEqual(m.c['1','2'].ub, 1.21)

    def test_10_to_20_delta_theta_bounds(self):
        md = ModelData.read(os.path.join(current_dir, 'transmission_test_instances', 'test_instances', 'delta_theta_bounds_10_to_20.m'))
        m, _ = create_psv_acopf_model(md)
        self.assertAlmostEqual(m.s['1','2'].lb, 0.81 * math.sin(math.radians(10)))
        self.assertAlmostEqual(m.s['1','2'].ub, 1.21 * math.sin(math.radians(20)))
        self.assertAlmostEqual(m.c['1','2'].ub, 1.21 * math.cos(math.radians(10)))
        self.assertAlmostEqual(m.c['1','2'].lb, 0.81 * math.cos(math.radians(20)))


if __name__ == '__main__':
    unittest.main()
