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
from egret.models.acopf import *
from egret.data.model_data import ModelData
from parameterized import parameterized
from egret.parsers.matpower_parser import create_ModelData
import numpy.testing as npt
import json

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


if __name__ == '__main__':
    unittest.main()
