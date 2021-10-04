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
import json
import pyomo.environ as pyo
from pyomo.opt import SolverFactory, TerminationCondition
import egret.models.acpf as acpf
import egret.models.acopf as acopf
import egret.model_library.transmission.tx_utils as tx_utils
from egret.data.model_data import ModelData
from parameterized import parameterized
from egret.parsers.matpower_parser import create_ModelData

current_dir = os.path.dirname(os.path.abspath(__file__))
case_names = ['pglib_opf_case14_ieee', 'pglib_opf_case30_ieee', 'pglib_opf_case118_ieee']
test_cases = [os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', '{}.m'.format(i)) for i in case_names]

def _get_power_and_voltage_from_md(md):
    d = dict()
    gen_attrs = md.attributes(element_type='generator')
    d['pg'] = {k:v for k,v in gen_attrs['pg'].items()}
    d['qg'] = {k:v for k,v in gen_attrs['qg'].items()}
    bus_attrs = md.attributes(element_type='bus')
    d['vm'] = {k:v for k,v in bus_attrs['vm'].items()}
    d['va'] = {k:v for k,v in bus_attrs['va'].items()}
    return d

class TestRepeatedSolves(unittest.TestCase):
    show_output = True

    @classmethod
    def setUpClass(self):
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib_opf import get_pglib_opf
            get_pglib_opf(download_dir)

    @parameterized.expand(test_cases)
    def test_acpf_model(self, test_case):
        md_orig = create_ModelData(test_case)
        md_orig.data['system']['reference_bus_angle'] = 30

        md_before_solve = md_orig.clone()
        md_after_acpf, m, results = acpf.solve_acpf(md_before_solve, "ipopt",
                             acpf_model_generator=acpf.create_psv_acpf_model,
                             solver_tee=False, return_results=True, return_model=True)
        self.assertTrue(results.solver.termination_condition \
                        == TerminationCondition.optimal)

        # get the soln for pg, qg, vm, va
        d1 = _get_power_and_voltage_from_md(md_after_acpf)

        # now solve it again using the model_data solution and see if we get the same answer
        md_after_acpf2, m, results = acpf.solve_acpf(md_after_acpf, "ipopt",
                             acpf_model_generator=acpf.create_psv_acpf_model,
                             solver_tee=False, return_results=True, return_model=True)
        self.assertTrue(results.solver.termination_condition \
                        == TerminationCondition.optimal)

        d2 = _get_power_and_voltage_from_md(md_after_acpf2)

        # compare the results
        for k in d1['pg']:
            self.assertAlmostEqual(1e-3*d1['pg'][k], 1e-3*d2['pg'][k], places=4)
            self.assertAlmostEqual(1e-3*d1['qg'][k], 1e-3*d2['qg'][k], places=4)
        for k in d1['vm']:
            self.assertAlmostEqual(1e-3*d1['vm'][k], 1e-3*d2['vm'][k], places=4)
            self.assertAlmostEqual(1e-3*d1['va'][k], 1e-3*d2['va'][k], places=4)

    @parameterized.expand(test_cases)
    def test_acopf_acpf(self, test_case):
        md_orig = create_ModelData(test_case)
        md_orig.data['system']['reference_bus_angle'] = 30

        # solve the acopf problem
        md_before_solve = md_orig.clone()
        md_after_acopf, m, results = acopf.solve_acopf(md_before_solve, "ipopt",
                             acopf_model_generator=acopf.create_psv_acopf_model,
                             solver_tee=False, return_results=True, return_model=True)
        self.assertTrue(results.solver.termination_condition \
                        == TerminationCondition.optimal)

        # get the soln for pg, qg, vm, va
        d1 = _get_power_and_voltage_from_md(md_after_acopf)
        
        # now solve the acpf problem using the model_data solution and see if we get the same answer
        md_after_acpf, m, results = acpf.solve_acpf(md_after_acopf, "ipopt",
                             acpf_model_generator=acpf.create_psv_acpf_model,
                             solver_tee=False, return_results=True, return_model=True)
        self.assertTrue(results.solver.termination_condition \
                        == TerminationCondition.optimal)

        d2 = _get_power_and_voltage_from_md(md_after_acpf)

        # compare the results
        for k in d1['pg']:
            self.assertAlmostEqual(1e-3*d1['pg'][k], 1e-3*d2['pg'][k], places=4)
            self.assertAlmostEqual(1e-3*d1['qg'][k], 1e-3*d2['qg'][k], places=4)
        for k in d1['vm']:
            self.assertAlmostEqual(1e-3*d1['vm'][k], 1e-3*d2['vm'][k], places=4)
            self.assertAlmostEqual(1e-3*d1['va'][k], 1e-3*d2['va'][k], places=4)

if __name__ == '__main__':
    unittest.main()
