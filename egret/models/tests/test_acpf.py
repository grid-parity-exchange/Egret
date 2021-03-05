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
import egret.model_library.transmission.tx_utils as tx_utils
from egret.data.model_data import ModelData
from parameterized import parameterized
from egret.parsers.matpower_parser import create_ModelData

current_dir = os.path.dirname(os.path.abspath(__file__))
#case_names = ['pglib_opf_case3_lmbd','pglib_opf_case30_ieee','pglib_opf_case300_ieee',]
case_names = ['pglib_opf_case14_ieee', 'pglib_opf_case30_ieee', 'pglib_opf_case118_ieee', 'pglib_opf_case3012wp_k']
test_cases = [os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', '{}.m'.format(i)) for i in case_names]
soln_cases = [os.path.join(current_dir, 'transmission_test_instances', 'acpf_solution_files', '{}_acpf_solution.json'.format(i)) for i in case_names]

class TestPSVACPF(unittest.TestCase):
    show_output = True

    @classmethod
    def setUpClass(self):
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib_opf import get_pglib_opf
            get_pglib_opf(download_dir)

    @parameterized.expand(zip(test_cases, soln_cases))
    def test_acpf_model(self, test_case, soln_case, include_kwargs=False):
        # md_soln = ModelData.read(soln_case)
        md_dict = create_ModelData(test_case)
        kwargs = {}
        if include_kwargs:
            kwargs = {'include_feasibility_slack':True}
        md, m, results = acpf.solve_acpf(md_dict, "ipopt",
                                      acpf_model_generator=acpf.create_psv_acpf_model,
                                      solver_tee=True, return_results=True, return_model=True, **kwargs)

        self.assertTrue(results.solver.termination_condition == TerminationCondition.optimal)        
        # check the solution
        with open(soln_case, 'r') as fd:
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
            self.assertAlmostEqual(1e-3*pg1, 1e-3*pg2, places=4)
            self.assertAlmostEqual(1e-3*qg1, 1e-3*qg2, places=4)

        bus_attrs = md.attributes(element_type='bus')
        for k in bus_attrs['vm'].keys():
            self.assertAlmostEqual(bus_attrs['vm'][k], soln['vm'][k], places=3)
            self.assertAlmostEqual(bus_attrs['va'][k], soln['va'][k], places=3)

if __name__ == '__main__':
    unittest.main()
