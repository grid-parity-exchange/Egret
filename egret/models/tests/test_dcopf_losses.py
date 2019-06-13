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
from egret.models.dcopf_losses import *
from egret.data.model_data import ModelData
from parameterized import parameterized
from egret.parsers.matpower_parser import create_ModelData

current_dir = os.path.dirname(os.path.abspath(__file__))
case_names = ['pglib_opf_case3_lmbd','pglib_opf_case30_ieee','pglib_opf_case300_ieee']
test_cases = [os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', '{}.m'.format(i)) for i in case_names]
btheta_soln_cases = [os.path.join(current_dir, 'transmission_test_instances', 'dcopf_losses_solution_files', '{}_btheta_dcopf_losses_solution.json'.format(i)) for i in case_names]
ptdf_soln_cases = [os.path.join(current_dir, 'transmission_test_instances', 'dcopf_losses_solution_files', '{}_ptdf_dcopf_losses_solution.json'.format(i)) for i in case_names]

class TestBThetaDCOPF(unittest.TestCase):
    show_output = True

    @classmethod
    def setUpClass(self):
        download_dir = os.path.join(current_dir, 'transmission_test_instances')
        if not os.path.exists(os.path.join(download_dir, 'pglib-opf-master')):
            from egret.thirdparty.get_pglib import get_pglib
            get_pglib(download_dir)

    @parameterized.expand(zip(test_cases, btheta_soln_cases))
    def test_btheta_losses_dcopf_model(self, test_case, soln_case, include_kwargs=False):
        dcopf_losses_model = create_btheta_losses_dcopf_model

        md_soln = ModelData()
        md_soln.read_from_json(soln_case)

        md_dict = create_ModelData(test_case)

        kwargs = {}
        if include_kwargs:
            kwargs = {'include_feasibility_slack': 'True'}
        md, results = solve_dcopf_losses(md_dict, "ipopt", dcopf_losses_model_generator=dcopf_losses_model, solver_tee=False, return_results=True, **kwargs)

        self.assertTrue(results.solver.termination_condition == TerminationCondition.optimal)
        comparison = math.isclose(md.data['system']['total_cost'], md_soln.data['system']['total_cost'], rel_tol=1e-4)
        self.assertTrue(comparison)


    @parameterized.expand(zip(test_cases, ptdf_soln_cases))
    def test_ptdf_losses_dcopf_model(self, test_case, soln_case, include_kwargs=False):
        dcopf_losses_model = create_ptdf_losses_dcopf_model

        md_soln = ModelData()
        md_soln.read_from_json(soln_case)

        md_dict = create_ModelData(test_case)

        from egret.models.acopf import solve_acopf
        md_dict, _, _ = solve_acopf(md_dict, "ipopt", solver_tee=False, return_model=True, return_results=True)

        kwargs = {}
        if include_kwargs:
            kwargs = {'include_feasibility_slack': 'True'}
        md, results = solve_dcopf_losses(md_dict, "ipopt", dcopf_losses_model_generator=dcopf_losses_model, solver_tee=False, return_results=True, **kwargs)

        self.assertTrue(results.solver.termination_condition == TerminationCondition.optimal)
        comparison = math.isclose(md.data['system']['total_cost'], md_soln.data['system']['total_cost'], rel_tol=1e-4)
        self.assertTrue(comparison)

if __name__ == '__main__':
     unittest.main()
