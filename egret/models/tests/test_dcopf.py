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

current_dir = os.path.dirname(os.path.abspath(__file__))
case_names = ['pglib_opf_case3_lmbd','pglib_opf_case30_ieee','pglib_opf_case300_ieee']#,'pglib_opf_case3012wp_k','pglib_opf_case13659_pegase']
test_cases = [os.path.join(current_dir, 'transmission_test_instances', '{}.json'.format(i)) for i in case_names]
soln_cases = [os.path.join(current_dir, 'transmission_test_instances', 'dcopf_solution_files', '{}_dcopf_solution.json'.format(i)) for i in case_names]

class TestBThetaDCOPF(unittest.TestCase):
    show_output = True

    @parameterized.expand(zip(test_cases, soln_cases))
    def test_btheta_dcopf_model(self, test_case, soln_case):
        dcopf_model = create_btheta_dcopf_model

        md_soln = ModelData()
        md_soln.read_from_json(soln_case)

        md_dict = ModelData()
        md_dict.read_from_json(test_case)

        md, results = solve_dcopf(md_dict, "gurobi", dcopf_model_generator=dcopf_model, solver_tee=False, return_results=True)

        self.assertTrue(results.solver.termination_condition == TerminationCondition.optimal)
        comparison = math.isclose(md.data['system']['total_cost'], md_soln.data['system']['total_cost'], rel_tol=1e04)
        self.assertTrue(comparison)


if __name__ == '__main__':
     unittest.main()