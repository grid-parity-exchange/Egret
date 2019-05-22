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

current_dir = os.path.dirname(os.path.abspath(__file__))
case_names = ['pglib_opf_case3_lmbd','pglib_opf_case30_ieee','pglib_opf_case300_ieee','pglib_opf_case3012wp_k','pglib_opf_case13659_pegase']
test_cases = [os.path.join(current_dir, 'transmission_test_instances', '{}.json'.format(i)) for i in case_names]
soln_cases = [os.path.join(current_dir, 'transmission_test_instances', 'acopf_solution_files', '{}_acopf_solution.json'.format(i)) for i in case_names]


def _test_acopf_model(acopf_model):
    for test_case, soln_case in zip(test_cases, soln_cases):

        md_soln = ModelData()
        md_soln.read_from_json(soln_case)

        md_dict = ModelData()
        md_dict.read_from_json(test_case)

        md, results = solve_acopf(md_dict, "ipopt", acopf_model_generator=acopf_model, solver_tee=False, return_results=True)

        assert results.solver.termination_condition == TerminationCondition.optimal
        assert math.isclose(md.data['system']['total_cost'], md_soln.data['system']['total_cost'], rel_tol=1e04)


@unittest.skipUnless(SolverFactory('ipopt').available(), "Ipopt solver is not available")


def test_all_acopf_models():
    _test_acopf_model(create_psv_acopf_model)
    _test_acopf_model(create_rsv_acopf_model)
    _test_acopf_model(create_riv_acopf_model)


def test_psv_acopf_model():
    _test_acopf_model(create_psv_acopf_model)


def test_rsv_acopf_model():
    _test_acopf_model(create_rsv_acopf_model)


def test_riv_acopf_model():
    _test_acopf_model(create_riv_acopf_model)


# if __name__ == '__main__':
#     test_all_acopf_models()
