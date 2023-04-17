#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
unit commitment tester
'''
import json
import os
import math

import unittest
import pytest
from pyomo.opt import SolverFactory, TerminationCondition
from pyomo.core.plugins.transform.relax_integrality \
        import RelaxIntegrality
from pyomo.solvers.plugins.solvers.persistent_solver import PersistentSolver
from egret.models.unit_commitment import *
from egret.data.model_data import ModelData
from egret.common.solver_interface import _set_options

current_dir = os.path.dirname(os.path.abspath(__file__))
test_cases = [os.path.join(current_dir,'uc_test_instances', 'test_case_{}.json'.format(i)) for i in range(1,6)]
test_int_objvals = [4201915.017320504, 5454367.7670904165, 5999272.361123627, 5461120.3231092375, 6062406.32677043]

## for FP comparisons
rel_tol = 1e-08

solver_list = ['xpress_persistent', 'gurobi_persistent', 'cplex_persistent', 'gurobi', 'cplex']
test_solver = None
comm_mip_avail = False
for solver in solver_list:
    try:
        if SolverFactory(solver).available():
            test_solver = solver
            comm_mip_avail = True
            break
    except:
        continue
if test_solver is None: 
    for solver in ['cbc', 'glpk']:
        try:
            if SolverFactory(solver).available():
                test_solver = solver
                break
        except:
            continue
if test_solver is None:
    raise RuntimeError("No MIP/LP solver found for unit commitment tests")

def _test_uc_model(uc_model, relax=False, test_objvals=test_int_objvals):

    for test_case, ref_objval in zip(test_cases, test_objvals):
    
        md_dict = json.load(open(test_case,'r'))
        md = ModelData(md_dict)
        
        model = uc_model(md, relaxed=relax, slack_type=SlackType.NONE)
        opt = SolverFactory(test_solver)
        _set_options(opt, mipgap=0.0)

        if isinstance(opt, PersistentSolver):
            opt.set_instance(model)
        
        result = opt.solve(model, tee=False)

        assert result.solver.termination_condition == TerminationCondition.optimal
        assert math.isclose(ref_objval, result.problem.upper_bound, rel_tol=rel_tol)

def _make_get_dcopf_uc_model(network):
    def get_dcopf_uc_model(model_data, relaxed=False, **kwargs):
        return create_tight_unit_commitment_model(model_data,
                                network_constraints=network,
                                relaxed=relaxed,
                                **kwargs)
    return get_dcopf_uc_model

## definitely skip MIP tests if we don't have one of gurobi or cplex available
@unittest.skipUnless(comm_mip_avail, "Neither Gurobi or CPLEX solver is available")
def test_int_all_uc_models():
    _test_uc_model(create_tight_unit_commitment_model)
    _test_uc_model(create_compact_unit_commitment_model)
    _test_uc_model(create_KOW_unit_commitment_model)
    _test_uc_model(create_ALS_unit_commitment_model)
    _test_uc_model(create_MLR_unit_commitment_model)
    _test_uc_model(create_random1_unit_commitment_model)
    _test_uc_model(create_random2_unit_commitment_model)
    _test_uc_model(create_OAV_unit_commitment_model)
    _test_uc_model(create_OAV_tighter_unit_commitment_model)
    _test_uc_model(create_OAV_original_unit_commitment_model)
    _test_uc_model(create_OAV_up_downtime_unit_commitment_model)
    _test_uc_model(create_CA_unit_commitment_model)
    _test_uc_model(create_CHP_unit_commitment_model)
    _test_uc_model(create_super_tight_unit_commitment_model)

def test_tight_uc_model():
    lp_obj_list = [4194720.23424, 5441076.85034, 5988496.92621, 5453617.47912, 6055376.54656]
    _test_uc_model(create_tight_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_compact_uc_model():
    lp_obj_list = [4194304.94748, 5440720.727, 5988068.23178, 5453218.02764, 6055020.46427]
    _test_uc_model(create_compact_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_KOW_uc_model():
    lp_obj_list = [4193749.67682, 5440148.79074, 5987686.94763, 5452888.22712, 6054163.40576]
    _test_uc_model(create_KOW_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_ALS_uc_model():
    lp_obj_list = [4193603.40346, 5439977.63794, 5987392.27642, 5452580.38476, 6054545.74347]
    _test_uc_model(create_ALS_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_MLR_uc_model():
    lp_obj_list = [4193700.64155, 5440122.0449, 5987617.01183, 5452837.51833, 6054088.71399]
    _test_uc_model(create_MLR_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_random1_uc_model():
    lp_obj_list = [4194304.94748, 5440720.727, 5988068.23178, 5453218.02764, 6055020.46427]
    _test_uc_model(create_random1_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_random2_uc_model():
    lp_obj_list = [4194686.42109, 5441087.41223, 5988465.58558, 5453619.48855, 6055360.5608] 
    _test_uc_model(create_random2_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_OAV_uc_model():
    lp_obj_list = [4190770.57777, 5436680.81342, 5984071.37653, 5449824.53072, 6051451.70067]
    _test_uc_model(create_OAV_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_OAV_tighter_uc_model():
    lp_obj_list = [4190774.76258, 5436685.6315, 5984097.794, 5449825.81448, 6051485.86608]
    _test_uc_model(create_OAV_tighter_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_OAV_original_uc_model():
    lp_obj_list = [4186901.74384, 5428888.70061, 5975676.69077, 5443849.68783, 6041296.59018]
    _test_uc_model(create_OAV_original_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_OAV_up_downtime_uc_model():
    lp_obj_list = [4190745.01259, 5436634.52576, 5984052.06305, 5449795.75874, 6051432.92077]
    _test_uc_model(create_OAV_up_downtime_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_CA_uc_model():
    lp_obj_list = [4185855.30972, 5423650.80043, 5965411.93718, 5439434.94733, 6029118.03019]
    _test_uc_model(create_CA_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

# This test sometimes times out Travis, so skip if using cbc or glpk
@unittest.skipUnless(comm_mip_avail, "Neither Gurobi or CPLEX solver is available")
def test_CHP_uc_model():
    lp_obj_list = [4195062.031337061, 5441117.29937132, 5989325.835633724, 5454266.5169105325, 6055624.773215021]
    _test_uc_model(create_CHP_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_super_tight_uc_model():
    lp_obj_list = [4194720.234239616, 5441076.850336192, 5988496.926211054, 5453617.479118543, 6055376.546557944]
    _test_uc_model(create_super_tight_unit_commitment_model, relax=True, test_objvals=lp_obj_list)

def test_uc_runner():
    test_names = ['tiny_uc_{}'.format(i) for i in range(1,12+1)]
    for test_name in test_names:
        input_json_file_name = os.path.join(current_dir, 'uc_test_instances', test_name+'.json')
        md_in = ModelData(input_json_file_name)
        md_results = solve_unit_commitment(md_in, solver=test_solver, mipgap=0.0, slack_type=SlackType.BUS_BALANCE )

        reference_json_file_name = os.path.join(current_dir, 'uc_test_instances', test_name+'_results.json')
        md_reference = ModelData(reference_json_file_name)
        assert math.isclose(md_reference.data['system']['total_cost'], md_results.data['system']['total_cost'], rel_tol=rel_tol)

def test_uc_transmission_models():

    ## the network tests can optionally specify some kwargs so we can pass them into solve_unit_commitment
    tc_networks = {'btheta_power_flow': [dict()],
                   'ptdf_power_flow':[
                       {'ptdf_options': {'lazy':False}},
                       dict(),
                      ],
                   'power_balance_constraints':[dict()],
                   }
    no_network = 'copperplate_power_flow'
    test_names = ['tiny_uc_tc'] + ['tiny_uc_tc_{}'.format(i) for i in range(2,11+1)]
    ## based on tiny_uc, tiny_uc_tc_2 has an interface, tiny_uc_tc_3 has a relaxed interface, tiny_uc_tc_4 has a relaxed flow limit; tiny_uc_tc_7 has a HVDC line

    for test_name in test_names:
        input_json_file_name = os.path.join(current_dir, 'uc_test_instances', test_name+'.json')
        md_in = ModelData(input_json_file_name)
        for tc in tc_networks:
            for kwargs in tc_networks[tc]:

                md_results = solve_unit_commitment(md_in, solver=test_solver, mipgap=0.0,
                                                    slack_type=SlackType.TRANSMISSION_LIMITS, uc_model_generator = _make_get_dcopf_uc_model(tc), **kwargs)
                reference_json_file_name = os.path.join(current_dir, 'uc_test_instances', test_name+'_results.json')
                md_reference = ModelData(reference_json_file_name)
                assert math.isclose(md_reference.data['system']['total_cost'], md_results.data['system']['total_cost'], rel_tol=rel_tol)

    ## test copperplate
    test_name = 'tiny_uc_1'
    md_in = ModelData(os.path.join(current_dir, 'uc_test_instances', 'tiny_uc_tc_2.json'))
    md_results = solve_unit_commitment(md_in, solver=test_solver, mipgap=0.0, uc_model_generator = _make_get_dcopf_uc_model(no_network))
    reference_json_file_name = os.path.join(current_dir, 'uc_test_instances', test_name+'_results.json')
    md_reference = ModelData(reference_json_file_name)
    assert math.isclose(md_reference.data['system']['total_cost'], md_results.data['system']['total_cost'], rel_tol=rel_tol)

def test_uc_relaxation():
    test_name = 'tiny_uc_tc'
    input_json_file_name = os.path.join(current_dir, 'uc_test_instances', test_name+'.json')

    md_in = ModelData(input_json_file_name)

    md_results = solve_unit_commitment(md_in, solver=test_solver, slack_type=SlackType.TRANSMISSION_LIMITS, relaxed=True)
    reference_json_file_name = os.path.join(current_dir, 'uc_test_instances', test_name+'_relaxed_results.json')
    md_reference = ModelData(reference_json_file_name)
    assert math.isclose(md_reference.data['system']['total_cost'], md_results.data['system']['total_cost'], rel_tol=rel_tol)

def test_uc_lazy_ptdf_thresholding():
    test_name = 'tiny_uc_tc'
    input_json_file_name = os.path.join(current_dir, 'uc_test_instances', test_name+'.json')

    md_in = ModelData(input_json_file_name)

    tc_sol_fn = test_name + '_relaxed_results.json'
    ntc_sol_fn = test_name + '_relaxed_unconstrained_results.json'

    tc_cost = json.load(open(os.path.join(current_dir, 'uc_test_instances', tc_sol_fn), 'r'))['system']['total_cost']
    ntc_cost = json.load(open(os.path.join(current_dir, 'uc_test_instances', ntc_sol_fn), 'r'))['system']['total_cost']

    # NOTE: This test case has one congested branch, one end of which has bus_kv = 300, and
    #       the other end has bus_kv = 200, so these test cases attempt to capture all the
    #       logic of the thresholding
    ptdf_sol_options = [
                        (tc_cost,  {'branch_kv_threshold':100, 'kv_threshold_type':'one' }),
                        (tc_cost,  {'branch_kv_threshold':100, 'kv_threshold_type':'both'}),
                        (tc_cost,  {'branch_kv_threshold':200, 'kv_threshold_type':'one' }),
                        (tc_cost,  {'branch_kv_threshold':200, 'kv_threshold_type':'both'}),
                        (tc_cost,  {'branch_kv_threshold':201, 'kv_threshold_type':'one' }), 
                        (ntc_cost, {'branch_kv_threshold':201, 'kv_threshold_type':'both'}),
                        (tc_cost,  {'branch_kv_threshold':300, 'kv_threshold_type':'one' }), 
                        (ntc_cost, {'branch_kv_threshold':300, 'kv_threshold_type':'both'}),
                        (ntc_cost, {'branch_kv_threshold':301, 'kv_threshold_type':'one' }),
                        (ntc_cost, {'branch_kv_threshold':301, 'kv_threshold_type':'both'}),
                       ]
    for c, ptdf_opt in ptdf_sol_options:
        md_results = solve_unit_commitment(md_in, solver=test_solver, relaxed=True, slack_type=SlackType.TRANSMISSION_LIMITS, ptdf_options=ptdf_opt)
        assert math.isclose(c, md_results.data['system']['total_cost'], rel_tol=rel_tol)

def test_uc_ptdf_termination():
    test_name = 'tiny_uc_tc_3'
    input_json_file_name = os.path.join(current_dir, 'uc_test_instances', test_name+'.json')

    md_in = ModelData(input_json_file_name)

    kwargs = {'ptdf_options':{'lazy': True, 'rel_ptdf_tol':10.}}
    md_results, results = solve_unit_commitment(md_in, solver=test_solver, relaxed=True, slack_type=SlackType.TRANSMISSION_LIMITS, return_results=True, **kwargs)

    assert results.egret_metasolver['iterations'] == 1

def test_uc_tc_disconnected():
    test_name = 'tiny_uc_tc_disconnected'
    input_json_file_name = os.path.join(current_dir, 'uc_test_instances', test_name+'.json')

    md_in = ModelData(input_json_file_name)

    kwargs = {'network_constraints':'btheta_power_flow'}

    md_results, results = solve_unit_commitment(md_in, solver=test_solver, relaxed=True, slack_type=SlackType.TRANSMISSION_LIMITS, return_results=True, **kwargs)

    reference_json_file_name = os.path.join(current_dir, 'uc_test_instances', 'tiny_uc_tc_relaxed_results.json')
    md_reference = ModelData(reference_json_file_name)
    assert math.isclose(md_reference.data['system']['total_cost'], md_results.data['system']['total_cost'], rel_tol=rel_tol)

    kwargs = {'network_constraints':'power_balance_constraints'}

    md_results, results = solve_unit_commitment(md_in, solver=test_solver, relaxed=True, slack_type=SlackType.TRANSMISSION_LIMITS, return_results=True, **kwargs)

    reference_json_file_name = os.path.join(current_dir, 'uc_test_instances', 'tiny_uc_tc_relaxed_results.json')
    md_reference = ModelData(reference_json_file_name)
    assert math.isclose(md_reference.data['system']['total_cost'], md_results.data['system']['total_cost'], rel_tol=rel_tol)

    kwargs = {'network_constraints':'ptdf_power_flow'}

    with pytest.raises(RuntimeError, match="Network is not connected, cannot use PTDF formulation"):
        md_results, results = solve_unit_commitment(md_in, solver=test_solver, relaxed=True, slack_type=SlackType.TRANSMISSION_LIMITS, return_results=True, **kwargs)

def test_scuc():
    input_json_file_name = os.path.join(current_dir, 'uc_test_instances', 'test_scuc_masked.json')

    md_in = ModelData.read(input_json_file_name)
    md_results = solve_unit_commitment(md_in, solver=test_solver, slack_type=SlackType.TRANSMISSION_LIMITS, relaxed=True)

    md_baseline = ModelData.read(os.path.join(
                            current_dir, 'uc_test_instances', 'test_scuc_full_enforce_relaxed_sol.json'))
    assert math.isclose( md_results.data['system']['total_cost'], md_baseline.data['system']['total_cost'], rel_tol=rel_tol)

    ptdf_options = {'branch_kv_threshold':300, 'kv_threshold_type':'both'}
    md_results = solve_unit_commitment(md_in, solver=test_solver, mipgap=0.0, slack_type=SlackType.TRANSMISSION_LIMITS, ptdf_options=ptdf_options)
    md_baseline = ModelData.read(os.path.join(
                            current_dir, 'uc_test_instances', 'test_scuc_sparse_enforce_sol.json'))
    assert math.isclose( md_results.data['system']['total_cost'], md_baseline.data['system']['total_cost'], rel_tol=rel_tol)

