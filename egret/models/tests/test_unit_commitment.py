'''
unit commitment tester
'''
import json
import os
import math

import pytest
from pyomo.opt import SolverFactory, TerminationCondition
from egret.models.unit_commitment import *
from egret.data.model_data import ModelData

current_dir = os.path.dirname(os.path.abspath(__file__))

test_cases = [os.path.join(current_dir,'uc_test_instances', 'test_case_{}.json'.format(i)) for i in range(1,6)]
test_objvals = [4218219.415648284, 5476393.707647476, 6023692.988920635, 5484256.671628478, 6091360.072517988]

def _test_uc_model(uc_model):

    for test_case, ref_objval in zip(test_cases, test_objvals):
    
        md_dict = json.load(open(test_case,'r'))
        md = ModelData(md_dict)
        model = uc_model(md)
        opt = SolverFactory('cbc')
        result = opt.solve(model, tee=True)

        assert result.solver.termination_condition == TerminationCondition.optimal
        assert math.isclose(ref_objval, result.problem.upper_bound)

def test_tight_uc_model():
    _test_uc_model(create_tight_unit_commitment_model)

def test_compact_uc_model():
    _test_uc_model(create_compact_unit_commitment_model)
