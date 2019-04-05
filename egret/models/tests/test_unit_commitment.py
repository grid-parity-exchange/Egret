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
test_objvals = [4218219.415648285, 5476393.707647482, 6023692.988920629, 5484256.685824468, 6091360.072517976]

def test_tight_uc_model():

    for test_case, ref_objval in zip(test_cases, test_objvals):
    
        md_dict = json.load(open(test_case,'r'))
        md = ModelData(md_dict)
        model = create_tight_unit_commitment_model(md)
        opt = SolverFactory('glpk')
        result = opt.solve(model, tee=False)

        assert result.solver.termination_condition == TerminationCondition.optimal
        assert math.isclose(ref_objval, model.TotalCostObjective())
