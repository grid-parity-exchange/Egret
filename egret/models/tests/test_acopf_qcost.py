#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
Piecewise reactive power costs in the ACOPF cost objective.

`gen.declare_expression_qg_operating_cost` expresses a piecewise reactive
power cost in terms of the `delta_qg` step variables.  It used to reference
`delta_pg` instead, and it also lacked the `len(cleaned_values) > 1` guard
that the real-power counterpart has.  `delta_pg` is only declared when some
generator has a *piecewise real power* cost, so a generator whose real power
cost is polynomial while its reactive power cost is piecewise made model
construction fail with

    AttributeError: 'ConcreteModel' object has no attribute 'delta_pg'

The MATPOWER parser emits `q_cost` whenever the gencost matrix carries
reactive power rows, so this combination is reachable from a normal input
file.  ACOPF is a nonlinear model, so these tests only build it; they do not
require a solver.
'''
import copy
import os

import pyomo.environ as pe
import pytest

from egret.data.model_data import ModelData
from egret.models.acopf import (create_psv_acopf_model,
                                create_rsv_acopf_model,
                                create_atan_acopf_model,
                                create_riv_acopf_model)

current_dir = os.path.dirname(os.path.abspath(__file__))
test_case = os.path.join(current_dir, 'transmission_test_instances',
                         'test_instances', 'hvdc_test_case3.json')

POLYNOMIAL_P_COST = {'data_type': 'cost_curve', 'cost_curve_type': 'polynomial',
                     'values': {'0': 500., '1': 40., '2': 0.05}}
PIECEWISE_P_COST = {'data_type': 'cost_curve', 'cost_curve_type': 'piecewise',
                    'values': [(-200., 0.), (0., 1000.), (200., 20000.)]}
PIECEWISE_Q_COST = {'data_type': 'cost_curve', 'cost_curve_type': 'piecewise',
                    'values': [(-100., 500.), (0., 0.), (100., 600.)]}


def _make_model_data(p_cost, q_cost):
    md = ModelData.read(test_case)
    for gen in md.data['elements']['generator'].values():
        gen['p_cost'] = copy.deepcopy(p_cost)
        if q_cost is None:
            gen.pop('q_cost', None)
        else:
            gen['q_cost'] = copy.deepcopy(q_cost)
    return md


@pytest.mark.parametrize('acopf_model', [
    create_psv_acopf_model, create_rsv_acopf_model,
    create_atan_acopf_model, create_riv_acopf_model,
])
def test_piecewise_q_cost_with_polynomial_p_cost(acopf_model):
    '''
    A piecewise reactive power cost must not drag in the real power step
    variables, which do not exist for a polynomial real power cost.
    '''
    model, _md = acopf_model(_make_model_data(POLYNOMIAL_P_COST, PIECEWISE_Q_COST))
    assert hasattr(model, 'qg_operating_cost')


@pytest.mark.parametrize('acopf_model', [
    create_psv_acopf_model, create_riv_acopf_model,
])
def test_piecewise_q_cost_with_piecewise_p_cost(acopf_model):
    '''Both curves piecewise: the two step variable sets coexist.'''
    model, _md = acopf_model(_make_model_data(PIECEWISE_P_COST, PIECEWISE_Q_COST))
    assert hasattr(model, 'qg_operating_cost')


def test_no_q_cost_is_unaffected():
    '''Without q_cost the real power path is unchanged.'''
    model, _md = create_psv_acopf_model(_make_model_data(PIECEWISE_P_COST, None))
    assert hasattr(model, 'pg_operating_cost')
