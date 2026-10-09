#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
Reactive power balance with no slack.

`power_balance._add_egret_power_flow` dispatches on `slack_type` for the
reactive power mismatch, exactly as it does for the real power mismatch.
That branch had two defects and could not run at all:

* the `SlackType.NONE` case called a misspelled
  `_add_blank_q_load_mistmatch` (the helper is `_add_blank_q_load_mismatch`),
  raising `NameError`;
* `_add_reactive_limits(model)` was called without the `grid` argument its
  signature required, raising `TypeError`.

All built-in unit commitment network formulations hard-code
`reactive_power=False`, so the branch is currently unreachable through the
public API.  These tests exercise it directly, through a power balance
formulation registered for the purpose, so that it keeps working now that it
is repaired.
'''
import pyomo.environ as pe
import pytest

from egret.data.model_data import ModelData
from egret.model_library.unit_commitment import power_balance as lb_power_balance
from egret.model_library.unit_commitment.uc_model_generator import (UCFormulation,
                                                                    generate_model)
from egret.model_library.unit_commitment.uc_utils import SlackType, add_model_attr

# The combination that is unreachable through the public API: reactive power
# enabled, with no slack anywhere.
_REACTIVE_NO_SLACK_NAME = 'reactive_no_slack_power_flow'


@add_model_attr('power_balance', requires={'data_loader': None,
                                           'power_vars': None,
                                           'non_dispatchable_vars': None,
                                           'storage_service': None})
def _reactive_no_slack_power_flow(model, slack_type=SlackType.NONE):
    lb_power_balance._add_egret_power_flow(
        model,
        lb_power_balance._copperplate_approx_network_model,
        reactive_power=True,
        slack_type=SlackType.NONE)


# generate_model looks the power balance formulation up by name on the module
setattr(lb_power_balance, _REACTIVE_NO_SLACK_NAME, _reactive_no_slack_power_flow)


def _make_model_data():
    return {
        'system': {
            'baseMVA': 100.,
            'reference_bus': 'Bus1',
            'reference_bus_angle': 0.,
            'time_keys': ['1', '2'],
            'time_period_length_minutes': 60,
            'load_mismatch_cost': 1000.,
            'q_load_mismatch_cost': 1000.,
            'reserve_shortfall_cost': 100.,
        },
        'elements': {
            'bus': {'Bus1': {'base_kv': 230}},
            'load': {
                'Load1': {
                    'bus': 'Bus1',
                    'in_service': True,
                    'p_load': {'data_type': 'time_series', 'values': [100., 100.]},
                    'q_load': {'data_type': 'time_series', 'values': [30., 30.]},
                },
            },
            'branch': {},
            'interface': {},
            'zone': {},
            'area': {},
            'generator': {
                'GEN1': {
                    'generator_type': 'thermal',
                    'bus': 'Bus1',
                    'fuel': 'gas',
                    'in_service': True,
                    'p_min': 10.,
                    'p_max': 300.,
                    'q_min': -150.,
                    'q_max': 150.,
                    'ramp_up_60min': 300.,
                    'ramp_down_60min': 300.,
                    'min_up_time': 1,
                    'min_down_time': 1,
                    'initial_status': 1,
                    'initial_p_output': 100.,
                    'initial_q_output': 30.,
                    'p_cost': {'data_type': 'cost_curve',
                               'cost_curve_type': 'polynomial',
                               'values': {'0': 0., '1': 30., '2': 0.}},
                },
            },
            'storage': {},
        },
    }


def _formulation():
    return UCFormulation(status_vars='garver_3bin_vars',
                         power_vars='garver_power_vars',
                         reserve_vars='garver_power_avail_vars',
                         generation_limits='gentile_generation_limits',
                         ramping_limits='damcikurt_ramping',
                         production_costs='KOW_production_costs',
                         uptime_downtime='rajan_takriti_UT_DT',
                         startup_costs='KOW_startup_costs',
                         network_constraints=_REACTIVE_NO_SLACK_NAME)


def test_reactive_power_with_no_slack_builds():
    '''The SlackType.NONE reactive power branch must build without error.'''
    model = generate_model(ModelData(_make_model_data()), _formulation())
    assert model.reactive_power
    # the blank reactive mismatch parameters come from _add_blank_q_load_mismatch
    assert hasattr(model, 'LoadMismatchCostReactive')
    assert hasattr(model, 'ReactivePowerGenerated')
    assert hasattr(model, 'EnforceReactiveUpperLimit')
    assert hasattr(model, 'EnforceReactiveLowerLimit')


def test_blank_q_load_mismatch_helper_is_callable():
    '''Guard against the helper being renamed without the call site following.'''
    assert callable(lb_power_balance._add_blank_q_load_mismatch)


def test_reactive_limits_helper_accepts_model_only():
    '''_add_reactive_limits is called with the model alone; that must be legal.'''
    from egret.model_library.unit_commitment import generation_limits
    model = ModelData(_make_model_data())
    uc = generate_model(model, _formulation())
    # calling it again with just the model must not raise TypeError
    generation_limits._add_reactive_limits(uc)
