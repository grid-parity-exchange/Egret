#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
Interfaces that are slack-relaxed under SlackType.TRANSMISSION_LIMITS.

`load_params` builds the `InterfacesWithSlack` set from the interfaces that
either carry a `violation_penalty` or, under SlackType.TRANSMISSION_LIMITS,
have none.  The second branch used to append the stray variable `bn` left
over from the preceding branch loop instead of the current interface name,
so building such a model raised

    ValueError: Cannot add value <last branch name> to Set InterfacesWithSlack

Neither the test instances nor the test suite exercise interfaces together
with SlackType.TRANSMISSION_LIMITS, which is why this went unnoticed.
'''
import pyomo.environ as pe

from egret.data.model_data import ModelData
from egret.models.unit_commitment import create_tight_unit_commitment_model
from egret.model_library.unit_commitment.uc_utils import SlackType

INTERFACE_NAME = 'TIE'


def _make_model_data():
    return {
        'system': {
            'baseMVA': 100.,
            'reference_bus': 'Bus1',
            'reference_bus_angle': 0.,
            'time_keys': ['1', '2'],
            'time_period_length_minutes': 60,
            'load_mismatch_cost': 1000.,
            'reserve_shortfall_cost': 100.,
        },
        'elements': {
            'bus': {
                'Bus1': {'base_kv': 230},
                'Bus2': {'base_kv': 230},
            },
            'load': {
                'Load2': {
                    'bus': 'Bus2',
                    'in_service': True,
                    'p_load': {'data_type': 'time_series', 'values': [100., 100.]},
                },
            },
            'branch': {
                'Branch1': {
                    'from_bus': 'Bus1',
                    'to_bus': 'Bus2',
                    'resistance': 0.01,
                    'reactance': 0.05,
                    'rating_long_term': 300.,
                    'branch_type': 'line',
                    'angle_diff_min': -180.,
                    'angle_diff_max': 180.,
                    'in_service': True,
                },
            },
            'interface': {
                INTERFACE_NAME: {
                    'lines': ['Branch1'],
                    'line_orientation': [1],
                    'minimum_limit': -300.,
                    'maximum_limit': 300.,
                },
            },
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
                    'ramp_up_60min': 300.,
                    'ramp_down_60min': 300.,
                    'min_up_time': 1,
                    'min_down_time': 1,
                    'initial_status': 1,
                    'initial_p_output': 100.,
                    'p_cost': {'data_type': 'cost_curve',
                               'cost_curve_type': 'polynomial',
                               'values': {'0': 0., '1': 30., '2': 0.}},
                },
            },
            'storage': {},
        },
    }


def test_interface_without_penalty_is_slacked_with_transmission_limits():
    '''
    An interface that carries no violation_penalty must be added to
    InterfacesWithSlack when slack_type is TRANSMISSION_LIMITS.
    '''
    model = create_tight_unit_commitment_model(
        ModelData(_make_model_data()),
        network_constraints='ptdf_power_flow',
        slack_type=SlackType.TRANSMISSION_LIMITS)

    slacked = [str(i) for i in model.InterfacesWithSlack]
    assert INTERFACE_NAME in slacked, \
        'expected {} in InterfacesWithSlack, got {}'.format(INTERFACE_NAME, slacked)


def test_interface_without_penalty_is_not_slacked_with_bus_balance():
    '''Under the default BUS_BALANCE slacks, the interface set stays empty.'''
    model = create_tight_unit_commitment_model(
        ModelData(_make_model_data()),
        network_constraints='ptdf_power_flow',
        slack_type=SlackType.BUS_BALANCE)

    assert [str(i) for i in model.InterfacesWithSlack] == []


def test_interface_with_penalty_is_always_slacked():
    '''An interface carrying a violation_penalty is slacked under any slack type.'''
    md_dict = _make_model_data()
    md_dict['elements']['interface'][INTERFACE_NAME]['violation_penalty'] = 1000.

    model = create_tight_unit_commitment_model(
        ModelData(md_dict),
        network_constraints='ptdf_power_flow',
        slack_type=SlackType.BUS_BALANCE)

    assert INTERFACE_NAME in [str(i) for i in model.InterfacesWithSlack]
