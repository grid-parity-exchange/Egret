#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
Storage efficiency key names.

`docs/types_sheet.xlsx` documents the discharge efficiency attribute as
`discharge_efficiency`, and `egret.parsers.prescient_dat_parser` emits that
spelling, but `params.load_params` used to read a misspelled
`discharge_efficienty`.  Data written with the documented spelling was
therefore silently ignored and the discharge efficiency fell back to its
default of 1.0.  Both spellings must be honoured.
'''
import pyomo.environ as pe
import pytest

from egret.data.model_data import ModelData
from egret.models.unit_commitment import create_tight_unit_commitment_model

STORAGE_NAME = 'BATTERY'


def _make_model_data():
    md_dict = {
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
            },
            'load': {
                'Load1': {
                    'bus': 'Bus1',
                    'in_service': True,
                    'p_load': {'data_type': 'time_series', 'values': [100., 100.]},
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
                    'p_max': 200.,
                    'ramp_up_60min': 200.,
                    'ramp_down_60min': 200.,
                    'min_up_time': 1,
                    'min_down_time': 1,
                    'initial_status': 1,
                    'initial_p_output': 100.,
                    'p_cost': {'data_type': 'cost_curve',
                               'cost_curve_type': 'polynomial',
                               'values': {'0': 0., '1': 30., '2': 0.}},
                },
            },
            'storage': {
                STORAGE_NAME: {
                    'bus': 'Bus1',
                    'in_service': True,
                    'energy_capacity': 400.,
                    'initial_state_of_charge': 0.5,
                    'end_stage_of_charge': 0.5,
                    'minimum_state_of_charge': 0.,
                    'charge_efficiency': 0.9,
                    'max_charge_rate': 100.,
                    'min_charge_rate': 0.,
                    'max_discharge_rate': 100.,
                    'min_discharge_rate': 0.,
                    'ramp_up_input_60min': 100.,
                    'ramp_down_input_60min': 100.,
                    'ramp_up_output_60min': 100.,
                    'ramp_down_output_60min': 100.,
                    'retention_rate_60min': 1.,
                },
            },
        },
    }
    return md_dict


def _discharge_efficiency_in_model(eff_key):
    """
    Build a UC model with the discharge efficiency given under `eff_key`
    and return the resulting OutputEfficiencyEnergy value.
    """
    md_dict = _make_model_data()
    md_dict['elements']['storage'][STORAGE_NAME][eff_key] = 0.75
    model = create_tight_unit_commitment_model(
        ModelData(md_dict), network_constraints='copperplate_power_flow')
    return pe.value(model.OutputEfficiencyEnergy[STORAGE_NAME])


def test_documented_discharge_efficiency_is_used():
    '''The spelling documented in docs/types_sheet.xlsx must take effect.'''
    assert _discharge_efficiency_in_model('discharge_efficiency') == pytest.approx(0.75)


def test_legacy_misspelled_discharge_efficiency_still_works():
    '''The historical misspelling must keep working for existing datasets.'''
    assert _discharge_efficiency_in_model('discharge_efficienty') == pytest.approx(0.75)


def test_discharge_efficiency_defaults_to_one_when_absent():
    '''Without either key the documented default of 1.0 still applies.'''
    md_dict = _make_model_data()
    model = create_tight_unit_commitment_model(
        ModelData(md_dict), network_constraints='copperplate_power_flow')
    assert pe.value(model.OutputEfficiencyEnergy[STORAGE_NAME]) == pytest.approx(1.0)
