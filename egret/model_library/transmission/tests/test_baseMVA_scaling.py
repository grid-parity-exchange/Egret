#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

import pytest
import math

from egret.data.model_data import ModelData
from egret.model_library.transmission.tx_utils import \
        scale_ModelData_to_pu, unscale_ModelData_to_pu
from egret.models.unit_commitment import solve_unit_commitment
from egret.models.tests.test_unit_commitment import test_solver

def test_scale_unscale():
    md = ModelData.read(
            '../../../models/tests/uc_test_instances/'
            'test_scuc_full_enforce_relaxed_sol.json')

    ## do type conversions
    original_base_MVA = md.data['system']['baseMVA']
    md.data['system']['baseMVA'] = 1.

    scale_ModelData_to_pu(md, inplace=True)

    md.data['system']['baseMVA'] = original_base_MVA

    md_transformed = scale_ModelData_to_pu(md, inplace=False)

    unscale_ModelData_to_pu(md_transformed, inplace=True)

    assert md.data['system'] == md_transformed.data['system']
    for esn, esd in md.data['elements'].items():
        for en, ed in esd.items():
            assert ed == md_transformed.data['elements'][esn][en]

    for esn, esd in md_transformed.data['elements'].items():
        for en, ed in esd.items():
            assert ed == md.data['elements'][esn][en]

def test_scaling_spot_check():
    md = ModelData.read(
            '../../../models/tests/uc_test_instances/'
            'test_scuc_full_enforce_relaxed_sol.json')

    baseMVA = md.data['system']['baseMVA']

    md_scaled = scale_ModelData_to_pu(md, inplace=False)

    md_scaled_unscaled = unscale_ModelData_to_pu(md_scaled, inplace=False)

    ## commitment should be unchanged
    assert md.data['elements']['generator']['101_STEAM_3_t']['commitment']['values'][10] == \
        md_scaled.data['elements']['generator']['101_STEAM_3_t']['commitment']['values'][10] == \
        md_scaled_unscaled.data['elements']['generator']['101_STEAM_3_t']['commitment']['values'][10]

    ## as should production cost
    assert md.data['elements']['generator']['101_STEAM_3_t']['production_cost']['values'][10] == \
        md_scaled.data['elements']['generator']['101_STEAM_3_t']['production_cost']['values'][10] == \
        md_scaled_unscaled.data['elements']['generator']['101_STEAM_3_t']['production_cost']['values'][10]

    ## as should voltage angle
    assert md.data['elements']['bus']['Alber']['va']['values'][10] == \
        md_scaled.data['elements']['bus']['Alber']['va']['values'][10] == \
        md_scaled_unscaled.data['elements']['bus']['Alber']['va']['values'][10]

    ## pg should be scaled
    assert md.data['elements']['generator']['101_STEAM_3_t']['pg']['values'][10] == \
        md_scaled.data['elements']['generator']['101_STEAM_3_t']['pg']['values'][10]/baseMVA == \
        md_scaled_unscaled.data['elements']['generator']['101_STEAM_3_t']['pg']['values'][10]

    ## load should be scaled
    assert md.data['elements']['bus']['Alber']['pl']['values'][10] == \
        md_scaled.data['elements']['bus']['Alber']['pl']['values'][10]/baseMVA == \
        md_scaled_unscaled.data['elements']['bus']['Alber']['pl']['values'][10]

    ## load should be scaled
    assert md.data['elements']['load']['Alber']['p_load']['values'][10] == \
        md_scaled.data['elements']['load']['Alber']['p_load']['values'][10]/baseMVA == \
        md_scaled_unscaled.data['elements']['load']['Alber']['p_load']['values'][10]

    ## flows should be scaled
    assert md.data['elements']['branch']['A22']['pf']['values'][20] == \
        md_scaled.data['elements']['branch']['A22']['pf']['values'][20]/baseMVA == \
        md_scaled_unscaled.data['elements']['branch']['A22']['pf']['values'][20]

    ## contingency flows should also be scaled
    assert md.data['elements']['contingency']['A1']['monitored_branches']['values'][10]['A11']['pf'] == \
        md_scaled.data['elements']['contingency']['A1']['monitored_branches']['values'][10]['A11']['pf']/baseMVA == \
        md_scaled_unscaled.data['elements']['contingency']['A1']['monitored_branches']['values'][10]['A11']['pf']

    ## lmp should be inversly scaled
    assert md.data['elements']['bus']['Alber']['lmp']['values'][10] == \
        md_scaled.data['elements']['bus']['Alber']['lmp']['values'][10]*baseMVA == \
        md_scaled_unscaled.data['elements']['bus']['Alber']['lmp']['values'][10]

    ## reserve prices should be inversly scaled
    assert md.data['system']['reserve_price']['values'][18] == \
        md_scaled.data['system']['reserve_price']['values'][18]*baseMVA == \
        md_scaled_unscaled.data['system']['reserve_price']['values'][18]

    ## shortfall price should be inversly scaled
    assert md.data['system']['reserve_shortfall_cost'] == \
        md_scaled.data['system']['reserve_shortfall_cost']*baseMVA == \
        md_scaled_unscaled.data['system']['reserve_shortfall_cost']

def test_scaling_solve():
    md = ModelData.read(
            '../../../models/tests/uc_test_instances/tiny_uc_7.json')

    assert md.data['system']['baseMVA'] == 1.
    mdo_unscaled = solve_unit_commitment(md, test_solver, relaxed=True)

    md.data['system']['baseMVA'] = 100.
    mdo_scaled = solve_unit_commitment(md, test_solver, relaxed=True)

    assert math.isclose(mdo_scaled.data['system']['total_cost'], mdo_unscaled.data['system']['total_cost'])
