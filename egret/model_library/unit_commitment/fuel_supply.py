#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## file for all the ancillary services
from pyomo.environ import *
import math

from .uc_utils import add_model_attr, uc_time_helper 
from .status_vars import _is_relaxed

## For now, this function will depend on having one of the production costs with
## PiecewiseProduction variables, and start-up costs which use auxillary variables
@add_model_attr('fuel_supply', requires = {'data_loader': None,
                                           'fuel_consumption':['fuel_consumption_model'],
                                           })
def fuel_supply_model(model):
    '''
    Defines a fuel supply model for generators. For now, just considers
    'instantaneous' supply, i.e., like NG
    '''
    md = model.model_data

    TimeMapper = uc_time_helper(model.TimePeriods)

    ## instantaneous fuel supply model
    inst_fuel_supply_attrs = md.attributes(element_type='fuel_supply', fuel_supply_type='instantaneous')
    thermal_gen_attrs = md.attributes(element_type='generator', generator_type='thermal')

    model.InstantaneousFuelSupplies = Set(initialize=inst_fuel_supply_attrs['names'])

    def gens_using_inst_fuel_supply_init(m, f):
        for g in m.FuelSupplyGenerators:
            if thermal_gen_attrs['fuel_supply'][g] == f:
                yield g
    model.ThermalGeneratorsUsingInstFuelSupply = Set(model.InstantaneousFuelSupplies, initialize=gens_using_inst_fuel_supply_init,)

    def gens_using_inst_fuel_supply_for_aux_init(m, f):
        for g in m.DualFuelGenerators:
            if 'aux_fuel_supply' in thermal_gen_attrs:
                if g in thermal_gen_attrs['aux_fuel_supply'] and thermal_gen_attrs['aux_fuel_supply'][g] == f:
                    yield g
    model.DualFuelGeneratorsUsingInstFuelSupplyForAuxFuel = Set(model.InstantaneousFuelSupplies, initialize=gens_using_inst_fuel_supply_for_aux_init,)

    model.InstFuelSupply = Param(model.InstantaneousFuelSupplies, model.TimePeriods, initialize=TimeMapper(inst_fuel_supply_attrs['fuel_available']))

    def total_fuel_consumed_expr(m, f, t):
        return sum(m.PrimaryFuelConsumed[g,t] for g in m.ThermalGeneratorsUsingInstFuelSupply[f]) + sum(m.AuxiliaryFuelConsumed[g,t] for g in m.DualFuelGeneratorsUsingInstFuelSupplyForAuxFuel[f])
    model.TotalFuelConsumedAtInstFuelSupply = Expression(model.InstantaneousFuelSupplies, model.TimePeriods, rule=total_fuel_consumed_expr)

    def total_fuel_consumed_rule(m, f, t):
        if m.ThermalGeneratorsUsingInstFuelSupply[f] or m.DualFuelGeneratorsUsingInstFuelSupplyForAuxFuel[f]:
            return m.TotalFuelConsumedAtInstFuelSupply[f,t] <= m.InstFuelSupply[f,t]
        else:
            if t == m.TimePeriods.first():
                print('WARNING: no generators attached to fuel_supply {}'.format(f))
            return Constraint.Feasible
    model.FuelLimitConstr = Constraint(model.InstantaneousFuelSupplies, model.TimePeriods, rule=total_fuel_consumed_rule)
    ## end instantaneous fuel supply model
