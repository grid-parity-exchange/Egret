#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

from pyomo.environ import *
import math

from .uc_utils import add_model_attr, build_uc_time_mapping
from .status_vars import _is_relaxed

@add_model_attr('dual_fuel', requires = {'data_loader': None,
                                         'status_vars': ['garver_3bin_vars','garver_2bin_vars','garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                         'power_vars': None,
                                         })
def dual_fuel_model(model):
    '''
    Defines extra variables and constraints for tracking fuel consumption for a 
    dual-fuel generator
    '''

    md = model.model_data

    system = md.data['system']
    time_keys = system['time_indices']
    TimeMapper = build_uc_time_mapping(time_keys)
    
    dual_fuel_attrs = md.attributes(element_type='generator', generator_type='thermal', aux_fuel_capable=True)

    relaxed = _is_relaxed(model)
    
    ## Binary or Unit Variables for unit status: on, start, stop for both fuel types

    def binary_or_unit_interval(m, g, t):
        if relaxed or \
                ('aux_fuel_blending' in dual_fuel_attrs \
                        and g in dual_fuel_attrs['aux_fuel_blending'] \
                        and dual_fuel_attrs['aux_fuel_blending'][g] ):
            return UnitInterval
        return Binary

    model.UnitOnPriFuel = Var(model.DualFuelGenerators, model.TimePeriods, within=binary_or_unit_interval)
    model.UnitOnAuxFuel = Var(model.DualFuelGenerators, model.TimePeriods, within=binary_or_unit_interval)

    model.UnitStartPriFuel = Var(model.DualFuelGenerators, model.TimePeriods, within=binary_or_unit_interval)
    model.UnitStartAuxFuel = Var(model.DualFuelGenerators, model.TimePeriods, within=binary_or_unit_interval)

    model.UnitStopPriFuel = Var(model.DualFuelGenerators, model.TimePeriods, within=binary_or_unit_interval)
    model.UnitStopAuxFuel = Var(model.DualFuelGenerators, model.TimePeriods, within=binary_or_unit_interval)

    ## Linking constraints for status to normal status variables

    def linking_unit_on(m, g, t):
        return m.UnitOnPriFuel[g,t] + m.UnitOnAuxFuel[g,t] == m.UnitOn[g,t]
    model.UnitOnLink = Constraint(model.DualFuelGenerators, model.TimePeriods, rule=linking_unit_on)

    def linking_unit_start(m, g, t):
        return m.UnitStartPriFuel[g,t] + m.UnitStartAuxFuel[g,t] == m.UnitStart[g,t]
    model.UnitStartLink = Constraint(model.DualFuelGenerators, model.TimePeriods, rule=linking_unit_start)

    def linking_unit_stop(m, g, t):
        return m.UnitStopPriFuel[g,t] + m.UnitStopAuxFuel[g,t] == m.UnitStop[g,t]
    model.UnitStopLink = Constraint(model.DualFuelGenerators, model.TimePeriods, rule=linking_unit_stop)

    model.UnitSwitchOperating = Param(model.DualFuelGenerators, within=Boolean, default=False, initialize=dual_fuel_attrs.get('aux_fuel_online_switching')) 


    ## These enforce that if the generator starts with a certain fuel type it
    ## must maintain it until shutdown.
    def pri_logical_constr(m, g, t):
        if value(m.UnitSwitchOperating[g]):
            return Constraint.Skip
        elif t==value(m.InitialTime):
            if value(m.UnitOnT0[g]):
                pri_fuel_init = 1 - dual_fuel_attrs['aux_fuel_supply_initial'][g]
                return m.UnitOnPriFuel[g,t] - pri_fuel_init == m.UnitStartPriFuel[g,t] - m.UnitStopPriFuel[g,t]
            else:
                return Constraint.Skip
        else:
            return m.UnitOnPriFuel[g,t] - m.UnitOnPriFuel[g,t-1] == m.UnitStartPriFuel[g,t] - m.UnitStopPriFuel[g,t]
    model.PriLogicalConst = Constraint(model.DualFuelGenerators, model.TimePeriods, rule=pri_logical_constr)

    def aux_logical_constr(m, g, t):
        if value(m.UnitSwitchOperating[g]):
            return Constraint.Skip
        elif t==value(m.InitialTime):
            if value(m.UnitOnT0[g]):
                aux_fuel_init = dual_fuel_attrs['aux_fuel_supply_initial'][g]
                return m.UnitOnAuxFuel[g,t] - aux_fuel_init == m.UnitStartAuxFuel[g,t] - m.UnitStopAuxFuel[g,t]
            else:
                return Constraint.Skip
        else:
            return m.UnitOnAuxFuel[g,t] - m.UnitOnAuxFuel[g,t-1] == m.UnitStartAuxFuel[g,t] - m.UnitStopAuxFuel[g,t]
    model.AuxLogicalConst = Constraint(model.DualFuelGenerators, model.TimePeriods, rule=pri_logical_constr)

    ## Producion Costs


