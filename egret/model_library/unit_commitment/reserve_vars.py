#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## system variables and constraints
from pyomo.environ import *
import math

from .uc_utils import add_model_attr, get_linear_expr
component_name = 'reserve_vars'

def check_reserve_requirement(model):
    system = model.model_data.data['system']
    return ('reserve_requirement' in system)

def _add_zero_reserve_hooks(model):

    def max_power_available_above_min_rule(m, g, t):
        return m.PowerGeneratedAboveMinimum[g,t]
    model.MaximumPowerAvailableAboveMinimum = Expression(model.ThermalGenerators, model.TimePeriods, rule=max_power_available_above_min_rule)
    model._get_maximum_power_available_above_minimum_lists = model._get_power_generated_above_minimum_lists

    def max_power_available_rule(m,g,t):
        return m.PowerGenerated[g,t]
    model.MaximumPowerAvailable = Expression(model.ThermalGenerators, model.TimePeriods, rule=max_power_available_rule)
    model._get_maximum_power_available_lists = model._get_power_generated_lists

    model.ReserveProvided = Param(model.ThermalGenerators, model.TimePeriods, default=0.0)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'power_vars': None,
                                            })
def garver_power_avail_vars(model):

    '''
    These never appear in Garver's paper, but they are an adaption of the 
    idea from the Carrion-Arroyo paper for maximum power available 
    to consider maximum power available over minimum
    '''

    ## only add reserves if the user specified them
    if not check_reserve_requirement(model):
        _add_zero_reserve_hooks(model)
        return

    # amount of power produced by each generator above minimum, at each time period.
    def garver_power_bounds_rule(m, g, t):
        return (0, m.MaximumPowerOutput[g,t]-m.MinimumPowerOutput[g,t])

    # maximum power output above minimum for each generator, at each time period.
    model.MaximumPowerAvailableAboveMinimum = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=garver_power_bounds_rule)
    model._get_maximum_power_available_above_minimum_lists = lambda m,g,t : ([m.MaximumPowerAvailableAboveMinimum[g,t]], [1.])
    
    ## Note: thes only get used in system balance constraints
    linear_expr = get_linear_expr(model.UnitOn)
    model._get_maximum_power_available_lists = lambda m,g,t : ([m.MaximumPowerAvailableAboveMinimum[g,t], m.UnitOn[g,t]], [1., m.MinimumPowerOutput[g,t]])
    def maximum_power_avaiable_expr_rule(m, g, t):
        return linear_expr(*m._get_maximum_power_available_lists(m,g,t))
    model.MaximumPowerAvailable = Expression(model.ThermalGenerators, model.TimePeriods, rule=maximum_power_avaiable_expr_rule)

    # m.MinimumPowerOutput[g] * m.UnitOn[g, t] <= m.PowerGenerated[g,t] <= m.MaximumPowerAvailable[g, t] <= m.MaximumPowerOutput[g] * m.UnitOn[g, t]
    # BK -- first <= now handled by bounds
    
    linear_expr = get_linear_expr(model.PowerGeneratedAboveMinimum)
    def enforce_generator_output_limits_rule_part_b(m, g, t):
        return (None, linear_expr(
                linear_vars=[m.PowerGeneratedAboveMinimum[g,t], m.MaximumPowerAvailableAboveMinimum[g,t]],
                linear_coefs=[1,-1]
                                        ), 0)
    model.EnforceGeneratorOutputLimitsPartB = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_generator_output_limits_rule_part_b)

    ## BK -- for reserve pricing
    _get_reserves_provided_lists = lambda m,g,t : ([m.MaximumPowerAvailableAboveMinimum[g,t], m.PowerGeneratedAboveMinimum[g,t]], [1., -1.])
    def reserve_provided_expr_rule(m, g, t):
        #return m.MaximumPowerAvailableAboveMinimum[g,t] - m.PowerGeneratedAboveMinimum[g,t]
        return linear_expr( *_get_reserves_provided_lists(m,g,t) )
    model.ReserveProvided = Expression(model.ThermalGenerators, model.TimePeriods, rule=reserve_provided_expr_rule) 

    return

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'power_vars': None,
                                            })
def MLR_reserve_vars(model):
    '''
    Reserves provided variables as in

    G. Morales-Espana, J. M. Latorre, and A. Ramos. Tight and compact MILP
    formulation for the thermal unit commitment problem. IEEE Transactions on
    Power Systems, 28(4):4897–4908, 2013.

    '''

    ## only add reserves if the user specified them
    if not check_reserve_requirement(model):
        _add_zero_reserve_hooks(model)
        return

    # amount of power produced by each generator above minimum, at each time period.
    def garver_power_bounds_rule(m, g, t):
        return (0, m.MaximumPowerOutput[g,t]-m.MinimumPowerOutput[g,t])

    # variable for reserves offered
    model.ReserveProvided = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=garver_power_bounds_rule)

    linear_expr = get_linear_expr(model.UnitOn)

    ## Note: thes only get used in system balance constraints
    def _get_maximum_power_available_lists(m,g,t):
        linear_vars, linear_coefs = m._get_power_generated_lists(m,g,t)
        linear_vars.append(m.ReserveProvided[g,t])
        linear_coefs.append(1.)
        return linear_vars, linear_coefs

    def maximum_power_avaiable_expr_rule(m, g, t):
        return linear_expr(*_get_maximum_power_available_lists(m,g,t))
    model.MaximumPowerAvailable = Expression(model.ThermalGenerators, model.TimePeriods, rule=maximum_power_avaiable_expr_rule)

    model._get_maximum_power_available_lists = _get_maximum_power_available_lists

    ## Note: thes only get used in system balance constraints
    def _get_maximum_power_available_above_minimum_lists(m,g,t):
        linear_vars, linear_coefs = m._get_power_generated_above_minimum_lists(m,g,t)
        linear_vars.append(m.ReserveProvided[g,t])
        linear_coefs.append(1.)
        return linear_vars, linear_coefs

    ## make this an expression so it propogates nicely with the term above
    def maximum_power_available_above_minimum_expr_rule(m, g, t):
        return linear_expr(*_get_maximum_power_available_above_minimum_lists(m,g,t))
    model.MaximumPowerAvailableAboveMinimum = Expression(model.ThermalGenerators, model.TimePeriods, rule=maximum_power_available_above_minimum_expr_rule)

    model._get_maximum_power_available_above_minimum_lists = _get_maximum_power_available_above_minimum_lists

    return


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'power_vars': None,
                                            })
def CA_power_avail_vars(model):
    '''
    MaximumPowerAvailable var in [0, MaximumPowerOutput], as in

    Carrion, M. and Arroyo, J. (2006) A Computationally Efficient Mixed-Integer
    Liner Formulation for the Thermal Unit Commitment Problem. IEEE Transactions
    on Power Systems, Vol. 21, No. 3, Aug 2006.
    '''

    ## only add reserves if the user specified them
    if not check_reserve_requirement(model):
        _add_zero_reserve_hooks(model)
        return

    # amount of power produced by each generator, at each time period.
    def power_bounds_rule(m, g, t):
        return (0, m.MaximumPowerOutput[g,t])
    model.MaximumPowerAvailable = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=power_bounds_rule) 
    model._get_maximum_power_available_lists = lambda m,g,t : ([m.MaximumPowerAvailable[g,t]], [1.])

    # this is useful for the damci_kurt ramping inequality
    # I think this seemlessly handles many cases when this expression is on the LHS,
    # including a variant of the state transition formulation
    linear_expr = get_linear_expr(model.UnitOn)

    model._get_maximum_power_available_above_minimum_lists = lambda m,g,t : ([m.MaximumPowerAvailable[g,t], m.UnitOn[g,t]], [1., -m.MinimumPowerOutput[g,t]])
    def maximum_power_available_above_minimum_expr_rule(m, g, t):
        return linear_expr(*m._get_maximum_power_available_above_minimum_lists(m,g,t))
    model.MaximumPowerAvailableAboveMinimum = Expression(model.ThermalGenerators, model.TimePeriods, rule=maximum_power_available_above_minimum_expr_rule)

    # m.MinimumPowerOutput[g] * m.UnitOn[g, t] <= m.PowerGenerated[g,t] <= m.MaximumPowerAvailable[g, t] <= m.MaximumPowerOutput[g] * m.UnitOn[g, t]
    # BK -- first <= now handled by bounds
    def enforce_generator_output_limits_rule_part_b(m, g, t):
       return m.PowerGenerated[g,t] <= m.MaximumPowerAvailable[g, t]

    model.EnforceGeneratorOutputLimitsPartB = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_generator_output_limits_rule_part_b)

    ## BK -- for reserve pricing
    def reserve_provided_expr_rule(m, g, t):
        return m.MaximumPowerAvailable[g,t] - m.PowerGenerated[g,t]
    model.ReserveProvided = Expression(model.ThermalGenerators, model.TimePeriods, rule=reserve_provided_expr_rule) 
