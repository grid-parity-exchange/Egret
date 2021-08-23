#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## file for power and reserve variables
from pyomo.environ import *
import math

from .uc_utils import add_model_attr, get_linear_expr
component_name = 'power_vars'

def _add_reactive_power_vars(model):

    def reactive_power_bounds_rule(m,g,t):
        return (m.MinimumReactivePowerOutput[g,t], m.MaximumReactivePowerOutput[g,t])
    model.ReactivePowerGenerated = Var(model.ThermalGenerators, model.TimePeriods, within=Reals, bounds=reactive_power_bounds_rule)

def _add_power_generated_startup_shutdown(model):
    assert model.InitialTime == 1

    # first, discover if we have startup/shutdown
    # curves in the model
    model_has_startup_shutdown_curves = False
    for s in model.StartupCurve.values():
        if len(s) > 0:
            model_has_startup_shutdown_curves = True
            break
    if not model_has_startup_shutdown_curves:
        for s in model.ShutdownCurve.values():
            if len(s) > 0:
                model_has_startup_shutdown_curves = True
                break

    if model_has_startup_shutdown_curves:
        # check the status vars to see if we're compatible
        # with startup/shutdown curves
        if model.status_vars not in ['garver_2bin_vars', 'garver_3bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars']:
            raise RuntimeError(f"Status variable formulation {model.status_vars} is not compatible with startup or shutdown curves")

        def power_generated_startup_shutdown_expr_rule(m, g, t):
            startup_curve = m.StartupCurve[g]
            shutdown_curve = m.ShutdownCurve[g]
            time_periods_before_startup = value(m.TimePeriodsBeforeStartup[g])
            time_periods_since_shutdown = value(m.TimePeriodsSinceShutdown[g])

            future_startup_past_shutdown_production = 0.
            future_startup_power_index = time_periods_before_startup + m.NumTimePeriods - t
            if future_startup_power_index <= len(startup_curve):
                future_startup_past_shutdown_production += startup_curve.at(future_startup_power_index)

            past_shutdown_power_index = time_periods_since_shutdown + t
            if past_shutdown_power_index <= len(shutdown_curve):
                future_startup_past_shutdown_production += shutdown_curve.at(past_shutdown_power_index)

            linear_vars, linear_coefs = m._get_power_generated_lists(m,g,t)
            for startup_idx in range(1, min( len(startup_curve)+1, m.NumTimePeriods+1-t )):
                linear_vars.append(m.UnitStart[g,t+startup_idx])
                linear_coefs.append(startup_curve.at(startup_idx))
            for shutdown_idx in range(1, min( len(shutdown_curve)+1, t+1 )):
                linear_vars.append(m.UnitStop[g,t-shutdown_idx+1])
                linear_coefs.append(shutdown_curve.at(shutdown_idx))
            linear_expr = get_linear_expr(m.UnitOn, m.UnitStart, m.UnitStop)
            return linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs, constant=future_startup_past_shutdown_production)

        model.PowerGeneratedStartupShutdown = Expression(model.ThermalGenerators, model.TimePeriods,
                                                         rule=power_generated_startup_shutdown_expr_rule)

    else:
        ## if we're here, then we can use 1-bin models
        ## and no need to do the additional work
        def power_generated_expr_rule(m, g, t):
            linear_vars, linear_coefs = m._get_power_generated_lists(m,g,t)
            linear_expr = get_linear_expr(m.UnitOn)
            return linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs)
        model.PowerGeneratedStartupShutdown = Expression(model.ThermalGenerators, model.TimePeriods,
                                                         rule=power_generated_expr_rule)

## garver/ME power variables (above minimum)
@add_model_attr(component_name, requires = {'data_loader': None, 'status_vars': None})
def garver_power_vars(model):
    '''
    The main variable representing generator output is PowerGeneratedAboveMinimum,
    which is exactly what it says. Originally proposed in

    L. L. Garver. Power generation scheduling by integer programming-development
    of theory. Power Apparatus and Systems, Part III. Transactions of the
    American Institute of Electrical Engineers, 81(3): 730–734, April 1962. ISSN
    0097-2460.
    '''

    # NOTE: this should work with any formulation of the status_vars and data_loader currently

    # amount of power produced by each generator above minimum, at each time period.
    def garver_power_bounds_rule(m, g, t):
        return (0, m.MaximumPowerOutput[g,t]-m.MinimumPowerOutput[g,t])

    model.PowerGeneratedAboveMinimum = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=garver_power_bounds_rule) 

    model._get_power_generated_above_minimum_lists = lambda m,g,t : ([m.PowerGeneratedAboveMinimum[g,t]], [1.])
    model._get_negative_power_generated_above_minimum_lists = lambda m,g,t : ([m.PowerGeneratedAboveMinimum[g,t]], [-1.])
    
    linear_expr = get_linear_expr(model.UnitOn)

    ## Note: these only get used in system balance constraints
    def power_generated_expr_rule(m, g, t):
        #return m.PowerGeneratedAboveMinimum[g,t] + m.MinimumPowerOutput[g,t]*m.UnitOn[g,t]
        return linear_expr( constant=0., linear_vars=[m.PowerGeneratedAboveMinimum[g,t], m.UnitOn[g,t]], \
                                linear_coefs=[1., m.MinimumPowerOutput[g,t]] )
    model.PowerGenerated = Expression(model.ThermalGenerators, model.TimePeriods, rule=power_generated_expr_rule)

    model._get_power_generated_lists = lambda m,g,t : ([m.PowerGeneratedAboveMinimum[g,t], m.UnitOn[g,t]], [1., m.MinimumPowerOutput[g,t]])
    model._get_negative_power_generated_lists = lambda m,g,t : ([m.PowerGeneratedAboveMinimum[g,t], m.UnitOn[g,t]], [-1., -m.MinimumPowerOutput[g,t]])

    _add_power_generated_startup_shutdown(model)

    return

## carrion arroyo power variables (above minimum)
@add_model_attr(component_name, requires = {'data_loader': None, 'status_vars': None})
def basic_power_vars(model):
    '''
    The main power variable represents to the total output of the generator, which
    is the variable that has typically been used in the literature, especially
    prior to Morales-Espana et. al. (2013).
    '''

    # amount of power produced by each generator, at each time period.
    def power_bounds_rule(m, g, t):
        return (0, m.MaximumPowerOutput[g,t])
    model.PowerGenerated = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=power_bounds_rule) 

    model._get_power_generated_lists = lambda m,g,t : ([m.PowerGenerated[g,t]], [1.])
    model._get_negative_power_generated_lists = lambda m,g,t : ([m.PowerGenerated[g,t]], [-1.])
    
    # This allows for automatic substitution in some places, such as PiecewisePorductionSum in production_costs
    def power_generated_expr_rule(m, g, t):
        return m.PowerGenerated[g,t] - m.MinimumPowerOutput[g,t]*m.UnitOn[g,t]
    model.PowerGeneratedAboveMinimum = Expression(model.ThermalGenerators, model.TimePeriods, rule=power_generated_expr_rule)

    model._get_power_generated_above_minimum_lists = lambda m,g,t : ([m.PowerGenerated[g,t], m.UnitOn[g,t]], [1., -m.MinimumPowerOutput[g,t]])
    model._get_negative_power_generated_above_minimum_lists = lambda m,g,t : ([m.PowerGenerated[g,t], m.UnitOn[g,t]], [-1., m.MinimumPowerOutput[g,t]])

    _add_power_generated_startup_shutdown(model)

    return
