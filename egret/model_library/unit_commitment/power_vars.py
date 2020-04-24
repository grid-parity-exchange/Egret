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

from .uc_utils import add_model_attr 
component_name = 'power_vars'

def _add_reactive_power_vars(model):

    def reactive_power_bounds_rule(m,g,t):
        return (m.MinimumReactivePowerOutput[g,t], m.MaximumReactivePowerOutput[g,t])
    model.ReactivePowerGenerated = Var(model.ThermalGenerators, model.TimePeriods, within=Reals, bounds=reactive_power_bounds_rule)

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
    
    ## Note: these only get used in system balance constraints
    def power_generated_expr_rule(m, g, t):
        return m.PowerGeneratedAboveMinimum[g,t] + m.MinimumPowerOutput[g,t]*m.UnitOn[g,t]
    model.PowerGenerated = Expression(model.ThermalGenerators, model.TimePeriods, rule=power_generated_expr_rule)

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

    
    # This allows for automatic substitution in some places, such as PiecewisePorductionSum in production_costs
    def power_generated_expr_rule(m, g, t):
        return m.PowerGenerated[g,t] - m.MinimumPowerOutput[g,t]*m.UnitOn[g,t]
    model.PowerGeneratedAboveMinimum = Expression(model.ThermalGenerators, model.TimePeriods, rule=power_generated_expr_rule)

    return


@add_model_attr(component_name, requires = {'data_loader': None, 'status_vars': None})
def rescaled_power_vars(model):
    '''
    Rescalling the power variables as in

    Linfeng Yang, Chen Zhang, Jinbao Jian, Ke Meng, Yan Xu, and Zhaoyang Dong.
    A novel projected two-binary-variable formulation for unit commitment in
    power systems. Applied energy, 187:732–745, 2017.
    '''
    
    model.UnitPowerGeneratedAboveMinimum = Var(model.ThermalGenerators, model.TimePeriods, within=UnitInterval)

    def power_generated_above_min_expr_rule(m, g, t):
        return (m.MaximumPowerOutput[g,t] - m.MinimumPowerOutput[g,t])*m.UnitPowerGeneratedAboveMinimum[g,t]
    model.PowerGeneratedAboveMinimum = Expression(model.ThermalGenerators, model.TimePeriods, rule=power_generated_above_min_expr_rule)

    def power_generated_expr_rule(m, g, t):
        return m.PowerGeneratedAboveMinimum[g,t] + m.MinimumPowerOutput[g,t]*m.UnitOn[g,t]
    model.PowerGenerated = Expression(model.ThermalGenerators, model.TimePeriods, rule=power_generated_expr_rule)
