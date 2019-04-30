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

from .uc_utils import add_model_attr 
component_name = 'power_balance'

#TODO: this doesn't check if storage_services is added first, 
#      but this will only happen when there are storage_services!
@add_model_attr(component_name, requires = {'data_loader': None,
                                            'power_vars': None,
                                            'non_dispatchable_vars': None
                                            })
def power_balance_constraints(model):
    '''
    adds the demand and network constraints to the model
    '''

    # system variables
    # amount of power flowing along each line, at each time period
    def line_power_bounds_rule(m, l, t):
       return (-m.ThermalLimit[l], m.ThermalLimit[l])
    model.LinePower = Var(model.TransmissionLines, model.TimePeriods, bounds=line_power_bounds_rule)
    
    # voltage angles at the buses (S) (lock the first bus at 0) in radians
    model.Angle = Var(model.Buses, model.TimePeriods, within=Reals, bounds=(-3.14159265,3.14159265))
    
    def fix_first_angle_rule(m,t):
        first_bus = list(sorted(m.Buses))[0]
        return m.Angle[first_bus,t] == 0.0
    model.FixFirstAngle = Constraint(model.TimePeriods, rule=fix_first_angle_rule)

    def line_power_rule(m, l, t):
        return m.LinePower[l,t] == (m.Angle[m.BusFrom[l], t] - m.Angle[m.BusTo[l], t]) / m.Impedence[l]
    model.CalculateLinePower = Constraint(model.TransmissionLines, model.TimePeriods, rule=line_power_rule)
    
    def interface_from_limit_rule(m,i,t):
        return sum(m.LinePower[l,t] for l in m.InterfaceLines[i]) <= m.InterfaceFromLimit[i]
    model.InterfaceFromLimitConstr = Constraint(model.Interfaces, model.TimePeriods, rule=interface_from_limit_rule)

    def interface_to_limit_rule(m,i,t):
        return sum(m.LinePower[l,t] for l in m.InterfaceLines[i]) >= -m.InterfaceToLimit[i]
    model.InterfaceToLimitConstr = Constraint(model.Interfaces, model.TimePeriods, rule=interface_to_limit_rule)
    
    
    #####################################################
    # load "shedding" can be both positive and negative #
    #####################################################
    model.LoadGenerateMismatch = Var(model.Buses, model.TimePeriods, within=Reals)
    model.posLoadGenerateMismatch = Var(model.Buses, model.TimePeriods, within=NonNegativeReals) # load shedding
    model.negLoadGenerateMismatch = Var(model.Buses, model.TimePeriods, within=NonNegativeReals) # over generation
    
    def define_pos_neg_load_generate_mismatch_rule(m, b, t):
        return m.posLoadGenerateMismatch[b, t] - m.negLoadGenerateMismatch[b, t] == m.LoadGenerateMismatch[b, t]
    model.DefinePosNegLoadGenerateMismatch = Constraint(model.Buses, model.TimePeriods, rule = define_pos_neg_load_generate_mismatch_rule)

    # the following constraints are necessarily, at least in the case of CPLEX 12.4, to prevent
    # the appearance of load generation mismatch component values in the range of *negative* e-5.
    # what these small negative values do is to cause the optimal objective to be a very large negative,
    # due to obviously large penalty values for under or over-generation. JPW would call this a heuristic
    # at this point, but it does seem to work broadly. we tried a single global constraint, across all
    # buses, but that failed to correct the problem, and caused the solve times to explode.
    
    def pos_load_generate_mismatch_tolerance_rule(m, b):
       return sum((m.posLoadGenerateMismatch[b,t] for t in m.TimePeriods)) >= 0.0
    model.PosLoadGenerateMismatchTolerance = Constraint(model.Buses, rule=pos_load_generate_mismatch_tolerance_rule)
    
    def neg_load_generate_mismatch_tolerance_rule(m, b):
       return sum((m.negLoadGenerateMismatch[b,t] for t in m.TimePeriods)) >= 0.0
    model.NegLoadGenerateMismatchTolerance = Constraint(model.Buses, rule=neg_load_generate_mismatch_tolerance_rule)
    
    
    # Power balance at each node (S)
    def power_balance(m, b, t):
        # bus b, time t (S)
        if m.storage_services:
            return sum((1 - m.GeneratorForcedOutage[g,t]) * m.PowerGenerated[g, t] for g in m.ThermalGeneratorsAtBus[b]) \
                + sum(m.PowerOutputStorage[s, t] for s in m.StorageAtBus[b])\
                - sum(m.PowerInputStorage[s, t] for s in m.StorageAtBus[b])\
                + sum(m.NondispatchablePowerUsed[g, t] for g in m.NondispatchableGeneratorsAtBus[b]) \
                + sum(m.LinePower[l,t] for l in m.LinesTo[b]) \
                - sum(m.LinePower[l,t] for l in m.LinesFrom[b]) \
                + m.LoadGenerateMismatch[b,t] \
                == m.Demand[b, t] 
        else:
            return sum((1 - m.GeneratorForcedOutage[g,t]) * m.PowerGenerated[g, t] for g in m.ThermalGeneratorsAtBus[b]) \
                + sum(m.NondispatchablePowerUsed[g, t] for g in m.NondispatchableGeneratorsAtBus[b]) \
                + sum(m.LinePower[l,t] for l in m.LinesTo[b]) \
                - sum(m.LinePower[l,t] for l in m.LinesFrom[b]) \
                + m.LoadGenerateMismatch[b,t] \
                == m.Demand[b, t] 
    
    model.PowerBalance = Constraint(model.Buses, model.TimePeriods, rule=power_balance)

    return
## end power_balance_constraints
