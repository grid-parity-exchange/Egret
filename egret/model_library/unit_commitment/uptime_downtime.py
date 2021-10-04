#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## file for minimum uptime/downtime constraints
from pyomo.environ import *
import math

from .uc_utils import add_model_attr, get_linear_expr
component_name = 'uptime_downtime'

def _add_initial(model):

    # constraint due to initial conditions.
    def enforce_up_time_constraints_initial(m, g):
        if value(m.InitialTimePeriodsOnLine[g]) == 0:
            return
        for t in range(m.TimePeriods.first(), value(m.InitialTimePeriodsOnLine[g])+m.TimePeriods.first()):
            if m.status_vars == 'ALS_state_transition_vars':
                m.UnitStayOn[g,t].value = 1
                m.UnitStayOn[g,t].fix()
            else:
                m.UnitOn[g,t].value = 1
                m.UnitOn[g,t].fix()
    
    model.EnforceUpTimeConstraintsInitial = BuildAction(model.ThermalGenerators, rule=enforce_up_time_constraints_initial)

    # constraint due to initial conditions.
    def enforce_down_time_constraints_initial(m, g):
        if value(m.InitialTimePeriodsOffLine[g]) == 0:
            return
        for t in range(m.TimePeriods.first(), value(m.InitialTimePeriodsOffLine[g])+m.TimePeriods.first()):
            if m.status_vars == 'ALS_state_transition_vars':
                m.UnitStayOn[g,t].value = 0
                m.UnitStayOn[g,t].fix()
            else:
                m.UnitOn[g,t].value = 0
                m.UnitOn[g,t].fix()
    
    model.EnforceDownTimeConstraintsInitial = BuildAction(model.ThermalGenerators, rule=enforce_down_time_constraints_initial)

def _add_fixed_and_initial(model):

    # Fixed commitment constraints
    def enforce_fixed_commitments_rule(m,g,t):
        if value(m.FixedCommitment[g,t]) is not None:
            if m.status_vars == 'ALS_state_transition_vars':
                m.UnitStayOn[g,t].value = value(m.FixedCommitment[g,t])
                m.UnitStayOn[g,t].fix()
            else:
                m.UnitOn[g,t].value = value(m.FixedCommitment[g,t])
                m.UnitOn[g,t].fix()
    model.EnforceFixedCommitments = BuildAction(model.ThermalGenerators, model.TimePeriods, rule=enforce_fixed_commitments_rule)

    _add_initial(model)

def _3bin_logic(model):
    
    linear_expr = get_linear_expr(model.UnitOn, model.UnitStart, model.UnitStop)

    initial_time = value(model.InitialTime)
    def logical_rule(m,g,t):
        if t==initial_time:
            return (linear_expr(
                linear_vars=[m.UnitOn[g,t], m.UnitStart[g,t], m.UnitStop[g,t]],
                linear_coefs=[1., -1., 1.],
                ), m.UnitOnT0[g] )
        return (linear_expr(
            linear_vars=[m.UnitOn[g,t], m.UnitOn[g,t-1], m.UnitStart[g,t], m.UnitStop[g,t]],
            linear_coefs=[1., -1, -1., 1.],
            ), 0. )
    
    model.Logical = Constraint(model.ThermalGenerators, model.TimePeriods, rule=logical_rule)

def _2bin_logic(model):

    linear_expr = get_linear_expr(model.UnitOn, model.UnitStart)

    initial_time = value(model.InitialTime)

    def logical_rule(m,g,t):
        if t==initial_time:
            return (None, linear_expr([m.UnitOn[g,t], m.UnitStart[g,t]], [1.,-1.]), m.UnitOnT0[g])
            #return m.UnitOn[g, t] - m.UnitOnT0[g] <= m.UnitStart[g,t]
        return (None, linear_expr([m.UnitOn[g,t], m.UnitOn[g,t-1], m.UnitStart[g,t]], [1.,-1.,-1.]), 0.)
        #return m.UnitOn[g,t] - m.UnitOn[g,t-1] <= m.UnitStart[g,t]

    model.Logical = Constraint(model.ThermalGenerators, model.TimePeriods,rule=logical_rule)

def _ALS_logic(model):

    def logical_rule(m,g,t):
        if t==value(m.InitialTime):
            UnitStartT0 = 1 if (value(m.UnitOnT0State[g]) == 1) else 0
            UnitStayOnT0 = value(m.UnitOnT0[g]) - UnitStartT0
            return m.UnitStayOn[g, t] - UnitStayOnT0 ==  UnitStartT0 - m.UnitStop[g,t]
        return m.UnitStayOn[g,t] - m.UnitStayOn[g,t-1] == m.UnitStart[g,t-1] - m.UnitStop[g,t]
    
    model.Logical = Constraint(model.ThermalGenerators, model.TimePeriods,rule=logical_rule)

## NOTE: if we use the 3-bin variables, we always need the logic constraint
@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['CA_1bin_vars', 'garver_2bin_vars', 'garver_3bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars']
                                            })
def TKW_UT_DT(model):
    '''
    Uptime/downtime constraints from

    Samer Takriti, Benedikt Krasenbrink, and Lilian S-Y Wu. Incorporating fuel
    constraints and electricity spot prices into the stochastic unit commitment
    problem. Operations Research, 48(2):268–280, 2000.
    '''

    _add_fixed_and_initial(model)

    def enforce_up_time_constraints(m, g, t, t_prime):
        if t+1 <= t_prime <= t+value(m.ScaledMinimumUpTime[g])-1:
            if t == m.InitialTime:
                return m.UnitOn[g,t] - m.UnitOnT0[g] <= m.UnitOn[g,t_prime]
            else:
                return m.UnitOn[g,t] - m.UnitOn[g,t-1] <= m.UnitOn[g,t_prime]
        else:
            return Constraint.Skip
    model.UpTime = Constraint(model.ThermalGenerators, model.TimePeriods, model.TimePeriods, rule=enforce_up_time_constraints)

    def enforce_down_time_constraints(m, g, t, t_prime):
        if t+1 <= t_prime <= t+value(m.ScaledMinimumDownTime[g])-1:
            if t == m.InitialTime:
                return m.UnitOnT0[g] - m.UnitOn[g,t] <= 1 - m.UnitOn[g,t_prime]
            else:
                return m.UnitOn[g,t-1] - m.UnitOn[g,t] <= 1 - m.UnitOn[g,t_prime]
        else:
            return Constraint.Skip
    model.DownTime = Constraint(model.ThermalGenerators, model.TimePeriods, model.TimePeriods, rule=enforce_down_time_constraints)

    if model.status_vars in ['garver_3bin_vars', 'garver_3bin_relaxed_stop_vars']:
        _3bin_logic(model)
    if model.status_vars in ['garver_2bin_vars']:
        _2bin_logic(model)
    if model.status_vars in ['ALS_state_transition_vars']:
        _ALS_logic(model)


## NOTE: if we use the 3-bin variables, we always need the logic constraint
@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['CA_1bin_vars', 'garver_2bin_vars', 'garver_3bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars']
                                            })
def CA_UT_DT(model):
    '''
    Uptime/downtime constraints (21)--(26) from

    Carrion, M. and Arroyo, J. (2006) A Computationally Efficient Mixed-Integer
    Liner Formulation for the Thermal Unit Commitment Problem. IEEE Transactions
    on Power Systems, Vol. 21, No. 3, Aug 2006.
    '''

    _add_fixed_and_initial(model)

    # constraint for each time period after that not involving the initial condition.
    def enforce_up_time_constraints_subsequent(m, g, t):
       if t <= value(m.InitialTimePeriodsOnLine[g]):
          # handled by the EnforceUpTimeConstraintInitial constraint.
          return Constraint.Skip
       elif t <= (value(m.NumTimePeriods - m.ScaledMinimumUpTime[g]) + 1):
          # the right-hand side terms below are only positive if the unit was off in the previous time period but on in this one =>
          # the value is the minimum number of subsequent consecutive time periods that the unit is required to be on.
          if t == m.InitialTime:
             return sum(m.UnitOn[g, n] for n in range(t,t+value(m.ScaledMinimumUpTime[g]))) >= \
                    m.ScaledMinimumUpTime[g] * (m.UnitOn[g, t] - m.UnitOnT0[g])
          else:
             return sum(m.UnitOn[g, n] for n in range(t,t+value(m.ScaledMinimumUpTime[g]))) >= \
                    m.ScaledMinimumUpTime[g] * (m.UnitOn[g, t] - m.UnitOn[g, t-1])
       else:
          # handle the final (ScaledMinimumUpTime[g] - 1) time periods - if a unit is started up in
          # this interval, it must remain on-line until the end of the time span.
          if t == m.InitialTime: # can happen when small time horizons are specified
             return sum((m.UnitOn[g, n] - (m.UnitOn[g, t] - m.UnitOnT0[g])) for n in range(t,m.TimePeriods.last()+1)) >= 0.0
          else:
             return sum((m.UnitOn[g, n] - (m.UnitOn[g, t] - m.UnitOn[g, t-1])) for n in range(t, m.TimePeriods.last()+1)) >= 0.0
    
    model.UpTime = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_up_time_constraints_subsequent)

    # constraint for each time period after that not involving the initial condition.
    def enforce_down_time_constraints_subsequent(m, g, t):
       if t <= value(m.InitialTimePeriodsOffLine[g]):
          # handled by the EnforceDownTimeConstraintInitial constraint.
          return Constraint.Skip
       elif t <= (value(m.NumTimePeriods - m.ScaledMinimumDownTime[g]) + 1):
          # the right-hand side terms below are only positive if the unit was off in the previous time period but on in this one =>
          # the value is the minimum number of subsequent consecutive time periods that the unit is required to be on.
          if t == m.InitialTime:
             return sum((1 - m.UnitOn[g, n]) for n in range(t,t+value(m.ScaledMinimumDownTime[g]))) >= \
                    m.ScaledMinimumDownTime[g] * (m.UnitOnT0[g] - m.UnitOn[g, t])
          else:
             return sum((1 - m.UnitOn[g, n]) for n in range(t,t+value(m.ScaledMinimumDownTime[g]))) >= \
                    m.ScaledMinimumDownTime[g] * (m.UnitOn[g, t-1] - m.UnitOn[g, t])
       else:
          # handle the final (ScaledMinimumDownTime[g] - 1) time periods - if a unit is shut down in
          # this interval, it must remain off-line until the end of the time span.
          if t == m.InitialTime: # can happen when small time horizons are specified
             return sum(((1 - m.UnitOn[g, n]) - (m.UnitOnT0[g] - m.UnitOn[g, t])) for n in range(t, m.TimePeriods.last()+1)) >= 0.0
          else:
             return sum(((1 - m.UnitOn[g, n]) - (m.UnitOn[g, t-1] - m.UnitOn[g, t])) for n in range(t, m.TimePeriods.last()+1)) >= 0.0
    
    model.DownTime = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_down_time_constraints_subsequent)

    if model.status_vars in ['garver_3bin_vars', 'garver_3bin_relaxed_stop_vars']:
        _3bin_logic(model)
    if model.status_vars in ['garver_2bin_vars']:
        _2bin_logic(model)
    if model.status_vars in ['ALS_state_transition_vars']:
        _ALS_logic(model)



@add_model_attr(component_name, requires = {'data_loader': None, 
                                            'status_vars': ['garver_3bin_vars', 'garver_2bin_vars', 'garver_3bin_relaxed_stop_vars']
                                            })
def DEKT_UT_DT(model):
    '''
    3-bin uptime downtime constraints originally appearing in

    T.S. Dillon, et. al. Integer programming approach to 
    the problem of optimal unit commitment with 
    probabilistic reserve determination. IEEE Transactions
    Power on Applied Systems, vol PAS-97, no. 6, pp.2154-2166,
    Nov/Dec 1978

    and more recently appearing in 

    J.M. Arroyo and A.J. Conejo, "Optimal responce of a thermal
    unit to an electricity spot market," IEEE Transactions
    on Power Systems, vol. 15, no. 3, pp.1098-1104, Aug. 2000

    These were used as a baseline in the Carrion-Arroyo (2006)
    paper and in Ostrowski et. al (2012).
    '''

    _add_fixed_and_initial(model)

    # constraint for each time period after that not involving the initial condition.
    def enforce_up_time_constraints_subsequent(m, g, t):
        if t <= value(m.InitialTimePeriodsOnLine[g]):
            # handled by the EnforceUpTimeConstraintInitial constraint.
            return Constraint.Skip
        elif t <= (value(m.NumTimePeriods - m.ScaledMinimumUpTime[g]) + 1):
          # the right-hand side terms below are only positive if the unit was off in the previous time period but on in this one =>
          # the value is the minimum number of subsequent consecutive time periods that the unit is required to be on.
            return sum(m.UnitOn[g, n] for n in range(t, t+value(m.ScaledMinimumUpTime[g]))) >= \
                    m.ScaledMinimumUpTime[g] * m.UnitStart[g, t]
        else:
            return sum((m.UnitOn[g, n] - m.UnitStart[g,t]) for n in range(t, m.TimePeriods.last()+1)) >= 0.0
    
    model.UpTime = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_up_time_constraints_subsequent)

    # constraint for each time period after that not involving the initial condition.
    def enforce_down_time_constraints_subsequent(m, g, t):
        if t <= value(m.InitialTimePeriodsOffLine[g]):
            # handled by the EnforceDownTimeConstraintInitial constraint.
             return Constraint.Skip
        elif t <= (value(m.NumTimePeriods - m.ScaledMinimumDownTime[g]) + 1):
            # the right-hand side terms below are only positive if the unit was off in the previous time period but on in this one =>
            # the value is the minimum number of subsequent consecutive time periods that the unit is required to be on.
            return sum((1 - m.UnitOn[g, n]) for n in range(t, t+value(m.ScaledMinimumDownTime[g]))) >= \
                    m.ScaledMinimumDownTime[g] * m.UnitStop[g, t]
        else:
            # handle the final (ScaledMinimumDownTime[g] - 1) time periods - if a unit is shut down in
            # this interval, it must remain off-line until the end of the time span.
            return sum(((1 - m.UnitOn[g, n]) - m.UnitStop[g, t]) for n in range(t, m.TimePeriods.last()+1)) >= 0.0
    
    model.DownTime = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_down_time_constraints_subsequent)

    if model.status_vars in ['garver_3bin_vars', 'garver_3bin_relaxed_stop_vars']:
        _3bin_logic(model)

    elif model.status_vars in ['garver_2bin_vars']:
        _2bin_logic(model)


@add_model_attr(component_name, requires = {'data_loader': None, 'status_vars': ['garver_3bin_vars', 'garver_3bin_relaxed_stop_vars']})
def rajan_takriti_UT_DT(model):
    '''
    Uptime/downtime constraints (3) and (4) from

    D. Rajan and S. Takriti. Minimum up/down polytopes of the unit commitment
    problem with start-up costs. IBM Res. Rep, 2005.
    '''

    _add_fixed_and_initial(model)
    
    #######################
    # up-time constraints #
    #######################
    
    linear_expr = get_linear_expr(model.UnitOn, model.UnitStart, model.UnitStop)

    def uptime_rule(m,g,t):
        if t < value(m.ScaledMinimumUpTime[g]):
            return Constraint.Skip
        linear_vars = [m.UnitStart[g,i] for i in range(t-value(m.ScaledMinimumUpTime[g])+1, t+1)]
        linear_coefs = [1.]*len(linear_vars)
        linear_vars.append(m.UnitOn[g,t])
        linear_coefs.append(-1.)
        return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0.)

    model.UpTime = Constraint(model.ThermalGenerators, model.TimePeriods, rule=uptime_rule)

    #########################
    # down-time constraints #
    #########################

    def downtime_rule(m,g,t):
        if t < value(m.ScaledMinimumDownTime[g]):
            return Constraint.Skip
        linear_vars = [m.UnitStop[g,i] for i in range(t-value(m.ScaledMinimumDownTime[g])+1, t+1)]
        linear_coefs = [1.]*len(linear_vars)
        linear_vars.append(m.UnitOn[g,t])
        linear_coefs.append(1.)
        return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 1.)
    
    model.DownTime = Constraint(model.ThermalGenerators,model.TimePeriods,rule=downtime_rule)
    
    _3bin_logic(model)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_2bin_vars', 'garver_3bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars']})
def rajan_takriti_UT_DT_2bin(model):
    '''
    Uptime/downtime constraints (3) and (5) from

    D. Rajan and S. Takriti. Minimum up/down polytopes of the unit commitment
    problem with start-up costs. IBM Res. Rep, 2005.
    '''

    _add_fixed_and_initial(model)
    
    #######################
    # up-time constraints #
    #######################
    
    def uptime_rule(m,g,t):
        if t < value(m.ScaledMinimumUpTime[g]):
            return Constraint.Skip
        if m.status_vars in ['ALS_state_transition_vars']:
            return sum(m.UnitStart[g,i] for i in range(t-value(m.ScaledMinimumUpTime[g])+1, t)) <= m.UnitStayOn[g,t] 
        else:
            return sum(m.UnitStart[g,i] for i in range(t-value(m.ScaledMinimumUpTime[g])+1, t+1)) <= m.UnitOn[g,t] 
    
    model.UpTime = Constraint(model.ThermalGenerators, model.TimePeriods, rule=uptime_rule)
    
    
    #########################
    # down-time constraints #
    #########################
    
    ## equation (5) and Rajan-Takriti
    def downtime_rule(m,g,t):
        if t < value(m.ScaledMinimumDownTime[g]):
            return Constraint.Skip
        if t == value(m.ScaledMinimumDownTime[g]):
            return sum(m.UnitStart[g,i] for i in range(t-value(m.ScaledMinimumDownTime[g])+1,t+1)) <= 1 - m.UnitOnT0[g]
        else: 
            return sum(m.UnitStart[g,i] for i in range(t-value(m.ScaledMinimumDownTime[g])+1,t+1)) <= 1 - m.UnitOn[g,t-value(m.ScaledMinimumDownTime[g])] 
    
    model.DownTime = Constraint(model.ThermalGenerators,model.TimePeriods,rule=downtime_rule)

    if model.status_vars in ['garver_3bin_vars', 'garver_3bin_relaxed_stop_vars']:
        _3bin_logic(model)

    elif model.status_vars in ['garver_2bin_vars']:
        _2bin_logic(model)

    elif model.status_vars in ['ALS_state_transition_vars']:
        _ALS_logic(model)

    else:
        raise Exception("Couldn't find logic constraint for status_vars "+model.status_vars)
