#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## file for startup cost formulations
from pyomo.environ import *
import math

from .uc_utils import add_model_attr 
from .status_vars import _is_relaxed
component_name = 'startup_costs'

## NOTE: the add_startup_cost_var is necessary in a simulation context
@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars','garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            })
def KOW_startup_costs(model, add_startup_cost_var=True):
    '''
    Start-up cost formulation "Match" from

    Ben Knueven, Jim Ostrowski, and Jean-Paul Watson. A novel matching
    formulation for startup costs in unit commitment, 2017.
    URL http://www.optimization-online.org/DB_FILE/2017/03/5897.pdf.
    '''

    #begin ostrowski startup costs
    def ValidShutdownTimePeriods_generator(m,g):
        ## for speed, if we don't have different startups
        if len(m.ScaledStartupLags[g]) <= 1:
            return []
                                        ## adds the necessary index for starting-up after a shutdown before the time horizon began
        return (t for t in (list(m.TimePeriods)+([] if (value(m.UnitOnT0State[g]) >= 0) else [m.InitialTime + int(round(value(m.UnitOnT0State[g]/value(m.TimePeriodLengthHours))))])))
    model.ValidShutdownTimePeriods=Set(model.ThermalGenerators, initialize=ValidShutdownTimePeriods_generator)
    
    def ShutdownHotStartupPairs_generator(m,g):
        ## for speed, if we don't have different startups
        if len(m.ScaledStartupLags[g]) <= 1:
            return [] 
        return ((t_prime, t) for t_prime in m.ValidShutdownTimePeriods[g] for t in m.TimePeriods if (m.ScaledStartupLags[g].first() <= t - t_prime < m.ScaledStartupLags[g].last()))
    model.ShutdownHotStartupPairs = Set(model.ThermalGenerators, initialize=ShutdownHotStartupPairs_generator, dimen=2)
    
    # (g,t',t) will be an inidicator for g for shutting down at time t' and starting up at time t
    def StartupIndicator_domain_generator(m):
        return ((g,t_prime,t) for g in m.ThermalGenerators for t_prime,t in m.ShutdownHotStartupPairs[g]) 
    model.StartupIndicator_domain=Set(initialize=StartupIndicator_domain_generator, dimen=3)
    
    if _is_relaxed(model):
        model.StartupIndicator=Var(model.StartupIndicator_domain, within=UnitInterval)
    else:
        model.StartupIndicator=Var(model.StartupIndicator_domain, within=Binary)

    ############################################################
    # compute the per-generator, per-time period startup costs #
    ############################################################
    
    def startup_match_rule(m, g, t):
        return sum(m.StartupIndicator[g, t_prime, s] for (t_prime, s) in m.ShutdownHotStartupPairs[g] if s == t) <= m.UnitStart[g,t]
    model.StartupMatch = Constraint(model.ThermalGenerators, model.TimePeriods, rule=startup_match_rule)
    
    def GeneratorShutdownPeriods_generator(m):
        return ((g,t) for g in m.ThermalGenerators for t in m.ValidShutdownTimePeriods[g])
    model.GeneratorShutdownPeriods = Set(initialize=GeneratorShutdownPeriods_generator, dimen=2)
    
    def shutdown_match_rule(m, g, t):
        if t < m.InitialTime:
            begin_pairs = [(s, t_prime) for (s, t_prime) in m.ShutdownHotStartupPairs[g] if s == t]
            if not begin_pairs: ##if this is empty
                return Constraint.Feasible
            else:
                return sum(m.StartupIndicator[g, s, t_prime] for (s, t_prime) in begin_pairs) <= 1
        else:
            return sum(m.StartupIndicator[g, s, t_prime] for (s, t_prime) in m.ShutdownHotStartupPairs[g] if s == t) <= m.UnitStop[g,t]
    model.ShutdownMatch = Constraint(model.GeneratorShutdownPeriods, rule=shutdown_match_rule)

    if add_startup_cost_var:
        model.StartupCost = Var(model.SingleFuelGenerators, model.TimePeriods, within=Reals)
    
    def ComputeStartupCost2_rule(m,g,t):
        return m.StartupCost[g,t] == m.StartupCosts[g].last()*m.UnitStart[g,t] + \
                                      sum( (list(m.StartupCosts[g])[s-1] - m.StartupCosts[g].last()) * \
                                         sum( m.StartupIndicator[g,tp,t] for tp in m.ValidShutdownTimePeriods[g] \
                                           if (list(m.ScaledStartupLags[g])[s-1] <= t - tp < (list(m.ScaledStartupLags[g])[s])) ) \
                                         for s in m.StartupCostIndices[g] if s < len(m.StartupCostIndices[g]))
    
    model.ComputeStartupCosts=Constraint(model.SingleFuelGenerators, model.TimePeriods, rule=ComputeStartupCost2_rule)

    return


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_3bin_relaxed_stop_vars','garver_2bin_vars' , 'ALS_state_transition_vars'],
                                            })
def MLR_startup_costs(model, add_startup_cost_var=True):
    '''
    Start-up cost formulation in equations (2)--(3) from

    G. Morales-Espana, J. M. Latorre, and A. Ramos. Tight and compact MILP
    formulation for the thermal unit commitment problem. IEEE Transactions on
    Power Systems, 28(4):4897–4908, 2013.
    '''

    #begin latorre startup costs
    
    ############################################################
    # compute the per-generator, per-time period startup costs #
    ############################################################
    ## BK -- replace with delta's
    def startup_costs_index_set_generator(m):
        return ((g,s,t) for t in m.TimePeriods for g in m.ThermalGenerators for s in m.StartupCostIndices[g])
    
    model.StartupCostsIndexSet = Set(initialize=startup_costs_index_set_generator, dimen=3)
    
    if _is_relaxed(model):
        model.delta=Var(model.StartupCostsIndexSet, within=UnitInterval)
    else:
        model.delta=Var(model.StartupCostsIndexSet, within=Binary)

    def delta_eq_rule(m,g,t):
        return sum(m.delta[g,s,t] for s in m.StartupCostIndices[g])==m.UnitStart[g,t]
    
    model.delta_eq=Constraint(model.ThermalGenerators, model.TimePeriods, rule=delta_eq_rule)
    
    ## BK -- updated to reflect previous generator condition
    ## BK -- assumes initial time is 1
    def delta_ineq_rule(m,g,s,t):
        assert(m.InitialTime == 1)
        ## the last startup indicator doesn't have a lag
        if s == len(m.StartupCostIndices[g]):
            return Constraint.Skip
        this_lag = list(m.ScaledStartupLags[g])[s-1]
        next_lag = list(m.ScaledStartupLags[g])[s]
        ## this is negative if the generator has previously been off
        generator_t0_state = int(round(value(m.UnitOnT0State[g])/value(m.TimePeriodLengthHours)))
        
        ## equation (15) in ME
        if next_lag + generator_t0_state < t < next_lag:
            return m.delta[g,s,t] == 0

        if m.status_vars == 'garver_2bin_vars':
            if t == next_lag:
                return m.delta[g,s,t] <= m.UnitOnT0[g] - m.UnitOn[g,t-this_lag] + sum(m.UnitStart[g,t-i] for i in range(this_lag, next_lag))
            elif t > next_lag:
                return m.delta[g,s,t] <= m.UnitOn[g,t-next_lag] - m.UnitOn[g,t-this_lag] + sum(m.UnitStart[g,t-i] for i in range(this_lag, next_lag))
        else:
            ## equation (2) in ME
            if t >= next_lag:
                return m.delta[g,s,t] <= sum(m.UnitStop[g,t-i] for i in range(this_lag, next_lag))
    
        return Constraint.Skip
    
    model.delta_ineq=Constraint(model.StartupCostsIndexSet, rule=delta_ineq_rule)

    if add_startup_cost_var:
        model.StartupCost = Var(model.SingleFuelGenerators, model.TimePeriods, within=Reals)
    
    def ComputeStartupCost2_rule(m,g,t):
        return m.StartupCost[g,t] == sum(m.delta[g,s,t]*list(m.StartupCosts[g])[s-1] for s in m.StartupCostIndices[g])
    
    model.ComputeStartupCosts=Constraint(model.SingleFuelGenerators, model.TimePeriods, rule=ComputeStartupCost2_rule)

    return


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars', 'garver_3bin_relaxed_stop_vars', 'garver_2bin_vars', 'ALS_state_transition_vars'],
                                            })
def KOW_3bin_startup_costs(model, add_startup_cost_var=True):
    '''
    Start-up cost formulation "3bin" from

    Ben Knueven, Jim Ostrowski, and Jean-Paul Watson. A novel matching
    formulation for startup costs in unit commitment, 2017.
    URL http://www.optimization-online.org/DB_FILE/2017/03/5897.pdf.
    '''

    if add_startup_cost_var:
        model.StartupCost = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)
    
    ############################################################
    # compute the per-generator, per-time period startup costs #
    ############################################################
    
    def startup_costs_index_set_generator(m):
       return ((g,t,i) for t in m.TimePeriods for g in m.ThermalGenerators for i in m.StartupCostIndices[g])
    
    model.StartupCostsIndexSet = Set(initialize=startup_costs_index_set_generator, dimen=3)
    
    def compute_startup_costs_rule(m, g, t, i):
        # irios, Nov 18th: I had to change this because the last version didn't work with the last update in coopr.
        this_lag = list(m.ScaledStartupLags[g])[i-1]
        this_cost = list(m.StartupCosts[g])[i-1]
    
        startup_lags = list(m.ScaledStartupLags[g])
        startup_costs = list(m.StartupCosts[g])
        
        generator_t0_state = int(round(value(m.UnitOnT0State[g])/value(m.TimePeriodLengthHours)))
    
        # if the lag involves time periods preceding t=1, then we need to determine if the T0 
        # state is consistent with the lag - if not, we can skip generation of the constraint.
        assert(m.InitialTime == 1)
        if this_lag >= t:
            if generator_t0_state >= 0:
                # the unit has been on - we can't meet the target lag.
                return Constraint.Skip
    
            time_diff = this_lag - t + 1
            if (-generator_t0_state) < time_diff:
                # the generator has not been off for a sufficient number of time periods.
                return Constraint.Skip
         
        
        # can only "look back" in terms of UnitOn variable state (t-1) or this_lag time periods - whichever is smallest.
        # the rest of the time period is captured in the unit T0 state, and is handled in the logic above.
        if m.status_vars == 'garver_2bin_vars':
            ## This is just the equation below with the UnitStop's projected out
            return m.StartupCost[g, t] >= this_cost * m.UnitStart[g, t] \
                                        - sum((this_cost - startup_costs[j])*
                                               ( (m.UnitOn[g,t-startup_lags[j+1]] if t-startup_lags[j+1] > 0 else ( m.UnitOnT0[g] if t-startup_lags[j+1] == 0 else 0. ))
                                                - (m.UnitOn[g,t-startup_lags[j]] if t-startup_lags[j] > 0 else ( m.UnitOnT0[g] if t-startup_lags[j] == 0 else 0. ))  
                                                + sum(m.UnitStart[g, t-k] for k in range(startup_lags[j], startup_lags[j+1]) if k < t ) )
                                          for j in range(0, i-1) )
        else:
            return m.StartupCost[g, t] >= this_cost * m.UnitStart[g, t] \
                                        - sum(
                                               sum(
                                                    (this_cost - startup_costs[j])*m.UnitStop[g, t-k] 
                                               for k in range(startup_lags[j], startup_lags[j+1]) if k < t )
                                          for j in range(0, i-1) )
    
    model.ComputeStartupCosts = Constraint(model.StartupCostsIndexSet, rule=compute_startup_costs_rule)

    return


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            })
def CA_SHB_startup_costs(model, add_startup_cost_var=True):
    '''
    Tightened Carrion-Arroyo start-up costs from

    Matthias Silbernagl, Matthias Huber, and Rene Brandenberg. Improving
    accuracy and efficiency of start-up cost formulations in MIP unit commitment
    by modeling power plant temperatures. IEEE Trans. Power Syst.,
    31(4):2578–2586, 2016.

    equation (13)
    '''

    if add_startup_cost_var:
        model.StartupCost = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)
    
    ############################################################
    # compute the per-generator, per-time period startup costs #
    ############################################################
    
    def startup_costs_index_set_generator(m):
       return ((g,t,i) for t in m.TimePeriods for g in m.ThermalGenerators for i in m.StartupCostIndices[g])
    
    model.StartupCostsIndexSet = Set(initialize=startup_costs_index_set_generator, dimen=3)
    
    def compute_startup_costs_rule(m, g, t, i):
       # irios, Nov 18th: I had to change this because the last version didn't work with the last update in coopr.
       this_lag = list(m.ScaledStartupLags[g])[i-1]
       this_cost = list(m.StartupCosts[g])[i-1]
    
       startup_lags = list(m.ScaledStartupLags[g])
       startup_costs = list(m.StartupCosts[g])
       
       generator_t0_state = int(round(value(m.UnitOnT0State[g])/value(m.TimePeriodLengthHours)))
    
       # if the lag involves time periods preceding t=1, then we need to determine if the T0 
       # state is consistent with the lag - if not, we can skip generation of the constraint.
       assert(m.InitialTime == 1)
       if this_lag >= t:
    
          if generator_t0_state >= 0:
             # the unit has been on - we can't meet the target lag.
             return Constraint.Skip
    
          time_diff = this_lag - t + 1
          if (-generator_t0_state) < time_diff:
             # the generator has not been off for a sufficient number of time periods.
              return Constraint.Skip
        
       
       # can only "look back" in terms of UnitOn variable state (t-1) or this_lag time periods - whichever is smallest.
       # the rest of the time period is captured in the unit T0 state, and is handled in the logic above.
       return m.StartupCost[g, t] >= this_cost * m.UnitOn[g, t] - sum( (this_cost)*m.UnitOn[g, t - k] for k in range(1, min(t, startup_lags[0]+1))) \
                                                                - sum( sum( (this_cost - startup_costs[j])*m.UnitOn[g, t-k] for k in range(startup_lags[j]+1, startup_lags[j+1]+1) if k < t )\
                                                                        for j in range(0, i-1) )
    
    model.ComputeStartupCosts = Constraint(model.StartupCostsIndexSet, rule=compute_startup_costs_rule)

    return


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            })
def CA_startup_costs(model, add_startup_cost_var=True):
    '''
    Equations (12) and (13) from

    Carrion, M. and Arroyo, J. (2006) A Computationally Efficient Mixed-Integer
    Liner Formulation for the Thermal Unit Commitment Problem. IEEE Transactions
    on Power Systems, Vol. 21, No. 3, Aug 2006.
    '''

    if add_startup_cost_var:
        model.StartupCost = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)
    
    ############################################################
    # compute the per-generator, per-time period startup costs #
    ############################################################
    
    def startup_costs_index_set_generator(m):
       return ((g,t,i) for t in m.TimePeriods for g in m.ThermalGenerators for i in m.StartupCostIndices[g])
    
    model.StartupCostsIndexSet = Set(initialize=startup_costs_index_set_generator, dimen=3)
    
    def compute_startup_costs_rule(m, g, t, i):
       # irios, Nov 18th: I had to change this because the last version didn't work with the last update in coopr.
       this_lag = list(m.ScaledStartupLags[g])[i-1]
       this_cost = list(m.StartupCosts[g])[i-1]
    
       generator_t0_state = int(round(value(m.UnitOnT0State[g])/value(m.TimePeriodLengthHours)))
    
       # if the lag involves time periods preceding t=1, then we need to determine if the T0 
       # state is consistent with the lag - if not, we can skip generation of the constraint.
       assert(m.InitialTime == 1)
       if this_lag >= t:
    
          if generator_t0_state >= 0:
             # the unit has been on - we can't meet the target lag.
             return Constraint.Skip
    
          time_diff = this_lag - t + 1
          if (-generator_t0_state) < time_diff:
             # the generator has not been off for a sufficient number of time periods.
              return Constraint.Skip
    
       # can only "look back" in terms of UnitOn variable state (t-1) or this_lag time periods - whichever is smallest.
       # the rest of the time period is captured in the unit T0 state, and is handled in the logic above.
       return m.StartupCost[g, t] >= this_cost * (m.UnitOn[g, t] - sum(m.UnitOn[g, t - k] for k in range(1, min(t, this_lag+1))))
    
    model.ComputeStartupCosts = Constraint(model.StartupCostsIndexSet, rule=compute_startup_costs_rule)

    return

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars', 'garver_3bin_relaxed_stop_vars', 'garver_2bin_vars', 'ALS_state_transition_vars'],
                                            })
def ALS_startup_costs(model, add_startup_cost_var=True):
    '''
    Equation (19) from

    Semih Atakan, Guglielmo Lulli, and Suvrajeet Sen. A state transition MIP
    formulation for the unit commitment problem. IEEE Transactions on Power
    Systems, 33(1):736–748, 2018.

    '''
    model.StartupCostOverHot = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)
    
    ############################################################
    # compute the per-generator, per-time period startup costs #
    ############################################################
    
    def startup_costs_index_set_generator(m):
       return ((g,t,i) for t in m.TimePeriods for g in m.ThermalGenerators for i in m.StartupCostIndices[g])
    
    model.StartupCostsIndexSet = Set(initialize=startup_costs_index_set_generator, dimen=3)
    
    def compute_startup_costs_rule(m, g, t, i):

        ## we don't need these for the first startup cost
        if i == 1:
            return Constraint.Skip
        # irios, Nov 18th: I had to change this because the last version didn't work with the last update in coopr.
        this_lag = list(m.ScaledStartupLags[g])[i-1]
        this_cost = list(m.StartupCosts[g])[i-1]
    
        startup_lags = list(m.ScaledStartupLags[g])
        startup_costs = list(m.StartupCosts[g])
        
        generator_t0_state = int(round(value(m.UnitOnT0State[g])/value(m.TimePeriodLengthHours)))
    
        # if the lag involves time periods preceding t=1, then we need to determine if the T0 
        # state is consistent with the lag - if not, we can skip generation of the constraint.
        assert(m.InitialTime == 1)
        if this_lag >= t:
            if generator_t0_state >= 0:
                # the unit has been on - we can't meet the target lag.
                return Constraint.Skip
    
            time_diff = this_lag - t + 1
            if (-generator_t0_state) < time_diff:
                # the generator has not been off for a sufficient number of time periods.
                return Constraint.Skip
         
        
        # can only "look back" in terms of UnitOn variable state (t-1) or this_lag time periods - whichever is smallest.
        # the rest of the time period is captured in the unit T0 state, and is handled in the logic above.
        if t-this_lag == 0:
            return m.StartupCostOverHot[g, t] >= (this_cost - m.StartupCosts[g].first())*(
                                        m.UnitStart[g, t] \
                                        - sum( m.UnitStart[g,t-i] for i in range(m.ScaledStartupLags[g].first(), this_lag) ) \
                                        - m.UnitOnT0[g] )
        elif t-this_lag < 0:
            return m.StartupCostOverHot[g, t] >= (this_cost - m.StartupCosts[g].first())*(
                                        m.UnitStart[g, t] \
                                        - sum( m.UnitStart[g,t-i] for i in range(m.ScaledStartupLags[g].first(), min(t,this_lag)) ) )
        else:
            return m.StartupCostOverHot[g, t] >= (this_cost - m.StartupCosts[g].first())*(
                                        m.UnitStart[g, t] \
                                        - sum( m.UnitStart[g,t-i] for i in range(m.ScaledStartupLags[g].first(), this_lag) ) \
                                        - m.UnitOn[g,t-this_lag] )
    
    model.ComputeStartupCostsOverHot = Constraint(model.StartupCostsIndexSet, rule=compute_startup_costs_rule)

    if add_startup_cost_var:
        model.StartupCost = Var( model.ThermalGenerators, model.TimePeriods, within=Reals)

    def compute_startup_costs_expr_rule(m, g, t):
        return m.StartupCost[g,t] == m.StartupCostOverHot[g,t] + m.StartupCosts[g].first()*m.UnitStart[g,t]
    model.ComputeStartupCosts = Constraint( model.ThermalGenerators, model.TimePeriods, rule=compute_startup_costs_expr_rule)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars', 'garver_3bin_relaxed_stop_vars', 'garver_2bin_vars', 'ALS_state_transition_vars'],
                                            })
def YZJMXD_startup_costs(model, add_startup_cost_var=True):
    '''
    Equations (34) from

    Linfeng Yang, Chen Zhang, Jinbao Jian, Ke Meng, Yan Xu, and Zhaoyang Dong.
    A novel projected two-binary-variable formulation for unit commitment in
    power systems. Applied energy, 187:732–745, 2017.
    '''

    model.StartupCostOverHot = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)
    
    ############################################################
    # compute the per-generator, per-time period startup costs #
    ############################################################
    
    def startup_costs_index_set_generator(m):
       return ((g,t,i) for t in m.TimePeriods for g in m.ThermalGenerators for i in m.StartupCostIndices[g])
    
    model.StartupCostsIndexSet = Set(initialize=startup_costs_index_set_generator, dimen=3)
    
    def compute_startup_costs_rule(m, g, t, i):

        ## we don't need these for the first startup cost
        if i == 1:
            return Constraint.Skip
        # irios, Nov 18th: I had to change this because the last version didn't work with the last update in coopr.
        this_lag = list(m.ScaledStartupLags[g])[i-1]
        this_cost = list(m.StartupCosts[g])[i-1]
    
        startup_lags = list(m.ScaledStartupLags[g])
        startup_costs = list(m.StartupCosts[g])
        
        generator_t0_state = int(round(value(m.UnitOnT0State[g])/value(m.TimePeriodLengthHours)))
    
        # if the lag involves time periods preceding t=1, then we need to determine if the T0 
        # state is consistent with the lag - if not, we can skip generation of the constraint.
        assert(m.InitialTime == 1)
        if this_lag >= t:
            if generator_t0_state >= 0:
                # the unit has been on - we can't meet the target lag.
                return Constraint.Skip
    
            time_diff = this_lag - t + 1
            if (-generator_t0_state) < time_diff:
                # the generator has not been off for a sufficient number of time periods.
                return Constraint.Skip
         
        
        # can only "look back" in terms of UnitOn variable state (t-1) or this_lag time periods - whichever is smallest.
        # the rest of the time period is captured in the unit T0 state, and is handled in the logic above.
        return m.StartupCostOverHot[g, t] >= (this_cost - m.StartupCosts[g].first())*(
                                        m.UnitStart[g, t] \
                                        - sum( m.UnitOn[g,t-i] for i in range(1, min(this_lag+1,t)) ) ) 
    
    model.ComputeStartupCostsOverHot = Constraint(model.StartupCostsIndexSet, rule=compute_startup_costs_rule)

    if add_startup_cost_var:
        model.StartupCost = Var( model.ThermalGenerators, model.TimePeriods, within=Reals)

    def compute_startup_costs_expr_rule(m, g, t):
        return m.StartupCost[g,t] ==  m.StartupCostOverHot[g,t] + m.StartupCosts[g].first()*m.UnitStart[g,t]
    model.ComputeStartupCosts = Constraint( model.ThermalGenerators, model.TimePeriods, rule=compute_startup_costs_expr_rule)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars', 'garver_3bin_relaxed_stop_vars', 'garver_2bin_vars', 'ALS_state_transition_vars'],
                                            })
def KOW_3bin_startup_costs2(model, add_startup_cost_var=True):
    '''
    This is KOW_3bin_startup_costs but like Atakan et. al. and Yang et. al.,
    we eliminate the constraints associated with the hottest start and
    consider a hot start directly in the objective function
    '''

    model.StartupCostOverHot = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)
    
    ############################################################
    # compute the per-generator, per-time period startup costs #
    ############################################################
    
    def startup_costs_index_set_generator(m):
       return ((g,t,i) for t in m.TimePeriods for g in m.ThermalGenerators for i in m.StartupCostIndices[g])
    
    model.StartupCostsIndexSet = Set(initialize=startup_costs_index_set_generator, dimen=3)
    
    def compute_startup_costs_rule(m, g, t, i):
         ## we don't need these for the first startup cost
        if i == 1:
            return Constraint.Skip
        # irios, Nov 18th: I had to change this because the last version didn't work with the last update in coopr.
        this_lag = list(m.ScaledStartupLags[g])[i-1]
        this_cost = list(m.StartupCosts[g])[i-1]
    
        startup_lags = list(m.ScaledStartupLags[g])
        startup_costs = list(m.StartupCosts[g])
        
        generator_t0_state = int(round(value(m.UnitOnT0State[g])/value(m.TimePeriodLengthHours)))
    
        # if the lag involves time periods preceding t=1, then we need to determine if the T0 
        # state is consistent with the lag - if not, we can skip generation of the constraint.
        assert(m.InitialTime == 1)
        if this_lag >= t:
            if generator_t0_state >= 0:
                # the unit has been on - we can't meet the target lag.
                return Constraint.Skip
    
            time_diff = this_lag - t + 1
            if (-generator_t0_state) < time_diff:
                # the generator has not been off for a sufficient number of time periods.
                return Constraint.Skip
         
        
        # can only "look back" in terms of UnitOn variable state (t-1) or this_lag time periods - whichever is smallest.
        # the rest of the time period is captured in the unit T0 state, and is handled in the logic above.
        if m.status_vars == 'garver_2bin_vars':
            ## This is just the equation below with the UnitStop's projected out
            return m.StartupCostOverHot[g, t] >= (this_cost -m.StartupCosts[g].first()) * m.UnitStart[g, t] \
                                        - sum((this_cost - startup_costs[j])*
                                               ( (m.UnitOn[g,t-startup_lags[j+1]] if t-startup_lags[j+1] > 0 else ( m.UnitOnT0[g] if t-startup_lags[j+1] == 0 else 0. ))
                                                - (m.UnitOn[g,t-startup_lags[j]] if t-startup_lags[j] > 0 else ( m.UnitOnT0[g] if t-startup_lags[j] == 0 else 0. ))  
                                                + sum(m.UnitStart[g, t-k] for k in range(startup_lags[j], startup_lags[j+1]) if k < t ) )
                                          for j in range(0, i-1) )
        else:
            return m.StartupCostOverHot[g, t] >= (this_cost -m.StartupCosts[g].first())* m.UnitStart[g, t] \
                                        - sum(
                                               sum(
                                                    (this_cost - startup_costs[j])*m.UnitStop[g, t-k] 
                                               for k in range(startup_lags[j], startup_lags[j+1]) if k < t )
                                          for j in range(0, i-1) )
    
    model.ComputeStartupCostsOverHot = Constraint(model.StartupCostsIndexSet, rule=compute_startup_costs_rule)

    if add_startup_cost_var:
        model.StartupCost = Var( model.ThermalGenerators, model.TimePeriods, within=Reals)

    def compute_startup_costs_expr_rule(m, g, t):
        return m.StartupCost[g,t] ==  m.StartupCostOverHot[g,t] + m.StartupCosts[g].first()*m.UnitStart[g,t]
    model.ComputeStartupCosts = Constraint( model.ThermalGenerators, model.TimePeriods, rule=compute_startup_costs_expr_rule)

    return


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_3bin_relaxed_stop_vars','garver_2bin_vars' , 'ALS_state_transition_vars'],
                                            })
def MLR_startup_costs2(model, add_startup_cost_var=True):
    '''
    This is MLR startup_costs but like KOW_startup_costs, puts the coldest
    start as the coeffiencent on the start-up variable and gives the deltas
    negativate coefficients for the hot-starts, thus eliminating a variable.
    '''

    #begin latorre startup costs
    
    ############################################################
    # compute the per-generator, per-time period startup costs #
    ############################################################
    ## BK -- replace with delta's
    def startup_costs_index_set_generator(m):
        return ((g,s,t) for t in m.TimePeriods for g in m.ThermalGenerators for s in m.StartupCostIndices[g] if s < len(m.StartupCostIndices[g]))
    
    model.StartupCostsIndexSet = Set(initialize=startup_costs_index_set_generator, dimen=3)
    
    if _is_relaxed(model):
        model.delta=Var(model.StartupCostsIndexSet, within=UnitInterval)
    else:
        model.delta=Var(model.StartupCostsIndexSet, within=Binary)

    def delta_eq_rule(m,g,t):
        return sum(m.delta[g,s,t] for s in m.StartupCostIndices[g] if s < len(m.StartupCostIndices[g])) <= m.UnitStart[g,t]
    
    model.delta_eq=Constraint(model.ThermalGenerators, model.TimePeriods, rule=delta_eq_rule)
    
    ## BK -- updated to reflect previous generator condition
    ## BK -- assumes initial time is 1
    def delta_ineq_rule(m,g,s,t):
        ## the last startup indicator doesn't have a lag
        assert(m.InitialTime == 1)
        if s == len(m.StartupCostIndices[g]):
            return Constraint.Skip
        this_lag = list(m.ScaledStartupLags[g])[s-1]
        next_lag = list(m.ScaledStartupLags[g])[s]
        ## this is negative if the generator has previously been off
        generator_t0_state = int(round(value(m.UnitOnT0State[g])/value(m.TimePeriodLengthHours)))
        
        ## equation (15) in ME
        if next_lag + generator_t0_state < t < next_lag:
            return m.delta[g,s,t] == 0

        if m.status_vars == 'garver_2bin_vars':
            if t == next_lag:
                return m.delta[g,s,t] <= m.UnitOnT0[g] - m.UnitOn[g,t-this_lag] + sum(m.UnitStart[g,t-i] for i in range(this_lag, next_lag))
            elif t > next_lag:
                return m.delta[g,s,t] <= m.UnitOn[g,t-next_lag] - m.UnitOn[g,t-this_lag] + sum(m.UnitStart[g,t-i] for i in range(this_lag, next_lag))
        else:
            ## equation (2) in ME
            if t >= next_lag:
                return m.delta[g,s,t] <= sum(m.UnitStop[g,t-i] for i in range(this_lag, next_lag))
    
        return Constraint.Skip
    
    model.delta_ineq=Constraint(model.StartupCostsIndexSet, rule=delta_ineq_rule)
    
    if add_startup_cost_var:
        model.StartupCost = Var( model.SingleFuelGenerators, model.TimePeriods, within=Reals)

    def ComputeStartupCost2_rule(m,g,t):
        return m.StartupCost[g,t] ==  m.StartupCosts[g].last()*m.UnitStart[g,t] + \
                sum(m.delta[g,s,t]*(list(m.StartupCosts[g])[s-1] - m.StartupCosts[g].last()) for s in m.StartupCostIndices[g] if s < len(m.StartupCostIndices[g]))
    model.ComputeStartupCosts=Constraint(model.SingleFuelGenerators, model.TimePeriods, rule=ComputeStartupCost2_rule)

    return

