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

from .uc_utils import add_model_attr, get_linear_expr
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
    time_period_list = list(model.TimePeriods)
    initial_time = model.TimePeriods.first()
    after_last_time = model.TimePeriods.last()+1
    def ValidShutdownTimePeriods_generator(m,g):
        ## for speed, if we don't have different startups
        if len(m.ScaledStartupLags[g]) <= 1:
            return []
        ## adds the necessary index for starting-up after a shutdown before the time horizon began
        unit_on_t0_state = value(m.UnitOnT0State[g])
        if unit_on_t0_state >= 0:
            return time_period_list
        else:
            return time_period_list+[initial_time + int(round(unit_on_t0_state/value(m.TimePeriodLengthHours)))]
    model.ValidShutdownTimePeriods=Set(model.ThermalGenerators, initialize=ValidShutdownTimePeriods_generator)
    
    def ShutdownHotStartupPairs_generator(m,g):
        ## for speed, if we don't have different startups
        if len(m.ScaledStartupLags[g]) <= 1:
            return
        first_lag = m.ScaledStartupLags[g].first()
        last_lag = m.ScaledStartupLags[g].last()
        for t_prime in m.ValidShutdownTimePeriods[g]:
            t_first = first_lag+t_prime
            t_last = last_lag+t_prime
            if t_first < initial_time:
                t_first = initial_time
            if t_last > after_last_time:
                t_last = after_last_time
            for t in range(t_first, t_last):
                yield (t_prime, t)
    model.ShutdownHotStartupPairs = Set(model.ThermalGenerators, initialize=ShutdownHotStartupPairs_generator, dimen=2)
    
    # (g,t',t) will be an inidicator for g for shutting down at time t' and starting up at time t
    def StartupIndicator_domain_generator(m):
        return ((g,t_prime,t) for g in m.ThermalGenerators for t_prime,t in m.ShutdownHotStartupPairs[g]) 
    model.StartupIndicator_domain=Set(initialize=StartupIndicator_domain_generator, dimen=3)

    if _is_relaxed(model):
        model.StartupIndicator=Var(model.StartupIndicator_domain, within=UnitInterval)
    else:
        model.StartupIndicator=Var(model.StartupIndicator_domain, within=Binary)

    linear_expr = get_linear_expr(model.UnitStart, model.UnitStop)

    ############################################################
    # compute the per-generator, per-time period startup costs #
    ############################################################

    def GeneratorShutdownPeriods_generator(m):
        return ((g,t) for g in m.ThermalGenerators for t in m.ValidShutdownTimePeriods[g])
    model.GeneratorShutdownPeriods = Set(initialize=GeneratorShutdownPeriods_generator, dimen=2)

    model.ShutdownsByStartups = Set(model.ThermalGenerators, model.TimePeriods)
    model.StartupsByShutdowns = Set(model.GeneratorShutdownPeriods)
    for g,t_p,t in model.StartupIndicator_domain:
        model.ShutdownsByStartups[g,t].add(t_p)
        model.StartupsByShutdowns[g,t_p].add(t)

    def startup_match_rule(m, g, t):
        linear_vars = list(m.StartupIndicator[g, t_prime, t] for t_prime in m.ShutdownsByStartups[g,t])
        linear_coefs = [1.]*len(linear_vars)
        linear_vars.append(m.UnitStart[g,t])
        linear_coefs.append(-1.)
        return (None, linear_expr(linear_coefs=linear_coefs, linear_vars=linear_vars), 0.)
    model.StartupMatch = Constraint(model.ThermalGenerators, model.TimePeriods, rule=startup_match_rule)

    def shutdown_match_rule(m, g, t):
        begin_times = m.StartupsByShutdowns[g,t]
        if not begin_times: ##if this is empty
            return Constraint.Feasible
        linear_vars = list(m.StartupIndicator[g, t, t_p] for t_p in begin_times)
        linear_coefs = [1.]*len(linear_vars)
        if t < initial_time:
            return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 1.)
        else:
            linear_vars.append(m.UnitStop[g,t])
            linear_coefs.append(-1.)
            return (None, linear_expr(linear_coefs=linear_coefs, linear_vars=linear_vars), 0.)
    model.ShutdownMatch = Constraint(model.GeneratorShutdownPeriods, rule=shutdown_match_rule)

    if add_startup_cost_var:
        model.StartupCost = Var(model.SingleFuelGenerators, model.TimePeriods, within=Reals)

    linear_expr = get_linear_expr()
    def ComputeStartupCost2_rule(m,g,t):
        startup_lags = m.ScaledStartupLags[g]
        startup_costs = m.StartupCosts[g]
        last_startup_cost = startup_costs.last()

        linear_vars = [m.StartupCost[g,t], m.UnitStart[g,t]]
        linear_coefs = [-1., last_startup_cost]

        for tp in m.ShutdownsByStartups[g,t]:
            for s in m.StartupCostIndices[g]:
                this_lag = startup_lags[s]
                next_lag = startup_lags[s+1]
                if this_lag <= t - tp < next_lag:
                    linear_vars.append(m.StartupIndicator[g,tp,t])
                    linear_coefs.append(startup_costs[s] - last_startup_cost)
                    break

        return (linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0.)
        '''
        return m.StartupCost[g,t] == m.StartupCosts[g].last()*m.UnitStart[g,t] + \
                                      sum( (m.StartupCosts[g][s] - m.StartupCosts[g].last()) * \
                                         sum( m.StartupIndicator[g,tp,t] for tp in m.ValidShutdownTimePeriods[g] \
                                           if (m.ScaledStartupLags[g][s] <= t - tp < (m.ScaledStartupLags[g][s+1])) ) \
                                         for s in m.StartupCostIndices[g] if s < len(m.StartupCostIndices[g]))
        '''

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

    linear_expr = get_linear_expr(model.UnitStart)

    def delta_eq_rule(m,g,t):
        linear_vars = [m.delta[g,s,t] for s in m.StartupCostIndices[g]]
        linear_coefs = [1.]*len(linear_vars)
        linear_vars.append(m.UnitStart[g,t])
        linear_coefs.append(-1.)
        return (linear_expr( linear_vars, linear_coefs ), 0.)
    
    model.delta_eq=Constraint(model.ThermalGenerators, model.TimePeriods, rule=delta_eq_rule)

    if model.status_vars == 'garver_2bin_vars':
        linear_expr = get_linear_expr(model.UnitOn, model.UnitStart)
    else:
        linear_expr = get_linear_expr(model.UnitStop)
    
    ## BK -- updated to reflect previous generator condition
    ## BK -- assumes initial time is 1
    def delta_ineq_rule(m,g,s,t):
        assert(m.InitialTime == 1)
        ## the last startup indicator doesn't have a lag
        if s == len(m.StartupCostIndices[g]):
            return Constraint.Skip
        this_lag = m.ScaledStartupLags[g][s]
        next_lag = m.ScaledStartupLags[g][s+1]
        ## this is negative if the generator has previously been off
        generator_t0_state = int(round(value(m.UnitOnT0State[g])/value(m.TimePeriodLengthHours)))
        
        ## equation (15) in ME
        if next_lag + generator_t0_state < t < next_lag:
            return (m.delta[g,s,t], 0)

        if m.status_vars == 'garver_2bin_vars':
            linear_vars = [m.UnitStart[g,t-i] for i in range(this_lag, next_lag)]
            linear_coefs = [-1.]*len(linear_vars)
            linear_vars += [m.delta[g,s,t], m.UnitOn[g,t-this_lag]]
            linear_coefs += [1., 1.]
            if t == next_lag:
                return (None, linear_expr(linear_vars, linear_coefs), m.UnitOnT0[g])
                #return m.delta[g,s,t] <= m.UnitOnT0[g] - m.UnitOn[g,t-this_lag] + sum(m.UnitStart[g,t-i] for i in range(this_lag, next_lag))
            elif t > next_lag:
                linear_vars.append(m.UnitOn[g,t-next_lag])
                linear_coefs.apppend(-1.)
                return (None, linear_expr(linear_vars, linear_coefs), 0.)
                #return m.delta[g,s,t] <= m.UnitOn[g,t-next_lag] - m.UnitOn[g,t-this_lag] + sum(m.UnitStart[g,t-i] for i in range(this_lag, next_lag))
        else:
            ## equation (2) in ME
            if t >= next_lag:
                linear_vars = [m.UnitStop[g,t-i] for i in range(this_lag, next_lag)]
                linear_coefs = [-1.]*len(linear_vars)
                linear_vars.append(m.delta[g,s,t])
                linear_coefs.append(1.)
                return (None, linear_expr(linear_vars, linear_coefs), 0.)
    
        return Constraint.Skip
    
    model.delta_ineq=Constraint(model.StartupCostsIndexSet, rule=delta_ineq_rule)

    if add_startup_cost_var:
        model.StartupCost = Var(model.SingleFuelGenerators, model.TimePeriods, within=Reals)
    
    linear_expr = get_linear_expr()
    def ComputeStartupCost2_rule(m,g,t):
        linear_vars = [m.delta[g,s,t] for s in m.StartupCostIndices[g]]
        linear_coefs = [m.StartupCosts[g][s] for s in m.StartupCostIndices[g]]
        linear_vars.append(m.StartupCost[g,t])
        linear_coefs.append(-1.)
        return (linear_expr(linear_vars, linear_coefs),0.)
    
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
        this_lag = m.ScaledStartupLags[g][i]
        this_cost = m.StartupCosts[g][i]
    
        startup_lags = m.ScaledStartupLags[g]
        startup_costs = m.StartupCosts[g]
        
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
                                          for j in range(1, i) )
        else:
            return m.StartupCost[g, t] >= this_cost * m.UnitStart[g, t] \
                                        - sum(
                                               sum(
                                                    (this_cost - startup_costs[j])*m.UnitStop[g, t-k] 
                                               for k in range(startup_lags[j], startup_lags[j+1]) if k < t )
                                          for j in range(1, i) )
    
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
       this_lag = m.ScaledStartupLags[g][i]
       this_cost = m.StartupCosts[g][i]
    
       startup_lags = m.ScaledStartupLags[g]
       startup_costs = m.StartupCosts[g]
       
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
       return m.StartupCost[g, t] >= this_cost * m.UnitOn[g, t] - sum( (this_cost)*m.UnitOn[g, t - k] for k in range(1, min(t, startup_lags[1]+1))) \
                                                                - sum( sum( (this_cost - startup_costs[j])*m.UnitOn[g, t-k] for k in range(startup_lags[j]+1, startup_lags[j+1]+1) if k < t )\
                                                                        for j in range(1, i) )
    
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
       this_lag = m.ScaledStartupLags[g][i]
       this_cost = m.StartupCosts[g][i]
    
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
        this_lag = m.ScaledStartupLags[g][i]
        this_cost = m.StartupCosts[g][i]
    
        startup_lags = m.ScaledStartupLags[g]
        startup_costs = m.StartupCosts[g]
        
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
        this_lag = m.ScaledStartupLags[g][i]
        this_cost = m.StartupCosts[g][i]
    
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
        this_lag = m.ScaledStartupLags[g][i]
        this_cost = m.StartupCosts[g][i]
    
        startup_lags = m.ScaledStartupLags[g]
        startup_costs = m.StartupCosts[g]
        
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
                                          for j in range(1, i) )
        else:
            return m.StartupCostOverHot[g, t] >= (this_cost -m.StartupCosts[g].first())* m.UnitStart[g, t] \
                                        - sum(
                                               sum(
                                                    (this_cost - startup_costs[j])*m.UnitStop[g, t-k] 
                                               for k in range(startup_lags[j], startup_lags[j+1]) if k < t )
                                          for j in range(1, i) )
    
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
        this_lag = m.ScaledStartupLags[g][s]
        next_lag = m.ScaledStartupLags[g][s+1]
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
                sum(m.delta[g,s,t]*(m.StartupCosts[g][s] - m.StartupCosts[g].last()) for s in m.StartupCostIndices[g] if s < len(m.StartupCostIndices[g]))
    model.ComputeStartupCosts=Constraint(model.SingleFuelGenerators, model.TimePeriods, rule=ComputeStartupCost2_rule)

    return

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars','garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            })
def pochet_wolsey_startup_costs(model, add_startup_cost_var=True):
    '''
    Startup costs based on the network-flow formulation
    Pochet & Wosley (2006), Production planning by mixed integer
    programming. Springer Science & Business Media.
    '''

    if add_startup_cost_var:
        model.StartupCost = Var( model.SingleFuelGenerators, model.TimePeriods, within=Reals)

    def ShutdownTimePeriods_generator(m,g):
        ## adds the necessary index for starting-up after a shutdown before the time horizon began
        T0State = value(m.UnitOnT0State[g])
        UT = value(m.ScaledMinimumUpTime[g])
        DT = value(m.ScaledMinimumDownTime[g])
        if T0State > 0:
            first_off = max(0, UT-T0State)
            shutdown_time_periods = list(m.TimePeriods)[first_off:]

        elif T0State < 0:
            first_on = max(0, DT+T0State)
            shutdown_time_periods = [m.InitialTime+T0State] + list(m.TimePeriods)[first_on+UT:]
        return (t for t in shutdown_time_periods)
    model.ShutdownTimePeriods=Set(model.ThermalGenerators, initialize=ShutdownTimePeriods_generator)

    def StartupTimePeriods_generator(m,g):
        ## adds the necessary index for shutting-down after a startup before the time horizon began
        T0State = value(m.UnitOnT0State[g])
        UT = value(m.ScaledMinimumUpTime[g])
        DT = value(m.ScaledMinimumDownTime[g])
        if T0State > 0:
            first_off = max(0, UT-T0State)
            startup_time_periods = [m.InitialTime-T0State]+ list(m.TimePeriods)[first_off+DT:]

        elif T0State < 0:
            first_on = max(0, DT+T0State)
            startup_time_periods = list(m.TimePeriods)[first_on:]
        return (t for t in startup_time_periods) 
    model.StartupTimePeriods=Set(model.ThermalGenerators, initialize=StartupTimePeriods_generator)

    def ShutdownStartupPairs_generator(m,g):
        DT = value(m.ScaledMinimumDownTime[g])
        return ((t_prime, t) for t_prime in m.ShutdownTimePeriods[g] for t in (list(m.TimePeriods)+[m.TimePeriods.last()+DT]) if ( t >= t_prime + DT ))
    model.ShutdownStartupPairs = Set(model.ThermalGenerators, initialize=ShutdownStartupPairs_generator, dimen=2)

    def StartupShutdownPairs_generator(m,g):
        UT = value(m.ScaledMinimumUpTime[g])
        return ((t_prime, t) for t_prime in m.StartupTimePeriods[g] for t in (list(m.TimePeriods)+[m.TimePeriods.last()+UT]) if ( t >= t_prime + UT ))
    model.StartupShutdownPairs = Set(model.ThermalGenerators, initialize=StartupShutdownPairs_generator, dimen=2)

    # (g,t',t) will be an inidicator for g for shutting down at time t' and starting up at time t
    def ShutdownStartupIndicator_index_generator(m):
        return ((g,t_prime,t) for g in m.ThermalGenerators for t_prime,t in m.ShutdownStartupPairs[g]) 
    model.ShutdownStartupIndicator_index=Set(initialize=ShutdownStartupIndicator_index_generator, dimen=3)

    indicator_domain = (UnitInterval if _is_relaxed(model) else Binary)
    model.ShutdownStartupIndicator=Var(model.ShutdownStartupIndicator_index, within=indicator_domain)

    # (g,t',t) will be an inidicator for g for starting up at time t' and shutting down at time t
    def StartupShutdownIndicator_index_generator(m):
        return ((g,t_prime,t) for g in m.ThermalGenerators for t_prime,t in m.StartupShutdownPairs[g]) 
    model.StartupShutdownIndicator_index=Set(initialize=StartupShutdownIndicator_index_generator, dimen=3)

    model.StartupShutdownIndicator=Var(model.StartupShutdownIndicator_index, within=indicator_domain)

    def startup_in_arcs_rule(m, g, t):
        return sum(m.ShutdownStartupIndicator[g, t_prime, s] for (t_prime, s) in m.ShutdownStartupPairs[g] if s == t) == m.UnitStart[g,t]
    model.StartupInArcs = Constraint(model.ThermalGenerators, model.TimePeriods, rule=startup_in_arcs_rule)

    def startup_out_arcs_rule(m, g, t):
        return sum(m.StartupShutdownIndicator[g, s, t_prime] for (s, t_prime) in m.StartupShutdownPairs[g] if s == t) == m.UnitStart[g,t]
    model.StartupOutArcs = Constraint(model.ThermalGenerators, model.TimePeriods, rule=startup_out_arcs_rule)

    def shutdown_in_arcs_rule(m, g, t):
        return sum(m.StartupShutdownIndicator[g, t_prime, s] for (t_prime, s) in m.StartupShutdownPairs[g] if s == t) == m.UnitStop[g,t]
    model.ShutdownInArcs = Constraint(model.ThermalGenerators, model.TimePeriods, rule=shutdown_in_arcs_rule)

    def shutdown_out_arcs_rule(m, g, t):
        return sum(m.ShutdownStartupIndicator[g, s, t_prime] for (s, t_prime) in m.ShutdownStartupPairs[g] if s == t) == m.UnitStop[g,t]
    model.ShutdownOutArcs = Constraint(model.ThermalGenerators, model.TimePeriods, rule=shutdown_out_arcs_rule)

    def unit_on_arcs_rule(m, g, t):
        return sum(m.StartupShutdownIndicator[g, t_prime, s] for (t_prime, s) in m.StartupShutdownPairs[g] if t_prime <= t < s) == m.UnitOn[g,t]
    model.UnitOnArcs = Constraint(model.ThermalGenerators, model.TimePeriods, rule=unit_on_arcs_rule)

    def unit_init_arcs_rule(m, g):
        return sum(m.ShutdownStartupIndicator[g, s, t_prime] for (s, t_prime) in m.ShutdownStartupPairs[g] if s < m.InitialTime) \
             + sum(m.StartupShutdownIndicator[g, s, t_prime] for (s, t_prime) in m.StartupShutdownPairs[g] if s < m.InitialTime) \
             == 1
    model.UnitInitArcs = Constraint(model.ThermalGenerators, rule=unit_init_arcs_rule)

    def unit_fin_arcs_rule(m, g):
        return sum(m.ShutdownStartupIndicator[g, t_prime, s] for (t_prime, s) in m.ShutdownStartupPairs[g] if s > m.TimePeriods.last() ) \
             + sum(m.StartupShutdownIndicator[g, t_prime, s] for (t_prime, s) in m.StartupShutdownPairs[g] if s > m.TimePeriods.last() ) \
             == 1
    model.UnitFinArcs = Constraint(model.ThermalGenerators, rule=unit_fin_arcs_rule)

    def ComputeStartupCost_rule(m,g,t):
        return m.StartupCost[g,t] == m.StartupCosts[g].last()*m.UnitStart[g,t] + \
                                      sum( (m.StartupCosts[g][s] - m.StartupCosts[g].last()) * \
                                         sum( m.ShutdownStartupIndicator[g,tp,t] for tp in m.ShutdownTimePeriods[g] \
                                           if (m.StartupLags[g][s] <= t - tp < (m.StartupLags[g][s+1])) ) \
                                         for s in m.StartupCostIndices[g] if s < len(m.StartupCostIndices[g]))
    model.ComputeStartupCost=Constraint(model.ThermalGenerators, model.TimePeriods, rule=ComputeStartupCost_rule)
