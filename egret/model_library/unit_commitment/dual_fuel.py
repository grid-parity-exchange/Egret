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

    ## Production Costs

    def get_num_piecewise_segments(m, g, p_cost_n):
        return len(dual_fuel_attrs[p_cost_n][g]['values'])-1
    model.NumProductionSegments = Param(model.DualFuelGenerators, ['p_cost', 'aux_p_cost'],
                                        initialize=get_num_piecewise_segments)

    def get_piecewise_generator(p_cost_n):
        p_cost_attr = dual_fuel_attrs[p_cost_n]
        def piecewise_generator(m):
            return ((g,t,i) for g in m.DualFuelGenerators\
                    for t in m.TimePeriods for i in range(m.NumProductionSegments[g,p_cost_n]))
        return piecewise_generator

    model.PriProductionCostIndex = Set(initialize=get_piecewise_generator('p_cost'), dimen=3)
    model.AuxProductionCostIndex = Set(initialize=get_piecewise_generator('aux_p_cost'), dimen=3)
    

    def get_piecewise_bounds_func(p_cost_n):
        p_cost_attr = dual_fuel_attrs[p_cost_n]
        def piecewise_prod_bounds(m, g, t, i):
            return (0, p_cost_attr[g]['values'][i+1][0] - p_cost_attr[g]['values'][i][0])
        return piecewise_prod_bounds

    model.PriPiecewiseProduction = Var( model.PriProductionCostIndex,
                                        within=NonNegativeReals,
                                        bounds=get_piecewise_bounds_func('p_cost'))
    model.AuxPiecewiseProduction = Var( model.AuxProductionCostIndex,
                                        within=NonNegativeReals,
                                        bounds=get_piecewise_bounds_func('aux_p_cost'))
    

    def get_piecewise_constr_rule(p_cost_n):
        p_cost_attr = dual_fuel_attrs[p_cost_n]
        def piecewise_limits_rule(m, g, t, i):
            if p_cost_n == 'p_cost':
                PiecewiseProductionVar = m.PriPiecewiseProduction
                StatusVar = m.UnitOnPriFuel
            elif p_cost_n == 'aux_p_cost':
                PiecewiseProductionVar = m.AuxPiecewiseProduction
                StatusVar = m.UnitOnAuxFuel
            else:
                raise Exception("Unrecognized production cost attribute")
            cost_curve = p_cost_attr[g]['values']
            return PiecewiseProductionVar[g,t,i] <= (cost_curve[i+1][0] - cost_curve[i][0])*StatusVar[g,t]
        return piecewise_limits_rule

    model.PriPiecewiseLimitsConstr = Constraint(model.PriProductionCostIndex, rule=get_piecewise_constr_rule('p_cost'))
    model.AuxPiecewiseLimitsConstr = Constraint(model.AuxProductionCostIndex, rule=get_piecewise_constr_rule('aux_p_cost'))


    def piecewise_production_sum(m, g, t):
        return sum( m.PriPiecewiseProduction[g,t,i] for i in range(m.NumProductionSegments[g,'p_cost'])) \
             + sum( m.AuxPiecewiseProduction[g,t,i] for i in range(m.NumProductionSegments[g,'aux_p_cost'])) \
             == \
               m.PowerGeneratedAboveMinimum[g,t]
    model.DualFuelProductionPiecewiseSum = Constraint(model.DualFuelGenerators, model.TimePeriods, rule=piecewise_production_sum)


    def piecewise_production_cost_rule(m, g, t):
        p_cost_vals = dual_fuel_attrs['p_cost'][g]['values']
        aux_cost_vals = dual_fuel_attrs['aux_p_cost'][g]['values']
        return m.TimePeriodLengthHours * ( \
                sum( (p_cost_vals[i+1][1] - p_cost_vals[i][1]) / (p_cost_vals[i+1][0] - p_cost_vals[i][0])*m.PriPiecewiseProduction[g,t,i] \
                     for i in range(m.NumProductionSegments[g,'p_cost'])) \
              + sum( (aux_cost_vals[i+1][1] - aux_cost_vals[i][1]) / (aux_cost_vals[i+1][0] - aux_cost_vals[i][0])*m.AuxPiecewiseProduction[g,t,i] \
                     for i in range(m.NumProductionSegments[g,'aux_p_cost'])) \
                ) \
             == \
               m.ProductionCost[g,t]
    model.DualFuelProductionCostConstr = Constraint(model.DualFuelGenerators, model.TimePeriods, rule=piecewise_production_cost_rule)

    def no_load_cost_expr_rule(m, g, t):
        p_cost_vals = dual_fuel_attrs['p_cost'][g]['values']
        aux_cost_vals = dual_fuel_attrs['aux_p_cost'][g]['values']
        return m.TimePeriodLengthHours * (p_cost_vals[0][1]*m.UnitOnPriFuel[g,t] + aux_cost_vals[0][1]*m.UnitOnAuxFuel[g,t])
    model.NoLoadCost = Expression(model.DualFuelGenerators, model.TimePeriods, rule=no_load_cost_expr_rule)

    ## End Production Costs

    ## Startup Costs

    ## This is a modified version of the KOW "Matching" startup costs,
    ## with an indicator for each fuel type
    def ValidShutdownTimePeriods_generator(m,g):
        ## for speed, if we don't have different startups
        if len(m.ScaledStartupLags[g]) <= 1 and len(m.AuxScaledStartupLags[g]) <= 1:
            return iter(()) 
                                        ## adds the necessary index for starting-up after a shutdown before the time horizon began
        return (t for t in (list(m.TimePeriods)+([] if (value(m.UnitOnT0State[g]) >= 0) else [m.InitialTime + int(round(value(m.UnitOnT0State[g]/value(m.TimePeriodLengthHours))))])))
    model.ValidShutdownTimePeriodsD =Set(model.DualFuelGenerators, initialize=ValidShutdownTimePeriods_generator)

    def GeneratorShutdownPeriods_generator(m):
        return ((g,t) for g in m.DualFuelGenerators for t in m.ValidShutdownTimePeriodsD[g])
    model.GeneratorShutdownPeriodsD = Set(initialize=GeneratorShutdownPeriods_generator, dimen=2)
    
    def ShutdownHotStartupPairs_generator(m,g):
        ## for speed, if we don't have different startups
        if len(m.ScaledStartupLags[g]) <= 1:
            return [] 
        return ((t_prime, t) for t_prime in m.ValidShutdownTimePeriodsD[g] for t in m.TimePeriods if (m.ScaledStartupLags[g].first() <= t - t_prime < m.ScaledStartupLags[g].last()))
    model.PriShutdownHotStartupPairs = Set(model.DualFuelGenerators, initialize=ShutdownHotStartupPairs_generator, dimen=2)

    def AuxShutdownHotStartupPairs_generator(m,g):
        ## for speed, if we don't have different startups
        if len(m.ScaledStartupLags[g]) <= 1:
            return [] 
        return ((t_prime, t) for t_prime in m.ValidShutdownTimePeriodsD[g] for t in m.TimePeriods if (m.AuxScaledStartupLags[g].first() <= t - t_prime < m.AuxScaledStartupLags[g].last()))
    model.AuxShutdownHotStartupPairs = Set(model.DualFuelGenerators, initialize=AuxShutdownHotStartupPairs_generator, dimen=2)
    
    # (g,t',t) will be an inidicator for g for shutting down at time t' and starting up at time t
    def PriStartupIndicator_domain_generator(m):
        return ((g,t_prime,t) for g in m.DualFuelGenerators for t_prime,t in m.PriShutdownHotStartupPairs[g]) 
    model.PriStartupIndicator_domain=Set(initialize=PriStartupIndicator_domain_generator, dimen=3)

    # (g,t',t) will be an inidicator for g for shutting down at time t' and starting up at time t
    def AuxStartupIndicator_domain_generator(m):
        return ((g,t_prime,t) for g in m.DualFuelGenerators for t_prime,t in m.AuxShutdownHotStartupPairs[g]) 
    model.AuxStartupIndicator_domain=Set(initialize=AuxStartupIndicator_domain_generator, dimen=3)
    
    def binary_or_unit_interval_startup(m, g, t_p, t):
        return binary_or_unit_interval(m, g, t)

    model.PriStartupIndicator=Var(model.PriStartupIndicator_domain, within=binary_or_unit_interval_startup)
    model.AuxStartupIndicator=Var(model.AuxStartupIndicator_domain, within=binary_or_unit_interval_startup)

    ############################################################
    # compute the per-generator, per-time period startup costs #
    ############################################################
    
    def Pri_startup_match_rule(m, g, t):
        return sum(m.PriStartupIndicator[g, t_prime, s] for (t_prime, s) in m.PriShutdownHotStartupPairs[g] if s == t) \
                <= m.UnitStartPriFuel[g,t]
    model.PriStartupMatch = Constraint(model.DualFuelGenerators, model.TimePeriods, rule=Pri_startup_match_rule)

    def Aux_startup_match_rule(m, g, t):
        return sum(m.AuxStartupIndicator[g, t_prime, s] for (t_prime, s) in m.AuxShutdownHotStartupPairs[g] if s == t) \
                <= m.UnitStartAuxFuel[g,t]
    model.AuxStartupMatch = Constraint(model.DualFuelGenerators, model.TimePeriods, rule=Aux_startup_match_rule)
    
    def shutdown_match_rule(m, g, t):
        if t < m.InitialTime:
            p_begin_pairs = [(s, t_prime) for (s, t_prime) in m.PriShutdownHotStartupPairs[g] if s == t]
            a_begin_pairs = [(s, t_prime) for (s, t_prime) in m.AuxShutdownHotStartupPairs[g] if s == t] 
            if (not p_begin_pairs) and (not a_begin_pairs): ##if both are empty
                return Constraint.Feasible
            else:
                return sum(m.PriStartupIndicator[g, s, t_prime] for (s, t_prime) in p_begin_pairs) \
                        + sum(m.AuxStartupIndicator[g, s, t_prime] for (s, t_prime) in a_begin_pairs) <= 1
        else:
            return  sum(m.PriStartupIndicator[g, s, t_prime] for (s, t_prime) in m.PriShutdownHotStartupPairs[g] if s == t) \
                  + sum(m.AuxStartupIndicator[g, s, t_prime] for (s, t_prime) in m.AuxShutdownHotStartupPairs[g] if s == t) \
                  <= m.UnitStop[g,t]
    model.ShutdownMatchD = Constraint(model.GeneratorShutdownPeriodsD, rule=shutdown_match_rule)

    def ComputeStartupCost2_rule(m,g,t):
        return m.StartupCost[g,t] == m.StartupCosts[g].last()*m.UnitStartPriFuel[g,t] + \
                                      sum( (list(m.StartupCosts[g])[s-1] - m.StartupCosts[g].last()) * \
                                         sum( m.PriStartupIndicator[g,tp,t] for tp in m.ValidShutdownTimePeriodsD[g] \
                                           if (list(m.ScaledStartupLags[g])[s-1] <= t - tp < (list(m.ScaledStartupLags[g])[s])) ) \
                                         for s in m.StartupCostIndices[g] if s < len(m.StartupCostIndices[g])) \
                                   +  m.AuxStartupCosts[g].last()*m.UnitStartAuxFuel[g,t] + \
                                      sum( (list(m.AuxStartupCosts[g])[s-1] - m.AuxStartupCosts[g].last()) * \
                                         sum( m.AuxStartupIndicator[g,tp,t] for tp in m.ValidShutdownTimePeriodsD[g] \
                                           if (list(m.AuxScaledStartupLags[g])[s-1] <= t - tp < (list(m.AuxScaledStartupLags[g])[s])) ) \
                                         for s in m.AuxStartupCostIndices[g] if s < len(m.AuxStartupCostIndices[g])) 
    model.ComputeStartupCostsD=Constraint(model.DualFuelGenerators, model.TimePeriods, rule=ComputeStartupCost2_rule)
