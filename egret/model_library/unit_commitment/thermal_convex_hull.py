#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## file for large thermal generator convex hull formulations

import pyomo.environ as pe
from .uc_utils import add_model_attr
from .status_vars import _is_relaxed
from .ramping_limits import _ramp_up_not_needed, _ramp_down_not_needed 

## used lots of places in this file
value = pe.value

def _convex_gen_time(m,g,t):
    MinP = value(m.MinimumPowerOutput[g,t])
    if MinP != 0.:
        return False
    MinRunCost = value(m.MinimumProductionCost[g,t])
    if MinRunCost != 0.:
        return False
    UT = value(m.ScaledMinimumUpTime[g])
    if UT > 1:
        return False
    DT = value(m.ScaledMinimumDownTime[g])
    if DT > 1:
        return False

    ## all of the above conditionals must be false,
    ## (their opposite true),
    ## so just look for non-zero startup costs
    StartupCosts = list(m.StartupCosts[g])
    for cost in StartupCosts:
        if value(cost) != 0:
            return False

    ## if all the above conditions are false,
    ## so the generator has MinP == 0, MinRunCost == 0.,
    ## UT, DT <=1, and no start-up costs. Therefore
    ## it is **convex** in its represation at this
    ## time period, and if that is true in every time
    ## period, then it is a generator with a convex
    ## representation
    return True

def get_ramping_gens(m):
    '''
    Returns the ramping-constrained thermal generators
    '''
    ramping_gens = []

    for g in m.ThermalGenerators:
        all_on = all((value(m.FixedCommitment[g,t]) == 1) for t in m.TimePeriods)
        # in this case, we're in the convex hull because the generator is forced on
        if all_on:
            continue
        all_off = all((value(m.FixedCommitment[g,t]) == 0) for t in m.TimePeriods)
        # in this case, we're in the convex hull because the generator is forced off
        if all_off:
            continue

        cached_enforce_t1_ramp_rates = m.enforce_t1_ramp_rates

        ## we don't need these enforced the purposes of this check --
        ## if they're needed the ramping constaint functions would have
        ## already picked these up
        m.enforce_t1_ramp_rates = False
        # in this case, 3-bin is exact
        if all(_ramp_up_not_needed(m,g,t) for t in m.TimePeriods) \
                and all(_ramp_up_not_needed(m,g,t) for t in m.TimePeriods):
            continue
        m.enforce_t1_ramp_rates = cached_enforce_t1_ramp_rates

        # in this case, there are no nonconvex features of this generator
        if all(_convex_gen_time(m,g,t) for t in m.TimePeriods):
            continue

        # if we haven't continued, we can't rule out
        # that the current formulation does not describe
        # ths convex hull
        ramping_gens.append(g)

    return ramping_gens 

def ramping_polytope_block_rule(rp, g):
    '''
    Constructs a block for the ramping polytope as 
    presented in 

    B. Knueven, J. Ostrowski, J. Wang (2018). The Ramping Polytope
    and Cut-Generation for the Unit Commitment Problem. INFORMS 
    Journal on Computing 30(4), 739-749.
    '''

    ## assume the generator parameters
    ## live on the parent block
    model = rp.parent_block()
    value = pe.value
    quicksum = pe.quicksum

    relaxed = _is_relaxed(model)

    UT = value(model.ScaledMinimumUpTime[g])
    DT = value(model.ScaledMinimumDownTime[g])

    timeperiods = model.TimePeriods

    MaxP = { t : model.MaximumPowerOutput[g,t] for t in timeperiods }
    MinP = { t : model.MinimumPowerOutput[g,t] for t in timeperiods }
    RampUp = { t : model.ScaledNominalRampUpLimit[g,t] for t in timeperiods }
    RampDown = { t : model.ScaledNominalRampDownLimit[g,t] for t in timeperiods }
    StartupRamp = { t : model.ScaledStartupRampLimit[g,t] for t in timeperiods }
    ShutdownRamp = { t : model.ScaledShutdownRampLimit[g,t] for t in timeperiods }

    P0 = model.PowerGeneratedT0[g]
    is_on = value(model.UnitOnT0[g])

    pochet_wolsey_startup = model.startup_costs in ['pochet_wolsey_startup_costs']
    if pochet_wolsey_startup: 
        ## in this case, we have some "y" variables to link back to
        def get_y_expr(b, t, t_p):
            return model.StartupShutdownIndicator[g,t,t_p]
        rp.y = pe.Expression(model.StartupShutdownPairs[g], rule=get_y_expr)
        start_stop_pairs = set(model.StartupShutdownPairs[g])
    else:
        if is_on:
            start_stop_pairs = [(bt, et) for bt in timeperiods if (bt >= DT+timeperiods.first()) 
                                    for et in range(bt, timeperiods.last()+1) 
                                    if (bt+UT <= et) ] \
                           + [ (timeperiods.first()-1, t) for t in timeperiods ] \
                           + [ (t, timeperiods.last()+1) for t in timeperiods ] \
                           + [ (timeperiods.first()-1, timeperiods.last()+1) ]
        else:
            start_stop_pairs = [(bt, et) for bt in timeperiods for et in range(bt, timeperiods.last()+1) 
                                    if (bt+UT <= et) ] \
                           + [ (t, timeperiods.last()+1) for t in timeperiods ] \
                           + [ (timeperiods.first()-1, timeperiods.last()+1) ]

        start_stop_pairs = set(start_stop_pairs)
        #print('gen: ', g, ' start_stop_pairs:', sorted(start_stop_pairs))

        ## putting upper bounds on these may cause some dual-degeneracy 
        ## since there are implied bounds based on the UnitOn variables
        ## Not a big deal when in a model, but as a cut subproblem we 
        ## do not want this 
        if relaxed:
            rp.y = pe.Var(start_stop_pairs, within=pe.NonNegativeReals)
        else:
            rp.y = pe.Var(start_stop_pairs, within=pe.NonNegativeIntegers)

        # equivalents to these constraints will have already been added
        rp.downtime_y = pe.Constraint(timeperiods)
        rp.on_link_y = pe.Constraint(timeperiods)
        rp.start_link_y = pe.Constraint(timeperiods)
        rp.stop_link_y = pe.Constraint(timeperiods)
        for t in timeperiods:
            rp.downtime_y[t] = quicksum( rp.y[bt,et] for bt,et in start_stop_pairs if (bt <= t < et+DT) ) <= 1

            rp.on_link_y[t] = quicksum( rp.y[bt,et] for bt,et in start_stop_pairs if (bt <= t < et) ) \
                                 == model.UnitOn[g,t]

            rp.start_link_y[t] = quicksum( rp.y[bt,et] for bt,et in start_stop_pairs if (bt == t) ) \
                                 == model.UnitStart[g,t]

            rp.stop_link_y[t] = quicksum( rp.y[bt,et] for bt,et in start_stop_pairs if (et == t) ) \
                                 == model.UnitStop[g,t]

    rp.p_ints = pe.Var(timeperiods, start_stop_pairs, within=pe.NonNegativeReals)
    rp.r_ints = pe.Var(timeperiods, start_stop_pairs, within=pe.NonNegativeReals)

    rp.ramp_up = pe.Constraint(timeperiods, start_stop_pairs)
    rp.ramp_down = pe.Constraint(timeperiods, start_stop_pairs)
    rp.capacity = pe.Constraint(timeperiods, start_stop_pairs)

    for pair in start_stop_pairs:
        bt,et = pair
        for t in timeperiods:
            if bt > t:
                rp.capacity[t,pair] = rp.p_ints[t,pair] + rp.r_ints[t,pair] <= 0
            elif bt == t:
                rp.capacity[t,pair] = \
                        rp.p_ints[t,pair] + rp.r_ints[t,pair] <= ( StartupRamp[t] - MinP[t])*rp.y[pair]
            elif bt < t:
                if t < et:
                    if t == timeperiods.first():
                        rp.ramp_up[t,pair] = \
                            rp.p_ints[t,pair] + rp.r_ints[t,pair] \
                            - P0*rp.y[pair] \
                                        <= (RampUp[t] + 0  - MinP[t])*rp.y[pair]
                        rp.ramp_down[t,pair] = \
                            P0*rp.y[pair] \
                            - rp.p_ints[t,pair] \
                                        <= (RampDown[t] + MinP[t] - 0)*rp.y[pair]
                    else:
                        rp.ramp_up[t,pair] = \
                            rp.p_ints[t,pair] + rp.r_ints[t,pair] - rp.p_ints[t-1,pair] <= (RampUp[t] + MinP[t-1] - MinP[t])*rp.y[pair]
                        rp.ramp_down[t,pair] = \
                            rp.p_ints[t-1,pair] - rp.p_ints[t,pair] <= (RampDown[t] + MinP[t] - MinP[t-1])*rp.y[pair]
                    if t == et-1 and et <= timeperiods.last():
                        rp.capacity[t,pair] = \
                                rp.p_ints[t,pair] + rp.r_ints[t,pair] <= (ShutdownRamp[t] - MinP[t])*rp.y[pair]
                    else:
                        rp.capacity[t,pair] = \
                                rp.p_ints[t,pair] + rp.r_ints[t,pair] <= (MaxP[t] - MinP[t])*rp.y[pair]
                elif t >= et:
                    rp.capacity[t,pair] = rp.p_ints[t,pair] + rp.r_ints[t,pair] <= 0

    rp.p_link_p_ints = pe.Constraint(timeperiods)
    rp.r_link_r_ints = pe.Constraint(timeperiods)
    for t in timeperiods:
        rp.p_link_p_ints[t] = \
                quicksum( rp.p_ints[t,pair] for pair in start_stop_pairs ) == \
                    model.PowerGeneratedAboveMinimum[g,t]
        rp.r_link_r_ints[t] = \
                quicksum( rp.r_ints[t,pair] for pair in start_stop_pairs ) == \
                    model.ReserveProvided[g,t]

    if len(model.PowerGenerationPiecewisePoints[g,t]) > 2:
        l_range = range(len(model.PowerGenerationPiecewisePoints[g,t])-1)

        l_ints_time = [ (l,t) for l in l_range for t in timeperiods ]
        power_generation_piecewise_points = \
                { (idx,t) : val \
                    for idx,val in enumerate(model.PowerGenerationPiecewisePoints[g,t]) for t in timeperiods }
        rp.pl_ints = pe.Var(l_ints_time, start_stop_pairs, within=pe.NonNegativeReals)
        
        rp.pl_capacity = pe.Constraint(l_ints_time, start_stop_pairs)
        rp.pl_int_p_int_link = pe.Constraint(timeperiods, start_stop_pairs)
        for pair in start_stop_pairs:
            bt,et = pair
            for lt in l_ints_time:
                l,t = lt
                if bt > t:
                    rp.pl_capacity[lt, pair] = rp.pl_ints[lt, pair] <= 0
                else:
                    if t < et:
                        rp.pl_capacity[lt, pair] = \
                            rp.pl_ints[lt, pair] <= \
                            (power_generation_piecewise_points[l+1,t]
                                    - power_generation_piecewise_points[l,t])*rp.y[pair]
                    else:
                        rp.pl_capacity[lt, pair] = rp.pl_ints[lt, pair] <= 0

            for t in timeperiods:
                rp.pl_int_p_int_link[t,pair] = \
                        quicksum( rp.pl_ints[l,t,pair] for l in l_range ) == rp.p_ints[t,pair] 
         
        rp.pl_link = pe.Constraint(l_ints_time)
        for lt in l_ints_time:
            l,t = lt
            rp.pl_link[lt] = \
                    quicksum( rp.pl_ints[lt,pair] for pair in start_stop_pairs ) \
                        == model.PiecewiseProduction[g,t,l]

## we'll require that non-ramping constrained generators have a convex hull description
@add_model_attr('unit_convex_hull', requires = {'data_loader' : None,
                                                'status_vars' : ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                                'power_vars' : None,
                                                'reserve_vars' : None,
                                                'generation_limits' : ['gentile_generation_limits', 'pan_guan_gentile_generation_limits', 'pan_guan_gentile_KOW_generation_limits'],
                                                'production_costs' : ['KOW_production_costs_super_tight', 'KOW_production_costs_tightened', 'rescaled_KOW_production_costs_tightened', ],
                                                'startup_costs' : ['KOW_startup_costs', 'pochet_wolsey_startup_costs'],
                                                'ancillary_service' : None,
                                                })
def add_convex_hull_for_all_units(model):

    ## first, limit the scope
    if len(model.Storage) > 0:
        raise Exception("No convex hull description for storage is implemented")

    if model.nonbasic_reserves:
        raise Exception("No convex hull description for complex ancillary services is implemented")

    if len(model.DualFuelGenerators) > 0 :
        raise Exception("No convex hull description for dual fuel units is implemented")

    ramping_generators = get_ramping_gens(model)
    model.ramping_polytope = pe.Block(ramping_generators, rule=ramping_polytope_block_rule)
