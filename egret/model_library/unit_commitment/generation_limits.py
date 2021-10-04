#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## for generation limit constraints 
from pyomo.environ import *
import math

from .uc_utils import add_model_attr, get_linear_expr

component_name = 'generation_limits'

def _add_reactive_limits(model, grid):

    def reactive_upper_limit(m,g,t):
        return m.ReactivePowerGenerated[g,t] <= m.MaximumReactivePowerOutput[g,t]*m.UnitOn[g,t]
    model.EnforceReactiveUpperLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=reactive_upper_limit)

    def reactive_lower_limit(m,g,t):
        return m.MinimumReactivePowerOutput[g,t]*m.UnitOn[g,t] <= m.ReactivePowerGenerated[g,t]
    model.EnforceReactiveLowerLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=reactive_lower_limit)

def _CA_lower_limit(model):

    def enforce_generator_output_limits_rule_part_a(m, g, t):
       return m.MinimumPowerOutput[g,t] * m.UnitOn[g, t] <= m.PowerGenerated[g,t]
    
    model.EnforceGeneratorOutputLimitsPartA = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_generator_output_limits_rule_part_a)

def _get_initial_maximum_power_available_upperbound_lists(m,g,t):
    linear_vars, linear_coefs = m._get_maximum_power_available_lists(m,g,t)
    linear_vars.append(m.UnitOn[g,t])
    linear_coefs.append(-m.MaximumPowerOutput[g,t])
    return linear_vars, linear_coefs

def _get_initial_power_generated_upperbound_lists(m,g,t):
    linear_vars, linear_coefs = m._get_power_generated_lists(m,g,t)
    linear_vars.append(m.UnitOn[g,t])
    linear_coefs.append(-m.MaximumPowerOutput[g,t])
    return linear_vars, linear_coefs

def _MLR_generation_limits_uptime_1(model, tightened=False):

    linear_expr = get_linear_expr(model.UnitOn, model.UnitStart, model.UnitStop)

    ## equations (9), (10) in ME:
    def power_limit_from_start_rule(m,g,t):
        if value(m.ScaledMinimumUpTime[g]) > 1:
            return Constraint.Skip
        linear_vars, linear_coefs = _get_initial_maximum_power_available_upperbound_lists(m,g,t)
        linear_vars.append(m.UnitStart[g,t])
        linear_coefs.append(m.MaximumPowerOutput[g,t] - m.ScaledStartupRampLimit[g,t])
        if t == m.NumTimePeriods or not tightened:
            return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0)
        coef = max(value(m.ScaledStartupRampLimit[g,t]) - value(m.ScaledShutdownRampLimit[g,t]), 0)
        if coef != 0.:
            linear_vars.append(m.UnitStop[g,t+1])
            linear_coefs.append(coef)
        return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0)
    
    model.power_limit_from_start=Constraint(model.ThermalGenerators, model.TimePeriods, rule=power_limit_from_start_rule)
    
    def power_limit_from_stop_rule(m,g,t):
        if value(m.ScaledMinimumUpTime[g]) > 1:
            return Constraint.Skip
        if t == m.NumTimePeriods:
            return Constraint.Skip ## This case is handled above
        linear_vars, linear_coefs = _get_initial_maximum_power_available_upperbound_lists(m,g,t)
        linear_vars.append(m.UnitStop[g,t+1])
        linear_coefs.append(m.MaximumPowerOutput[g,t] - m.ScaledShutdownRampLimit[g,t])
        if not tightened:
            return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0)
        coef = max(value(m.ScaledShutdownRampLimit[g,t]) - value(m.ScaledStartupRampLimit[g,t]), 0)
        if coef != 0.:
            linear_vars.append(m.UnitStart[g,t])
            linear_coefs.append(coef)
        return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0)
    
    model.power_limit_from_stop=Constraint(model.ThermalGenerators, model.TimePeriods, rule=power_limit_from_stop_rule)

def _MLR_generation_limits(model, tightened=False):

    _MLR_generation_limits_uptime_1(model, tightened)

    linear_expr = get_linear_expr(model.UnitOn, model.UnitStart, model.UnitStop)
    
    ## equation (11) in ME:
    def power_limit_from_start_stop_rule(m,g,t):
        if value(m.ScaledMinimumUpTime[g]) <= 1:
            return Constraint.Skip
        linear_vars, linear_coefs = _get_initial_maximum_power_available_upperbound_lists(m,g,t)
        linear_vars.append(m.UnitStart[g,t])
        linear_coefs.append(m.MaximumPowerOutput[g,t] - m.ScaledStartupRampLimit[g,t])
        if t == m.NumTimePeriods: 
            return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0)
        linear_vars.append(m.UnitStop[g,t+1])
        linear_coefs.append(m.MaximumPowerOutput[g,t] - m.ScaledShutdownRampLimit[g,t])
        return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0)
    
    model.power_limit_from_start_stop=Constraint(model.ThermalGenerators,model.TimePeriods,rule=power_limit_from_start_stop_rule)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            'reserve_vars': None,
                                            })
def MLR_generation_limits(model):
    '''
    Equations (9)--(11) from

    G. Morales-Espana, J. M. Latorre, and A. Ramos. Tight and compact MILP
    formulation for the thermal unit commitment problem. IEEE Transactions on
    Power Systems, 28(4):4897–4908, 2013.

    with lower limits if required
    '''
    if model.power_vars in ['garver_power_vars',]:
        pass
    else:
        _CA_lower_limit(model)
    _MLR_generation_limits(model, False)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            'reserve_vars': None,
                                            })
def gentile_generation_limits(model):
    '''
    Equations (1) -- (5) from

    Gentile, G Morales-Espana, and A Ramos. A tight MIP formulation of the unit
    commitment problem with start-up and shut-down constraints. EURO Journal on
    Computational Optimization, 5(1–2):177– 201, 2017.

    with lower limits if required.
    '''
    if model.power_vars in ['garver_power_vars',]:
        pass
    else:
        _CA_lower_limit(model)
    _MLR_generation_limits(model, True)


## generate new time periods by looking forward
def _get_look_forward_periods(m,g,t,UT_end):
    expr = 0.
    p_max_gt = value(m.MaximumPowerOutput[g,t])
    if UT_end is not None:
        end = min(value(UT_end), value(m.NumTimePeriods)-t-1)
    else:
        end = value(m.NumTimePeriods)-t-1
    if end <= 0:
        return end
    ramping_tot = 0
    for i in range(1,end+1):
        shutdown_gi = value(m.ScaledShutdownRampLimit[g,t+i])
        ramping_tot += value(m.ScaledNominalRampDownLimit[g,t+i])
        if shutdown_gi + ramping_tot >= p_max_gt:
            ## the prior index what the correct one
            return i-1
    ## then we can go to the end
    return i

## generate new time periods by looking back
def _get_look_back_periods(m,g,t,UT_end):
    expr = 0.
    p_max_gt = value(m.MaximumPowerOutput[g,t])
    if UT_end is not None:
        end = min(value(UT_end), t-value(m.InitialTime))
    else:
        end = t-value(m.InitialTime)
    ramping_tot = 0
    if end <= 0:
        return end
    for i in range(1,end+1):
        startup_gi = value(m.ScaledStartupRampLimit[g,t-i])
        ramping_tot += value(m.ScaledNominalRampUpLimit[g,t-i])
        if startup_gi + ramping_tot >= p_max_gt:
            ## the prior index what the correct one
            return i-1
    ## then we can go to the end
    return i


def _pan_guan_generation_limits(model, include_UT_1=True):
   
    linear_expr = get_linear_expr(model.UnitOn, model.UnitStart, model.UnitStop)

    #add the stronger ramp-up based inequality, which is a variant of power_limit_from_start_stop
    def power_limit_from_start_stop_rule(m,g,t):
        if (not include_UT_1) and (value(m.ScaledMinimumUpTime[g]) <= 1):
            return Constraint.Skip
        ## time to ramp-up
        Start = m.UnitStart
        Pmax = m.MaximumPowerOutput[g,t]
        SU = m.ScaledStartupRampLimit
        RU = m.ScaledNominalRampUpLimit
        if t == m.NumTimePeriods: 
            ##### ^^^ in this case we can squeeze one more into the sum
            linear_vars, linear_coefs = _get_initial_maximum_power_available_upperbound_lists(m,g,t)
            for i in range(0, _get_look_back_periods(m,g,t,m.ScaledMinimumUpTime[g]-1)+1):
                linear_vars.append(Start[g,t-i])
                linear_coefs.append(Pmax - SU[g,t-i] - sum(RU[g,t-j] for j in range(1,i+1)))
            return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0)
        else:
            linear_vars, linear_coefs = _get_initial_maximum_power_available_upperbound_lists(m,g,t)
            for i in range(0, _get_look_back_periods(m,g,t,m.ScaledMinimumUpTime[g]-2)+1):
                linear_vars.append(Start[g,t-i])
                linear_coefs.append(Pmax - SU[g,t-i] - sum(RU[g,t-j] for j in range(1,i+1)))
            linear_vars.append(m.UnitStop[g,t+1])
            linear_coefs.append(Pmax - m.ScaledShutdownRampLimit[g,t])
            return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0)
    
    model.power_limit_from_start_stop=Constraint(model.ThermalGenerators,model.TimePeriods,rule=power_limit_from_start_stop_rule)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            'reserve_vars': None,
                                            })
def pan_guan_gentile_generation_limits(model):
    '''
    gentile_generation_limits for generators with UT==1 plus
    pan_guan_generation_limits     
    '''
    if model.power_vars in ['garver_power_vars',]:
        pass
    else:
        _CA_lower_limit(model)
    
    #add the strenghtened MLR generation limits fot UT==1 generators
    _MLR_generation_limits_uptime_1(model,True)
    _pan_guan_generation_limits(model,False)


def _KOW_generation_limits(model):

    linear_expr = get_linear_expr(model.UnitOn, model.UnitStart, model.UnitStop)

    ## We'll assume _MLR_generation_limits_uptime_1 and _pan_guan_generation_limits are included
    def max_power_limit_from_starts_rule(m,g,t):
        time_RU = _get_look_back_periods(m,g,t,None)
        if time_RU <= 0:
            return Constraint.Skip
        UT = value(m.ScaledMinimumUpTime[g])
        ## this case is handled better above
        if time_RU <= UT - 2 or t == m.NumTimePeriods:
            return Constraint.Skip
        Start = m.UnitStart
        Pmax = m.MaximumPowerOutput[g,t]
        SU = m.ScaledStartupRampLimit
        RU = m.ScaledNominalRampUpLimit

        linear_vars, linear_coefs = _get_initial_maximum_power_available_upperbound_lists(m,g,t)
        for i in range(0, min(time_RU, UT-1, t-value(m.InitialTime))+1):
            linear_vars.append(Start[g,t-i])
            linear_coefs.append(Pmax - SU[g,t-i] - sum(RU[g,t-j] for j in range(1,i+1)))
        return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0)
        '''
        return m.MaximumPowerAvailable[g,t] <= (m.MaximumPowerOutput[g,t]) *m.UnitOn[g,t] \
                                              - sum((m.MaximumPowerOutput[g,t] - m.ScaledStartupRampLimit[g,t-i] - sum(m.ScaledNominalRampUpLimit[g,t-j] for j in range(1,i+1)))*m.UnitStart[g,t-i] \
                                                                for i in range(0, min(time_RU, UT-1, t-value(m.InitialTime))+1) )
        '''
    model.max_power_limit_from_starts=Constraint(model.ThermalGenerators, model.TimePeriods, rule=max_power_limit_from_starts_rule)

    ## NOTE: it seems this tightening should really be done on the p^l variables, when they exist
    def power_limit_from_start_stops_rule(m,g,t):
        UT = value(m.ScaledMinimumUpTime[g])
        SD_time_limit = _get_look_forward_periods(m,g,t,UT-1)
        if SD_time_limit <= 0: ## this is handled by the _MLR_generation_limits or _pan_guan_generation_limits
                               ## and this is needed so this number isn't negative in the computation of SU_time_limit below
            return Constraint.Skip
        SU_time_limit = _get_look_back_periods(m,g,t,UT-2-SD_time_limit)

        Start = m.UnitStart
        Stop = m.UnitStop
        Pmax = m.MaximumPowerOutput[g,t]

        SU = m.ScaledStartupRampLimit
        SD = m.ScaledShutdownRampLimit
        RU = m.ScaledNominalRampUpLimit
        RD = m.ScaledNominalRampDownLimit

        linear_vars, linear_coefs = _get_initial_power_generated_upperbound_lists(m,g,t)
        for i in range(0, SD_time_limit+1):
            linear_vars.append(Stop[g,t+i+1])
            linear_coefs.append(Pmax - SD[g,t+i] - sum(RD[g,t+j] for j in range(1,i+1)))
        for i in range(0, SU_time_limit+1):
            linear_vars.append(Start[g,t-i])
            linear_coefs.append(Pmax - SU[g,t-i] - sum(RU[g,t-j] for j in range(1,i+1)))

        full_range = (UT >= max(SU_time_limit,0) + max(SD_time_limit,0) + 2)
        if not full_range: 
            i = SU_time_limit+1
            if (t-i) >= value(m.InitialTime):
                coef =  max((value(Pmax) - value(SD[g,t+SD_time_limit]) - sum(value(RD[g,t+j]) for j in range(1,SD_time_limit+1)) - \
                        (value(Pmax) - value(SU[g,t-i]) - sum(value(RU[g,t-j]) for j in range(1,i+1)))), 0)
                if coef != 0:
                    linear_vars.append(Start[g,t-i])
                    linear_coefs.append(coef)
        return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0)

    model.power_limit_from_start_stops=Constraint(model.ThermalGenerators, model.TimePeriods, rule=power_limit_from_start_stops_rule)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            'reserve_vars': None,
                                            })
def pan_guan_gentile_KOW_generation_limits(model):
    #TODO: the equation numbers on this will need updating
    '''
    pan_guan_gentile_generation_limits plus the generalized upper bound
    limits introduced in the text equations (32)--(34)
    
    '''
    
    if model.power_vars in ['garver_power_vars',]:
        pass
    else:
        _CA_lower_limit(model)
    #add the strenghtened MLR generation limits fot UT==1 generators
    _MLR_generation_limits_uptime_1(model,True)
    _pan_guan_generation_limits(model,False)
    _KOW_generation_limits(model)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            'reserve_vars': None,
                                            })
def pan_guan_generation_limits(model):
    '''
    Equation (28) from

    Kai Pan and Yongpei Guan. A polyhedral study of the integrated
    minimum-up/-down time and ramping polytope. arXiv preprint
    arXiv:1604.02184, 2016.

    as modified in the text for different up/down ramps
    '''

    if model.power_vars in ['garver_power_vars',]:
        pass
    else:
        _CA_lower_limit(model)
    _pan_guan_generation_limits(model,True)

def _OAV_generation_limits(model):

    # the following constraint encodes Constraint 19 defined in Carrion and Arroyo, as 
    # modified my MRL to include the UnitStop variable
    
    def enforce_max_available_ramp_down_rates_rule(m, g, t):

        if t == m.NumTimePeriods:
            return m.MaximumPowerAvailable[g, t] <= \
                     m.MaximumPowerOutput[g,t] * m.UnitOn[g, t]
        else:
            return m.MaximumPowerAvailable[g, t] <= \
                     m.MaximumPowerOutput[g,t] * (m.UnitOn[g, t] - m.UnitStop[g,t+1]) + \
                     m.ScaledShutdownRampLimit[g,t] * m.UnitStop[g, t+1]
    
    model.EnforceMaxAvailableRampDownRates = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_available_ramp_down_rates_rule)

# TODO: this needs to be changed to arroyo conejo
@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': ['basic_power_vars', 'garver_power_vars'],
                                            'reserve_vars': None,
                                            })
def OAV_generation_limits(model):
    '''
    Equations (5) from

    Ostrowski, J., et. al. Tight Mixed Integer Linear Programming Formulations
    for the Unit Commitment Problem. IEEE Transactions on Power Systems, 
    Vol. 27, No. 1, Feb 2012.

    Additionally, enforce_max_availabe_ramp_down_rates_rule is from the 
    online companion for MLR (equation (20)), and is included so that reserves do not
    exceed shutdown power in the time period the generator is turned off.
    This makes this formulation then identical to CA and MLR
    '''

    if model.power_vars in ['garver_power_vars',]:
        pass
    else:
        _CA_lower_limit(model)
    _OAV_generation_limits(model)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars','garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': ['basic_power_vars', 'garver_power_vars'],
                                            'reserve_vars': None,
                                            })
def OAV_generation_limits_enhanced(model):
    '''
    Equations (5) and (19) from

    Ostrowski, J., et. al. Tight Mixed Integer Linear Programming Formulations
    for the Unit Commitment Problem. IEEE Transactions on Power Systems, 
    Vol. 27, No. 1, Feb 2012.

    Additionally, enforce_max_availabe_ramp_down_rates_rule is from the 
    online companion for MLR (equation (20)), and is included so that reserves do not
    exceed shutdown power in the time period the generator is turned off.
    This makes this formulation then identical to CA and MLR
    '''

    if model.power_vars in ['garver_power_vars',]:
        pass
    else:
        _CA_lower_limit(model)
    _OAV_generation_limits(model)

    ##TODO: figure out if (19) can be handed to the solver as Lazy
    ##      or as UserCut

    # This is constraint 19 in OAV
    def _K(m,g,t):
        UT = value(m.ScaledMinimumUpTime[g])
        running_ramp = 0
        for k in range(1, UT+1):
            if k+t >= value(m.NumTimePeriods):
                return k-1
            SD = value(m.ScaledShutdownRampLimit[g,t+k])
            oP = value(m.MaximumPowerOutput[g,t])
            if SD + running_ramp >= oP:
                return k-1
            running_ramp += value(m.ScaledNominalRampDownLimit[g,t+k-1])
        return k-1 

    def oav_upper_bound_rule(m,g,t):
        K_t = _K(m,g,t)
        return m.PowerGenerated[g,t] <= m.MaximumPowerOutput[g,t]*m.UnitOn[g,t+K_t] \
                                      + sum((m.ScaledShutdownRampLimit[g,t+i-1] + sum(m.ScaledNominalRampDownLimit[g,t+j-1] for j in range(1,i)))*m.UnitStop[g,t+i] for i in range(1,K_t+1)) \
                                      - sum(m.MaximumPowerOutput[g,t+i]*m.UnitStart[g,t+i] for i in range(1,K_t+1))
    model.OAVUpperBound = Constraint(model.ThermalGenerators, model.TimePeriods, rule=oav_upper_bound_rule)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'power_vars': ['basic_power_vars', 'garver_power_vars'],
                                            'reserve_vars': None,
                                            })
def CA_generation_limits(model):
    '''
    Equations (16) and (17) and (20) from

    Carrion, M. and Arroyo, J. (2006) A Computationally Efficient Mixed-Integer
    Liner Formulation for the Thermal Unit Commitment Problem. IEEE Transactions
    on Power Systems, Vol. 21, No. 3, Aug 2006.

    Note that equation (20) is called a "ramping limit" but it really enforces
    an upper-bound constraint, which is why it is placed here.
    '''

    if model.power_vars in ['garver_power_vars',]:
        pass
    else:
        _CA_lower_limit(model)

    def enforce_max_capacity_rule(m, g, t):
        return m.MaximumPowerAvailable[g,t] <= m.MaximumPowerOutput[g,t]*m.UnitOn[g,t]
    model.EnforceMaxCapacity = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_capacity_rule)

    # the following constraint encodes Constraint 19 defined in Carrion and Arroyo.
    
    def enforce_max_available_ramp_down_rates_rule(m, g, t):
        #NOTE: As expressed in Carrion-Arroyo and subsequently here, this constraint does NOT consider ramp down from initial conditions to t=1!
        if t == value(m.NumTimePeriods):
            return Constraint.Skip
        else:
           return m.MaximumPowerAvailable[g, t] <= \
                  m.MaximumPowerOutput[g,t] * m.UnitOn[g, t+1] + \
                  m.ScaledShutdownRampLimit[g,t] * (m.UnitOn[g, t] - m.UnitOn[g, t+1])
    
    model.EnforceMaxAvailableRampDownRates = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_available_ramp_down_rates_rule)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'power_vars': ['rescaled_power_vars'],
                                            'reserve_vars': ['rescaled_power_avail_vars', ],
                                            })
def YZJMXD_generation_limits(model):
    '''
    Equation (23) from

    Linfeng Yang, Chen Zhang, Jinbao Jian, Ke Meng, Yan Xu, and Zhaoyang Dong.
    A novel projected two-binary-variable formulation for unit commitment in
    power systems. Applied energy, 187:732–745, 2017.
    '''

    def rescaled_generation_limits_rule(m, g, t):
        return m.UnitMaximumPowerAvailableAboveMinimum[g,t] <= m.UnitOn[g,t]
    model.EnforceGenerationLimits = Constraint(model.ThermalGenerators, model.TimePeriods, rule=rescaled_generation_limits_rule)

