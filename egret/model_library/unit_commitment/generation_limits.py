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

from .uc_utils import add_model_attr 

component_name = 'generation_limits'

def _add_reactive_limits(model, grid):

    def reactive_upper_limit(m,g,t):
        return m.ReactivePowerGenerated[g,t] <= m.MaximumReactivePowerOutput[g]*m.UnitOn[g,t]
    model.EnforceReactiveUpperLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=reactive_upper_limit)

    def reactive_lower_limit(m,g,t):
        return m.MinimumReactivePowerOutput[g]*m.UnitOn[g,t] <= m.ReactivePowerGenerated[g,t]
    model.EnforceReactiveLowerLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=reactive_lower_limit)

def _CA_lower_limit(model):

    def enforce_generator_output_limits_rule_part_a(m, g, t):
       return m.MinimumPowerOutput[g] * m.UnitOn[g, t] <= m.PowerGenerated[g,t]
    
    model.EnforceGeneratorOutputLimitsPartA = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_generator_output_limits_rule_part_a)


def _MLR_generation_limits_uptime_1(model, tightened=False):
    ## equations (9), (10) in ME:
    def power_limit_from_start_rule(m,g,t):
        if value(m.ScaledMinimumUpTime[g]) > 1:
            return Constraint.Skip
        if t == m.NumTimePeriods:
            return m.MaximumPowerAvailable[g,t] <= (m.MaximumPowerOutput[g]) * m.UnitOn[g,t] \
                                                - (m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g])*m.UnitStart[g,t]
        return m.MaximumPowerAvailable[g,t] <= (m.MaximumPowerOutput[g]) * m.UnitOn[g,t] \
                                                - (m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g])*m.UnitStart[g,t] \
                        - (max(value(m.ScaledStartupRampLimit[g] - m.ScaledShutdownRampLimit[g]), 0)*m.UnitStop[g,t+1] if tightened else 0.)
    					  
    
    model.power_limit_from_start=Constraint(model.ThermalGenerators, model.TimePeriods, rule=power_limit_from_start_rule)
    
    def power_limit_from_stop_rule(m,g,t):
        if value(m.ScaledMinimumUpTime[g]) > 1:
            return Constraint.Skip
        if t == m.NumTimePeriods:
            return Constraint.Skip ## This case is handled above
        return m.MaximumPowerAvailable[g,t] <= (m.MaximumPowerOutput[g]) * m.UnitOn[g,t] \
                                                - (m.MaximumPowerOutput[g] - m.ScaledShutdownRampLimit[g])*m.UnitStop[g,t+1] \
                        - (max(value(m.ScaledShutdownRampLimit[g] - m.ScaledStartupRampLimit[g]), 0)*m.UnitStart[g,t] if tightened else 0.)
    
    model.power_limit_from_stop=Constraint(model.ThermalGenerators, model.TimePeriods, rule=power_limit_from_stop_rule)

def _MLR_generation_limits(model, tightened=False):

    _MLR_generation_limits_uptime_1(model, tightened)
    
    ## equation (11) in ME:
    def power_limit_from_start_stop_rule(m,g,t):
        if value(m.ScaledMinimumUpTime[g]) <= 1:
            return Constraint.Skip
        if t == m.NumTimePeriods: 
            return m.MaximumPowerAvailable[g,t] <= (m.MaximumPowerOutput[g]) *m.UnitOn[g,t] \
                                                  - (m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g])*m.UnitStart[g,t]
        return m.MaximumPowerAvailable[g,t] <= (m.MaximumPowerOutput[g]) *m.UnitOn[g,t] \
                                              - (m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g])*m.UnitStart[g,t] \
                                              - (m.MaximumPowerOutput[g] - m.ScaledShutdownRampLimit[g])*m.UnitStop[g,t+1]
    
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

def _pan_guan_generation_limits(model, include_UT_1=True):

    #add the stronger ramp-up based inequality, which is a variant of power_limit_from_start_stop
    def power_limit_from_start_stop_rule(m,g,t):
        if (not include_UT_1) and (value(m.ScaledMinimumUpTime[g]) <= 1):
            return Constraint.Skip
        ## time to ramp-up
        if value(m.ScaledNominalRampUpLimit[g]) == 0.:
            ## if the generator can't ramp up, then how can
            ## it ever be above min power?
            return m.MaximumPowerAvailable[g,t] <= m.MinimumPowerOutput[g] * m.UnitOn[g,t]
        time_RU = int(math.floor( value(m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g])/value(m.ScaledNominalRampUpLimit[g]) ))
        if t == m.NumTimePeriods: 
            return m.MaximumPowerAvailable[g,t] <= (m.MaximumPowerOutput[g]) *m.UnitOn[g,t] \
                                                  - sum((m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g] - i*m.ScaledNominalRampUpLimit[g])*m.UnitStart[g,t-i] \
                                                                    for i in range(0, min(time_RU, value(m.ScaledMinimumUpTime[g])-1, t-value(m.InitialTime))+1) )
                                                                                                                            ##### ^^^ in this case we can squeeze one more into the sum
        return m.MaximumPowerAvailable[g,t] <= (m.MaximumPowerOutput[g]) *m.UnitOn[g,t] \
                                              - sum((m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g] - i*m.ScaledNominalRampUpLimit[g])*m.UnitStart[g,t-i] \
                                                                for i in range(0, min(time_RU, value(m.ScaledMinimumUpTime[g])-2, t-value(m.InitialTime))+1) ) \
                                              - (m.MaximumPowerOutput[g] - m.ScaledShutdownRampLimit[g])*m.UnitStop[g,t+1]
    
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

    ## We'll assume _MLR_generation_limits_uptime_1 and _pan_guan_generation_limits are included
    def max_power_limit_from_starts_rule(m,g,t):
        if value(m.ScaledNominalRampUpLimit[g]) == 0.:
            ## if the generator can't ramp up, then how can
            ## it ever be above min power?
            return m.MaximumPowerAvailable[g,t] <= m.MinimumPowerOutput[g] * m.UnitOn[g,t]
        time_RU = int(math.floor( value(m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g])/value(m.ScaledNominalRampUpLimit[g]) ))
        if time_RU <= 0:
            return Constraint.Skip
        UT = value(m.ScaledMinimumUpTime[g])
        ## this case is handled better above
        if time_RU <= UT - 2 or t == m.NumTimePeriods:
            return Constraint.Skip
        return m.MaximumPowerAvailable[g,t] <= (m.MaximumPowerOutput[g]) *m.UnitOn[g,t] \
                                              - sum((m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g] - i*m.ScaledNominalRampUpLimit[g])*m.UnitStart[g,t-i] \
                                                                for i in range(0, min(time_RU, UT-1, t-value(m.InitialTime))+1) )
    model.max_power_limit_from_starts=Constraint(model.ThermalGenerators, model.TimePeriods, rule=max_power_limit_from_starts_rule)

    ## NOTE: it seems this tightening should really be done on the p^l variables, when they exist
    def power_limit_from_start_stops_rule(m,g,t):
        if value(m.ScaledNominalRampUpLimit[g]) == 0.:
            ## if the generator can't ramp up, then how can
            ## it ever be above min power?
            return m.MaximumPowerAvailable[g,t] <= m.MinimumPowerOutput[g] * m.UnitOn[g,t]
        elif value(m.ScaledNominalRampDownLimit[g]) == 0.:
            ## if the generator can't ramp down, then how can
            ## it ever be above min power?
            return m.MaximumPowerAvailable[g,t] <= m.MinimumPowerOutput[g] * m.UnitOn[g,t]
        time_RU = int(math.floor( value(m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g])/value(m.ScaledNominalRampUpLimit[g]) ))
        time_RD = int(math.floor( value(m.MaximumPowerOutput[g] - m.ScaledShutdownRampLimit[g])/value(m.ScaledNominalRampDownLimit[g]) ))
        if time_RU+time_RD <= 0:
            return Constraint.Skip
        UT = value(m.ScaledMinimumUpTime[g])
        SD_time_limit = min(time_RD,value(m.NumTimePeriods)-t-1, UT-1)
        if SD_time_limit <= 0: ## this is handled by the _MLR_generation_limits or _pan_guan_generation_limits
                               ## and this is needed so this number isn't negative in the computation of SU_time_limit below
            return Constraint.Skip
        SU_time_limit = min(time_RU, UT-2-SD_time_limit, t-value(m.InitialTime))
        expr = m.MaximumPowerOutput[g]*m.UnitOn[g,t] \
                                        - sum((m.MaximumPowerOutput[g] - m.ScaledShutdownRampLimit[g] - i*m.ScaledNominalRampDownLimit[g])*m.UnitStop[g,t+i+1] \
                                               for i in range(0, SD_time_limit+1)) \
                                        - sum((m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g] - i*m.ScaledNominalRampUpLimit[g])*m.UnitStart[g,t-i] \
                                               for i in range(0, SU_time_limit+1))                                

        full_range = (UT >= max(SU_time_limit,0) + max(SD_time_limit,0) + 2)
        if not full_range: 
            i = SU_time_limit+1
            if (t-i) >= value(m.InitialTime):
                expr -= max(value( (m.MaximumPowerOutput[g] - m.ScaledShutdownRampLimit[g] - SD_time_limit*m.ScaledNominalRampDownLimit[g]) - (m.MaximumPowerOutput[g] - m.ScaledStartupRampLimit[g] - i*m.ScaledNominalRampUpLimit[g]) ),0)*m.UnitStart[g,t-i]
        return m.PowerGenerated[g,t] <= expr

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
                     m.MaximumPowerOutput[g] * m.UnitOn[g, t]
        else:
            return m.MaximumPowerAvailable[g, t] <= \
                     m.MaximumPowerOutput[g] * (m.UnitOn[g, t] - m.UnitStop[g,t+1]) + \
                     m.ScaledShutdownRampLimit[g] * m.UnitStop[g, t+1]
    
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
        SD = value(m.ScaledShutdownRampLimit[g])
        RD = value(m.ScaledNominalRampDownLimit[g])
        UT = value(m.ScaledMinimumUpTime[g])
        oP = value(m.MaximumPowerOutput[g])
        for k in range(1, UT+1):
            if SD + (k-1)*RD >= oP:
                return k-1
            if k+t >= value(m.NumTimePeriods):
                return k-1
        return k-1 

    def oav_upper_bound_rule(m,g,t):
        K_t = _K(m,g,t)
        return m.PowerGenerated[g,t] <= m.MaximumPowerOutput[g]*m.UnitOn[g,t+K_t] \
                                      + sum((m.ScaledShutdownRampLimit[g] + (i-1)*m.ScaledNominalRampDownLimit[g])*m.UnitStop[g,t+i] for i in range(1,K_t+1)) \
                                      - sum(m.MaximumPowerOutput[g]*m.UnitStart[g,t+i] for i in range(1,K_t+1))
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
        return m.MaximumPowerAvailable[g,t] <= m.MaximumPowerOutput[g]*m.UnitOn[g,t]
    model.EnforceMaxCapacity = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_capacity_rule)

    # the following constraint encodes Constraint 19 defined in Carrion and Arroyo.
    
    def enforce_max_available_ramp_down_rates_rule(m, g, t):
        #NOTE: As expressed in Carrion-Arroyo and subsequently here, this constraint does NOT consider ramp down from initial conditions to t=1!
        if t == value(m.NumTimePeriods):
            return Constraint.Skip
        else:
           return m.MaximumPowerAvailable[g, t] <= \
                  m.MaximumPowerOutput[g] * m.UnitOn[g, t+1] + \
                  m.ScaledShutdownRampLimit[g] * (m.UnitOn[g, t] - m.UnitOn[g, t+1])
    
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

