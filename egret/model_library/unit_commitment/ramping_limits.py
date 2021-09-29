#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## file for ramping constraints
from pyomo.environ import *
import math

from .uc_utils import add_model_attr, get_linear_expr
component_name = 'ramping_limits'

## TODO: FIXME: FINISH CONVERTING RAMPING CONSTRAINTS

generation_limits_w_startup_shutdown = ['MLR_generation_limits',
                                        'gentile_generation_limits',
                                        'pan_guan_gentile_generation_limits',
                                        'pan_guan_gentile_KOW_generation_limits',
                                       ]

## NOTE: These two functions control if constraints are even added to the model
##       based on some state data. They will have to be re-visited if we ever 
##       support presistence in a deeper way, just like the minimum up/down time
##       constraints
def _ramp_up_not_needed(m,g,t):
    if m.generation_limits not in generation_limits_w_startup_shutdown:
        return False
    if t == m.InitialTime:
        ## no need for ramping constraints if the generator is off, and
        ## we're enforcing startup/shutdown elsewhere
        if not value(m.UnitOnT0[g]):
            return True
        if value(m.ScaledNominalRampUpLimit[g,t]) >= value(m.MaximumPowerOutput[g,t] - m.PowerGeneratedT0[g]):
            ## the generator can get all the way to max at the first time period
            return True
        return False
    if value(m.ScaledNominalRampUpLimit[g,t]) >= value(m.MaximumPowerOutput[g,t] - m.MinimumPowerOutput[g,t-1]):
        return True
    return False

def _ramp_down_not_needed(m,g,t):
    if t == m.InitialTime:
        if not m.enforce_t1_ramp_rates:
            return True
        ## if the unit is off, we don't need ramp down constraints
        if not value(m.UnitOnT0[g]):
            return True
        if value(m.ScaledNominalRampDownLimit[g,t]) < value(m.PowerGeneratedT0[g] - m.MinimumPowerOutput[g,t]):
            return False
        ## if this and the opposite of the above condition are true, 
        ## we don't need an inital ramp-down inequality
        if value(m.ScaledShutdownRampLimit[g,t]) >= value(m.PowerGeneratedT0[g]):
            return True
        return False
    if m.generation_limits not in generation_limits_w_startup_shutdown:
        return False
    if value(m.ScaledNominalRampDownLimit[g,t]) >= value(m.MaximumPowerOutput[g,t-1] - m.MinimumPowerOutput[g,t]):
        return True
    return False

def _model_time_invariant(m):
    '''
    A test for if certain parameters important to 
    ramping constraints are time invariant
    '''
    ft = m.InitialTime
    for g in m.ThermalGenerators:
        pmax = value(m.MaximumPowerOutput[g,ft])
        pmin = value(m.MinimumPowerOutput[g,ft])
        SU = value(m.ScaledStartupRampLimit[g,ft])
        SD = value(m.ScaledShutdownRampLimit[g,ft])
        RU = value(m.ScaledNominalRampUpLimit[g,ft])
        RD = value(m.ScaledNominalRampDownLimit[g,ft])
        for t in m.TimePeriods:
            if t == ft:
                continue
            if pmax != value(m.MaximumPowerOutput[g,t]) or \
                    pmin != value(m.MinimumPowerOutput[g,t]) or \
                    SU != value(m.ScaledStartupRampLimit[g,t]) or \
                    SD != value(m.ScaledShutdownRampLimit[g,t]) or \
                    RU != value(m.ScaledNominalRampUpLimit[g,t]) or \
                    RD != value(m.ScaledNominalRampDownLimit[g,t]):
                return False
    return True

def _damcikurt_basic_ramping(model):

    linear_expr = get_linear_expr(model.MaximumPowerAvailableAboveMinimum, model.UnitOn, \
                                    model.UnitStart, model.UnitStop)
        
    ## NOTE: with the expression MaximumPowerAvailableAboveMinimum and PowerGeneratedAboveMinimum, 
    ##       these constraints are expressed as needed, there's no cancelation even though we end
    ##       up using these expressions
    def enforce_max_available_ramp_up_rates_rule(m, g, t):
        if _ramp_up_not_needed(m,g,t):
            return Constraint.Skip
        if t == m.InitialTime:
            # if the unit was on in t0, then it's m.PowerGeneratedT0[g] >= m.MinimumPowerOutput[g], and m.UnitOnT0 == 1
            # if not, then m.UnitOnT0[g] == 0 and so (m.PowerGeneratedT0[g] - m.MinimumPowerOutput[g]) * m.UnitOnT0[g] is 0
            ## assume m.MinimumPowerOutput[g,T0] == 0
            linear_vars_power_t, linear_coefs_power_t = m._get_maximum_power_available_above_minimum_lists(m, g, t)
            rhs_linear_vars = [m.UnitOn[g,t], m.UnitStart[g,t]]
            rhs_neg_linear_coefs = [-m.ScaledNominalRampUpLimit[g,t] - 0 + m.MinimumPowerOutput[g,t],
                                    -m.ScaledStartupRampLimit[g,t] + 0 + m.ScaledNominalRampUpLimit[g,t]]

            linear_vars = [*linear_vars_power_t, *rhs_linear_vars]
            linear_coefs = [*linear_coefs_power_t, *rhs_neg_linear_coefs]
            
            RHS = m.PowerGeneratedT0[g]
            return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), RHS)
        else:
            linear_vars_power_t, linear_coefs_power_t = m._get_maximum_power_available_above_minimum_lists(m, g, t)
            linear_vars_power_t_1, linear_coefs_power_t_1 = m._get_negative_power_generated_above_minimum_lists(m, g, t-1)

            rhs_linear_vars = [m.UnitOn[g,t], m.UnitStart[g,t]]
            rhs_neg_linear_coefs = [
                -m.ScaledNominalRampUpLimit[g,t] - m.MinimumPowerOutput[g,t-1] + m.MinimumPowerOutput[g,t],
                -m.ScaledStartupRampLimit[g,t] + m.MinimumPowerOutput[g,t-1] + m.ScaledNominalRampUpLimit[g,t]
                ]

            linear_vars = [*linear_vars_power_t, *linear_vars_power_t_1, *rhs_linear_vars]
            linear_coefs = [*linear_coefs_power_t, *linear_coefs_power_t_1, *rhs_neg_linear_coefs]

            return (None, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0)

    model.EnforceMaxAvailableRampUpRates = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_available_ramp_up_rates_rule)

    def enforce_ramp_down_limits_rule(m, g, t):
        if _ramp_down_not_needed(m,g,t):
            return Constraint.Skip
        if t == m.InitialTime:
            ## assume m.MinimumPowerOutput[g,T0] == 0
            ## TODO: figure out ScaledShutdownRampLimitT0[g]
            linear_vars_power_t, linear_coefs_power_t = m._get_power_generated_above_minimum_lists(m, g, t)

            lhs_linear_vars = [m.UnitStop[g,t]]
            lhs_neg_linear_coefs = [(m.ScaledShutdownRampLimitT0[g] - m.MinimumPowerOutput[g,t] - m.ScaledNominalRampDownLimit[g,t]) ]

            linear_vars = [*linear_vars_power_t, *lhs_linear_vars]
            linear_coefs = [*linear_coefs_power_t, *lhs_neg_linear_coefs]

            LHS = m.PowerGeneratedT0[g] - (m.ScaledNominalRampDownLimit[g,t] + m.MinimumPowerOutput[g,t] - 0)*m.UnitOnT0[g]

            return (LHS, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), None)
        else:
            linear_vars_power_t, linear_coefs_power_t = m._get_power_generated_above_minimum_lists(m, g, t)
            linear_vars_power_t_1, linear_coefs_power_t_1 = m._get_negative_power_generated_above_minimum_lists(m, g, t-1)

            lhs_linear_vars = [m.UnitOn[g,t-1], m.UnitStop[g,t]]
            lhs_neg_linear_coefs = [
                    m.ScaledNominalRampDownLimit[g,t] + m.MinimumPowerOutput[g,t] - m.MinimumPowerOutput[g,t-1],
                    m.ScaledShutdownRampLimit[g,t-1] - m.MinimumPowerOutput[g,t] - m.ScaledNominalRampDownLimit[g,t]
                    ]
            linear_vars = [*linear_vars_power_t, *linear_vars_power_t_1, *lhs_linear_vars]
            linear_coefs = [*linear_coefs_power_t, *linear_coefs_power_t_1, *lhs_neg_linear_coefs]

            return (0, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), None)

    model.EnforceScaledNominalRampDownLimits = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_ramp_down_limits_rule)

    return

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None, 
                                            'reserve_vars': None,
                                            'generation_limits':None,
                                            })
def damcikurt_ramping(model):
    '''
    Equations (3) and (18) from

    Pelin Damci-Kurt, Simge Kucukyavuz, Deepak Rajan, and Alper Atamturk. A
    polyhedral study of production ramping. Mathematical Programming,
    158(1-2):175–205, 2016.
    '''
    _damcikurt_basic_ramping(model)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None, 
                                            'reserve_vars': None,
                                            'generation_limits':None,
                                            })
def damcikurt_ramping_2period(model):
    '''
    Equations (3) and (18), plus equations (20) and (23), from

    Pelin Damci-Kurt, Simge Kucukyavuz, Deepak Rajan, and Alper Atamturk. A
    polyhedral study of production ramping. Mathematical Programming,
    158(1-2):175–205, 2016.
    '''
    _damcikurt_basic_ramping(model)

    if not _model_time_invariant(model):
        raise NotImplementedError("damcikurt_ramping_2period has not be extended to model time-varying minimum or maximum power")

    #TODO: fix for time-varying SU/SD, Pmin, Pmax
    def two_period_ramp_up_rule(m, g, t):
        if value(m.ScaledStartupRampLimit[g,t]) < value(m.MinimumPowerOutput[g,t] + m.ScaledNominalRampUpLimit[g,t]):
            return Constraint.Skip
        j = math.floor(min(value(m.NumTimePeriods)-t, value(m.ScaledStartupRampLimit[g,t] - m.MinimumPowerOutput[g,t])/value(m.ScaledNominalRampUpLimit[g,t])))
        if j > 1: ## j == 1 is handled above
            return m.MaximumPowerAvailableAboveMinimum[g,t+j] - m.PowerGeneratedAboveMinimum[g,t] <= j*m.ScaledNominalRampUpLimit[g,t]*m.UnitOn[g,t+j] \
                    + sum( min(value(m.ScaledStartupRampLimit[g,t] - m.MinimumPowerOutput[g,t] - i*m.ScaledNominalRampUpLimit[g,t]), \
                                value(m.MaximumPowerOutput[g,t] - m.MinimumPowerOutput[g,t] - j*m.ScaledNominalRampUpLimit[g,t]))*m.UnitStart[g,i] for i in range(1, j+1) ) 
        return Constraint.Skip
    model.EnforceTwoPeriodRampUpRule = Constraint(model.ThermalGenerators, model.TimePeriods, rule=two_period_ramp_up_rule)

    def two_period_ramp_down_rule(m, g, t):
        if value(m.ScaledShutdownRampLimit[g,t]) < value(m.MinimumPowerOutput[g,t] + m.ScaledNominalRampDownLimit[g,t]):
            return Constraint.Skip
        j = math.floor(min(value(m.NumTimePeriods)-t, value(m.ScaledShutdownRampLimit[g,t] - m.MinimumPowerOutput[g,t])/value(m.ScaledNominalRampDownLimit[g,t])))
        if j > 1: ## j == 1 is handled above
            return m.PowerGeneratedAboveMinimum[g,t] - m.PowerGeneratedAboveMinimum[g,t+j] <= j*m.ScaledNominalRampDownLimit[g,t]*m.UnitOn[g,t] \
                    + sum( min(value(m.ScaledShutdownRampLimit[g,t] - m.MinimumPowerOutput[g,t] - (j-i+1)*m.ScaledNominalRampDownLimit[g,t]), \
                                value(m.MaximumPowerOutput[g,t] - m.MinimumPowerOutput[g,t] - j*m.ScaledNominalRampDownLimit[g,t]))*m.UnitStop[g,t+i] for i in range(1, j+1) ) 
        return Constraint.Skip
    model.EnforceTwoPeriodRampDownRule = Constraint(model.ThermalGenerators, model.TimePeriods, rule=two_period_ramp_down_rule)


    
@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['ALS_state_transition_vars'],
                                            'power_vars': None, 
                                            'reserve_vars': None,
                                            'generation_limits':None,
                                            })
def ALS_damcikurt_ramping(model):
    '''
    Equations (20a) and (20b) from
    
    Semih Atakan, Guglielmo Lulli, and Suvrajeet Sen. A state transition MIP
    formulation for the unit commitment problem. IEEE Transactions on Power
    Systems, 33(1):736–748, 2018.

    which are modifications of the damcikurt ramping limits.
    '''

    ## TODO: Check math for these constrains under time-varying pmin/pmax
    if not _model_time_invariant(model):
        raise NotImplementedError("ALS_damcikurt_ramping has not be extended to model time-varying minimum or maximum power")


    ## NOTE: with the expression MaximumPowerAvailableAboveMinimum and PowerGeneratedAboveMinimum, 
    ##       these constraints are expressed as needed, there's no cancelation even though we end
    ##       up using these expressions
    def enforce_max_available_ramp_up_rates_rule(m, g, t):
        if _ramp_up_not_needed(m,g,t):
            return Constraint.Skip
        if t == m.InitialTime:
          # if the unit was on in t0, then it's m.PowerGeneratedT0[g] >= m.MinimumPowerOutput[g], and m.UnitOnT0 == 1
          # if not, then m.UnitOnT0[g] == 0 and so (m.PowerGeneratedT0[g] - m.MinimumPowerOutput[g]) * m.UnitOnT0[g] is 0
            return m.MaximumPowerAvailable[g, t] <= m.PowerGeneratedT0[g] - m.MinimumPowerOutput[g,t]*m.UnitOnT0[g] + \
                                                  (m.ScaledNominalRampUpLimit[g,t]+m.MinimumPowerOutput[g,t])*m.UnitStayOn[g,t] + \
    					      m.ScaledStartupRampLimit[g,t]*m.UnitStart[g,t] 
        else:
            return m.MaximumPowerAvailable[g, t] <= m.PowerGeneratedAboveMinimum[g, t-1] + \
                                                  (m.ScaledNominalRampUpLimit[g,t]+m.MinimumPowerOutput[g,t])*m.UnitStayOn[g,t] + \
    					        m.ScaledStartupRampLimit[g,t]*m.UnitStart[g,t] 
    
    model.EnforceMaxAvailableRampUpRates = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_available_ramp_up_rates_rule)
    
    def enforce_ramp_down_limits_rule(m, g, t):
        if _ramp_down_not_needed(m,g,t):
            return Constraint.Skip
        if t == m.InitialTime:
            return m.PowerGeneratedT0[g] - m.MinimumPowerOutput[g,t]*m.UnitOnT0[g] - m.PowerGeneratedAboveMinimum[g, t] <= \
                    m.ScaledNominalRampDownLimit[g,t]*m.UnitStayOn[g,t] + (m.ScaledShutdownRampLimit[g,t] - m.MinimumPowerOutput[g,t])*m.UnitStop[g,t]
        else:
            return m.PowerGeneratedAboveMinimum[g, t-1] - m.PowerGeneratedAboveMinimum[g, t] <= \
                        m.ScaledNominalRampDownLimit[g,t]*m.UnitStayOn[g,t] + (m.ScaledShutdownRampLimit[g,t] - m.MinimumPowerOutput[g,t])*m.UnitStop[g,t]
    
    model.EnforceScaledNominalRampDownLimits = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_ramp_down_limits_rule)

    return

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'power_vars': None,
                                            'reserve_vars': None, 
                                            'generation_limits':generation_limits_w_startup_shutdown,
                                            })
def MLR_ramping(model):
    '''
    Equations (12) and (13) from 

    G. Morales-Espana, J. M. Latorre, and A. Ramos. Tight and compact MILP
    formulation for the thermal unit commitment problem. IEEE Transactions on
    Power Systems, 28(4):4897–4908, 2013.

    with T0 ramp-down limit which is required to make this consistent with other
    formulataions for ramping.

    '''

    # TODO: ADJUST FOR DIFFERING MINIMUM POWER OUTPUTS
    if not _model_time_invariant(model):
        raise NotImplementedError("MLR_ramping has not be extended to model time-varying minimum or maximum power")

    # the following constraint encodes Constraint 12 defined in ME 
    
    def enforce_max_available_ramp_up_rates_rule(m, g, t):
        if _ramp_up_not_needed(m,g,t):
           return Constraint.Skip
        if t == m.InitialTime:
          # if the unit was on in t0, then it's m.PowerGeneratedT0[g] >= m.MinimumPowerOutput[g], and m.UnitOnT0 == 1
          # if not, then m.UnitOnT0[g] == 0 and so (m.PowerGeneratedT0[g] - m.MinimumPowerOutput[g]) * m.UnitOnT0[g] is 0
            return m.MaximumPowerAvailableAboveMinimum[g, t] <= m.PowerGeneratedT0[g] - m.MinimumPowerOutput[g,t]*m.UnitOnT0[g] + \
                                                  m.ScaledNominalRampUpLimit[g,t] 
        else:
            return m.MaximumPowerAvailableAboveMinimum[g, t] <= m.PowerGeneratedAboveMinimum[g, t-1] + m.ScaledNominalRampUpLimit[g,t]
    
    model.EnforceMaxAvailableRampUpRates = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_available_ramp_up_rates_rule)
    
    # the following constraint encodes Constraint 13 defined in ME
    
    def enforce_ramp_down_limits_rule(m, g, t):
        if _ramp_down_not_needed(m,g,t):
            return Constraint.Skip
        if t == m.InitialTime:
            return m.PowerGeneratedT0[g] - m.MinimumPowerOutput[g,t]*m.UnitOnT0[g] - m.PowerGeneratedAboveMinimum[g, t] <= \
                    m.ScaledNominalRampDownLimit[g,t] 
        else:
            return m.PowerGeneratedAboveMinimum[g, t-1] - m.PowerGeneratedAboveMinimum[g, t] <= \
                        m.ScaledNominalRampDownLimit[g,t]
    
    model.EnforceScaledNominalRampDownLimits = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_ramp_down_limits_rule)

    ## need this so we agree with the other ramping models when using MLR Ramping
    ## (i.e., can't shutdown at t=1 unless we're below ScaledShutdownRampLimit)
    def power_limit_t0_stop_rule(m,g):
        if not m.enforce_t1_ramp_rates:
            return Constraint.Skip
        else:
            t = m.InitialTime
            return m.PowerGeneratedT0[g] <= (m.MaximumPowerOutput[g,t])*m.UnitOnT0[g] \
                                            - (m.MaximumPowerOutput[g,t] - m.ScaledShutdownRampLimit[g,t])*m.UnitStop[g,t]
    model.power_limit_t0_stop = Constraint(model.ThermalGenerators,rule=power_limit_t0_stop_rule)

    return


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars','garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            'reserve_vars': None,
                                            'generation_limits':None,
                                            })
def arroyo_conejo_ramping(model):

    '''
    equations (17) and (18) from

    J.M. Arroyo and A.J. Conejo, Optimal Response of a Thermal Unit 
    to an Electricity Spot Market, IEEE Transactions on Power Systems
    Vol. 15, No. 3, Aug 2000
    '''

    # impose upper bounds on the maximum power available for each generator in each time period,
    # based on standard and start-up ramp limits.
    
    # the following constraint encodes Constraint 6 defined in OAV
    
    def enforce_max_available_ramp_up_rates_rule(m, g, t):
        if _ramp_up_not_needed(m,g,t):
            return Constraint.Skip
        # 4 cases, split by (t-1, t) unit status (RHS is defined as the delta from m.PowerGenerated[g, t-1])
        if t == m.InitialTime:
            return m.MaximumPowerAvailable[g, t] <= m.PowerGeneratedT0[g] + \
                                                   m.ScaledNominalRampUpLimit[g,t] * m.UnitOnT0[g] + \
                                                   m.ScaledStartupRampLimit[g,t] * m.UnitStart[g, t]
        else:
            return m.MaximumPowerAvailable[g, t] <= m.PowerGenerated[g, t-1] + \
                                                  m.ScaledNominalRampUpLimit[g,t] * m.UnitOn[g, t-1] + \
                                                  m.ScaledStartupRampLimit[g,t] * m.UnitStart[g,t]
    
    model.EnforceMaxAvailableRampUpRates = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_available_ramp_up_rates_rule)

    # the following constraint encodes Constraint 7 defined in OAV
    
    def enforce_ramp_down_limits_rule(m, g, t):
        if _ramp_down_not_needed(m,g,t):
            return Constraint.Skip
        if t == m.InitialTime:
            return m.PowerGeneratedT0[g] - m.PowerGenerated[g, t] <= \
                 m.ScaledNominalRampDownLimit[g,t] * m.UnitOn[g, t] + \
                 m.ScaledShutdownRampLimitT0[g]  * m.UnitStop[g, t]
        else:
            return m.PowerGenerated[g, t-1] - m.PowerGenerated[g, t] <= \
                 m.ScaledNominalRampDownLimit[g,t]  * m.UnitOn[g, t] + \
                 m.ScaledShutdownRampLimit[g,t-1]  * m.UnitStop[g, t]
    
    model.EnforceScaledNominalRampDownLimits = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_ramp_down_limits_rule)

def _OAV_enhanced(model):
    '''
    baseline for the OAV enhanced formulations
    '''

    ## TODO : TIME-VARYING Pmin/Pmax
    
    # the following constraint encodes Constraint 23 defined in OAV
    
    def enforce_ramp_up_limits_rule(m, g, t):
        if _ramp_up_not_needed(m,g,t):
            return Constraint.Skip
        if (value(m.ScaledNominalRampUpLimit[g,t]) > value(m.ScaledShutdownRampLimit[g,t] - m.MinimumPowerOutput[g,t])) \
                and (value(m.ScaledMinimumUpTime[g]) >= 2):
            if t == m.InitialTime:
                return m.MaximumPowerAvailable[g, t] <= m.PowerGeneratedT0[g] \
                                                       + m.ScaledNominalRampUpLimit[g,t] * m.UnitOn[g,t] \
                                                       + (m.ScaledStartupRampLimit[g,t]-m.ScaledNominalRampUpLimit[g,t]) * m.UnitStart[g, t] \
                                                       - m.MinimumPowerOutput[g,t]*m.UnitStop[g,t] \
                                                       - (m.ScaledNominalRampUpLimit[g,t]- m.ScaledShutdownRampLimit[g,t] + m.MinimumPowerOutput[g,t])*m.UnitStop[g,t+1]
            if t >= value(m.NumTimePeriods):
                return m.MaximumPowerAvailable[g, t] <= m.PowerGenerated[g,t-1] \
                                                       + m.ScaledNominalRampUpLimit[g,t] * m.UnitOn[g,t] \
                                                       + (m.ScaledStartupRampLimit[g,t]-m.ScaledNominalRampUpLimit[g,t]) * m.UnitStart[g, t] \
                                                       - m.MinimumPowerOutput[g,t]*m.UnitStop[g,t]
            else:
                return m.MaximumPowerAvailable[g, t] <= m.PowerGenerated[g,t-1] \
                                                       + m.ScaledNominalRampUpLimit[g,t] * m.UnitOn[g,t] \
                                                       + (m.ScaledStartupRampLimit[g,t]-m.ScaledNominalRampUpLimit[g,t]) * m.UnitStart[g, t] \
                                                       - m.MinimumPowerOutput[g,t]*m.UnitStop[g,t] \
                                                       - (m.ScaledNominalRampUpLimit[g,t]- m.ScaledShutdownRampLimit[g,t] + m.MinimumPowerOutput[g,t])*m.UnitStop[g,t+1]

        else: 
            if t == m.InitialTime:
               return m.MaximumPowerAvailable[g, t] <= m.PowerGeneratedT0[g] + \
                                                       m.ScaledNominalRampUpLimit[g,t] * m.UnitOnT0[g] + \
                                                       m.ScaledStartupRampLimit[g,t] * m.UnitStart[g, t]
            else:
               return m.MaximumPowerAvailable[g, t] <= m.PowerGenerated[g, t-1] + \
                                                       m.ScaledNominalRampUpLimit[g,t] * m.UnitOn[g, t-1] + \
                                                       m.ScaledStartupRampLimit[g,t] * m.UnitStart[g,t]
    
    model.EnforceMaxAvailableRampUpRates = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_ramp_up_limits_rule)

    # the following constraint encodes Constraint 7, 20, 21 defined in OAV
    
    def enforce_ramp_down_limits_rule(m, g, t):
        if _ramp_down_not_needed(m,g,t):
            return Constraint.Skip
        if t == m.InitialTime:
            ## equation 7
            if (value(m.ScaledNominalRampDownLimit[g,t]) <= value(m.ScaledStartupRampLimit[g,t] - m.MinimumPowerOutput[g,t])) \
                    or (value(m.ScaledMinimumUpTime[g]) < 2):
                return m.PowerGeneratedT0[g] - m.PowerGenerated[g, t] <= \
                    m.ScaledNominalRampDownLimit[g,t] * m.UnitOn[g, t] + \
                    m.ScaledShutdownRampLimit[g,t]  * m.UnitStop[g, t]
            elif value(m.ScaledMinimumUpTime[g]) < 3 or value(m.ScaledMinimumDownTime[g]) < 2: # now we can use equation 20
                return m.PowerGeneratedT0[g] - m.PowerGenerated[g, t] <= \
                    + m.ScaledNominalRampDownLimit[g,t] * m.UnitOn[g, t] \
                    + m.ScaledShutdownRampLimit[g,t]  * m.UnitStop[g, t] \
                    - (m.ScaledNominalRampDownLimit[g,t]+m.MinimumPowerOutput[g,t])*m.UnitStart[g,t]
            else: # we can use equation (21)
                return m.PowerGeneratedT0[g] - m.PowerGenerated[g, t] <= \
                    + m.ScaledNominalRampDownLimit[g,t] * m.UnitOn[g, t+1] \
                    + m.ScaledShutdownRampLimit[g,t] * m.UnitStop[g, t] \
                    + m.ScaledNominalRampDownLimit[g,t] * m.UnitStop[g,t+1] \
                    -(m.ScaledNominalRampDownLimit[g,t]+m.MinimumPowerOutput[g,t]) * m.UnitStart[g,t] \
                    - m.ScaledNominalRampDownLimit[g,t] * m.UnitStart[g,t+1]

        else:
            ## equation 7
            if (value(m.ScaledNominalRampDownLimit[g,t]) <= value(m.ScaledStartupRampLimit[g,t] - m.MinimumPowerOutput[g,t])) \
                    or (value(m.ScaledMinimumUpTime[g]) < 2):
                return m.PowerGenerated[g,t-1] - m.PowerGenerated[g, t] <= \
                    m.ScaledNominalRampDownLimit[g,t] * m.UnitOn[g, t] + \
                    m.ScaledShutdownRampLimit[g,t]  * m.UnitStop[g, t]
            elif value(m.ScaledMinimumUpTime[g]) < 3 or value(m.ScaledMinimumDownTime[g]) < 2 or t >= value(m.NumTimePeriods): # now we can use equation 20
                return m.PowerGenerated[g,t-1] - m.PowerGenerated[g, t] <= \
                    + m.ScaledNominalRampDownLimit[g,t] * m.UnitOn[g, t] \
                    + m.ScaledShutdownRampLimit[g,t]  * m.UnitStop[g, t] \
                    -(m.ScaledNominalRampDownLimit[g,t]-m.ScaledStartupRampLimit[g,t]+m.MinimumPowerOutput[g,t])*m.UnitStart[g,t-1] \
                    - (m.ScaledNominalRampDownLimit[g,t]+m.MinimumPowerOutput[g,t])*m.UnitStart[g,t]
            else: # we can use equation (21)
                return m.PowerGenerated[g,t-1] - m.PowerGenerated[g, t] <= \
                    + m.ScaledNominalRampDownLimit[g,t] * m.UnitOn[g, t+1] \
                    + m.ScaledShutdownRampLimit[g,t] * m.UnitStop[g, t] \
                    + m.ScaledNominalRampDownLimit[g,t] * m.UnitStop[g,t+1] \
                    -(m.ScaledNominalRampDownLimit[g,t]-m.ScaledStartupRampLimit[g,t]+m.MinimumPowerOutput[g,t])*m.UnitStart[g,t-1] \
                    -(m.ScaledNominalRampDownLimit[g,t]+m.MinimumPowerOutput[g,t]) * m.UnitStart[g,t] \
                    - m.ScaledNominalRampDownLimit[g,t] * m.UnitStart[g,t+1]
    
    model.EnforceScaledNominalRampDownLimits = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_ramp_down_limits_rule)

## TODO: These should really be refactored so we don't double- or triple-up on ramping limits
@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            'reserve_vars': None,
                                            'generation_limits':None,
                                            })
def OAV_ramping_enhanced(model):

    '''
    Equations (6),(7),(20),(21),(23) from 

    Ostrowski, J., et. al. Tight Mixed Integer Linear Programming Formulations
    for the Unit Commitment Problem. IEEE Transactions on Power Systems, 
    Vol. 27, No. 1, Feb 2012.

    We only add the strongest valid ramp-up or ramp-down equality we can,
    and discard the others
    '''

    ## TODO: check/figure out math for this case
    if not _model_time_invariant(model):
        raise NotImplementedError("OAV_ramping_enhanced has not be extended to model time-varying minimum or maximum power")

    _OAV_enhanced(model)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            'reserve_vars': None,
                                            'generation_limits':None,
                                            })
def OAV_ramping_enhanced_2period(model):

    '''
    This is OAV_ramping_enhanced plus the two-period ramping inequalities
    in equations (22) and (24) from

    Ostrowski, J., et. al. Tight Mixed Integer Linear Programming Formulations
    for the Unit Commitment Problem. IEEE Transactions on Power Systems, 
    Vol. 27, No. 1, Feb 2012.

    '''

    ## TODO: check/figure out math for this case
    if not _model_time_invariant(model):
        raise NotImplementedError("OAV_ramping_enhanced_2period has not be extended to model time-varying minimum or maximum power")

    #TODO: This isn't quite valid, needs debugging
    _OAV_enhanced(model)

    ## TODO: this shouldn't be necessary, and the MaximumPowerAvailable
    ##       should be on the LHS of these equations
    def OAV_two_period_ramp_up_rule(m,g,t):
        if 2*value(m.ScaledNominalRampUpLimit[g,t]) >= value(m.MaximumPowerOutput[g,t] - m.MinimumPowerOutput[g,t]):
            return Constraint.Skip
        if (value(m.ScaledNominalRampUpLimit[g,t]) > value(m.ScaledShutdownRampLimit[g,t] - m.MinimumPowerOutput[g,t])) \
                and (value(m.ScaledMinimumUpTime[g]) >= 2) and (value(m.ScaledMinimumDownTime[g]) >= 2) \
                and (t > value(m.InitialTime)):
            if t == value(m.InitialTime) + 1: ## t == 2
                return m.MaximumPowerAvailable[g,t] - m.PowerGeneratedT0[g] <= \
                          2 * m.ScaledNominalRampUpLimit[g,t] * m.UnitOn[g,t] \
                        - m.MinimumPowerOutput[g,t]*(m.UnitStop[g,t-1]+m.UnitStop[g,t]) \
                        + (m.ScaledStartupRampLimit[g,t] - m.ScaledNominalRampUpLimit[g,t])*m.UnitStart[g,t-1] \
                        + (m.ScaledStartupRampLimit[g,t] - 2*m.ScaledNominalRampUpLimit[g,t])*m.UnitStart[g,t]
            else:
                return m.MaximumPowerAvailable[g,t] - m.PowerGenerated[g,t-2] <= \
                          2 * m.ScaledNominalRampUpLimit[g,t] * m.UnitOn[g,t] \
                        - m.MinimumPowerOutput[g,t]*(m.UnitStop[g,t-1]+m.UnitStop[g,t]) \
                        + (m.ScaledStartupRampLimit[g,t] - m.ScaledNominalRampUpLimit[g,t])*m.UnitStart[g,t-1] \
                        + (m.ScaledStartupRampLimit[g,t] - 2*m.ScaledNominalRampUpLimit[g,t])*m.UnitStart[g,t]
        else:
            return Constraint.Skip
    model.OAVTwoPeriodRampUp = Constraint(model.ThermalGenerators, model.TimePeriods, rule=OAV_two_period_ramp_up_rule)

    ## NOTE: in the text this doesn't have any conditions on when it is valid,
    ##       so the valid conditions were inferred
    def OAV_two_period_ramp_down_rule(m,g,t):
        if 2*value(m.ScaledNominalRampDownLimit[g,t]) >= value(m.MaximumPowerOutput[g,t] - m.MinimumPowerOutput[g,t]):
            return Constraint.Skip
        if (value(m.ScaledNominalRampDownLimit[g,t]) > value(m.ScaledStartupRampLimit[g,t] - m.MinimumPowerOutput[g,t])) \
                and (value(m.ScaledMinimumUpTime[g]) >= 3) and (value(m.ScaledMinimumDownTime[g]) >= 2) \
                and (t > value(m.InitialTime)):
            if t == value(m.InitialTime) + 1: ## t == 2
                if not m.enforce_t1_ramp_rates:
                    return Constraint.Skip
                return m.PowerGeneratedT0[g] - m.PowerGenerated[g,t] <= \
                          2*m.ScaledNominalRampDownLimit[g,t]*m.UnitOn[g,t] \
                        + m.ScaledShutdownRampLimit[g,t]*m.UnitStop[g,t-1] \
                        + (m.ScaledShutdownRampLimit[g,t]+m.ScaledNominalRampDownLimit[g,t])*m.UnitStop[g,t] \
                        -(2*m.ScaledNominalRampDownLimit[g,t]+m.MinimumPowerOutput[g,t])*(m.UnitStart[g,t-1]+m.UnitStart[g,t])
            else:
                return m.PowerGenerated[g,t-2] - m.PowerGenerated[g,t] <= \
                          2*m.ScaledNominalRampDownLimit[g,t]*m.UnitOn[g,t] \
                        + m.ScaledShutdownRampLimit[g,t]*m.UnitStop[g,t-1] \
                        + (m.ScaledShutdownRampLimit[g,t]+m.ScaledNominalRampDownLimit[g,t])*m.UnitStop[g,t] \
                        - 2*m.ScaledNominalRampDownLimit[g,t]*m.UnitStart[g,t-2] \
                        -(2*m.ScaledNominalRampDownLimit[g,t]+m.MinimumPowerOutput[g,t])*(m.UnitStart[g,t-1]+m.UnitStart[g,t])
        else:
            return Constraint.Skip
                
    model.OAVTwoPeriodRampDown = Constraint(model.ThermalGenerators, model.TimePeriods, rule=OAV_two_period_ramp_down_rule)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'power_vars': None,
                                            'reserve_vars': None,
                                            'generation_limits':None,
                                            })
def CA_ramping_limits(model):

    '''
    Equations (18),(19) and (20) from

    Carrion, M. and Arroyo, J. (2006) A Computationally Efficient Mixed-Integer
    Liner Formulation for the Thermal Unit Commitment Problem. IEEE Transactions
    on Power Systems, Vol. 21, No. 3, Aug 2006.
    '''

    # impose upper bounds on the maximum power available for each generator in each time period,
    # based on standard and start-up ramp limits.
    
    # the following constraint encodes Constraint 18 defined in Carrion and Arroyo.
    
    def enforce_max_available_ramp_up_rates_rule(m, g, t):
       # 4 cases, split by (t-1, t) unit status (RHS is defined as the delta from m.PowerGenerated[g, t-1])
       # (0, 0) - unit staying off:   RHS = maximum generator output (degenerate upper bound due to unit being off)
       # (0, 1) - unit switching on:  RHS = startup ramp limit
       # (1, 0) - unit switching off: RHS = standard ramp limit minus startup ramp limit plus maximum power output (degenerate upper bound due to unit off)
       # (1, 1) - unit staying on:    RHS = standard ramp limit plus power generated in previous time period
        if _ramp_up_not_needed(m,g,t):
            return Constraint.Skip
        if t == m.InitialTime:
            return m.MaximumPowerAvailable[g, t] <= m.PowerGeneratedT0[g] + \
                                                   m.ScaledNominalRampUpLimit[g,t] * m.UnitOnT0[g] + \
                                                   m.ScaledStartupRampLimit[g,t] * (m.UnitOn[g, t] - m.UnitOnT0[g]) + \
                                                   m.MaximumPowerOutput[g,t] * (1 - m.UnitOn[g, t])
        else:
            return m.MaximumPowerAvailable[g, t] <= m.PowerGenerated[g, t-1] + \
                                                  m.ScaledNominalRampUpLimit[g,t] * m.UnitOn[g, t-1] + \
                                                  m.ScaledStartupRampLimit[g,t] * (m.UnitOn[g, t] - m.UnitOn[g, t-1]) + \
                                                  m.MaximumPowerOutput[g,t] * (1 - m.UnitOn[g, t])
    
    model.EnforceMaxAvailableRampUpRates = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_available_ramp_up_rates_rule)
    
    
    # the following constraint encodes Constraint 20 defined in Carrion and Arroyo.
    
    def enforce_ramp_down_limits_rule(m, g, t):
        # 4 cases, split by (t-1, t) unit status:
        # (0, 0) - unit staying off:   RHS = maximum generator output (degenerate upper bound)
        # (0, 1) - unit switching on:  RHS = standard ramp-down limit minus shutdown ramp limit plus maximum generator output - this is the strangest case.
        #NOTE: This may never be physically true, but if a generator has ScaledShutdownRampLimit >> MaximumPowerOutput, this constraint causes problems
        # (1, 0) - unit switching off: RHS = shutdown ramp limit
        # (1, 1) - unit staying on:    RHS = standard ramp-down limit
        if _ramp_down_not_needed(m,g,t):
            return Constraint.Skip
        if t == m.InitialTime:
            return m.PowerGeneratedT0[g] - m.PowerGenerated[g, t] <= \
                 m.ScaledNominalRampDownLimit[g,t] * m.UnitOn[g, t] + \
                 m.ScaledShutdownRampLimitT0[g]  * (m.UnitOnT0[g] - m.UnitOn[g, t]) + \
                 m.MaximumPowerOutput[g,t] * (1 - m.UnitOnT0[g])
        else:
            return m.PowerGenerated[g, t-1] - m.PowerGenerated[g, t] <= \
                 m.ScaledNominalRampDownLimit[g,t]  * m.UnitOn[g, t] + \
                 m.ScaledShutdownRampLimit[g,t-1]  * (m.UnitOn[g, t-1] - m.UnitOn[g, t]) + \
                 m.MaximumPowerOutput[g,t] * (1 - m.UnitOn[g, t-1])
    
    model.EnforceScaledNominalRampDownLimits = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_ramp_down_limits_rule)
