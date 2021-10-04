#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## file for production cost functions
from pyomo.environ import *
import math
from functools import lru_cache

from .uc_utils import add_model_attr, get_linear_expr
from .generation_limits import _get_look_back_periods, _get_look_forward_periods

component_name = 'production_costs'


## NOTE: for now we'll just consider all piecewise variables to represent
##       power above minimum. This is how it's done the the Carrion-Arroyo 
##       paper, as well as several others

# a function for use in piecewise linearization of the cost function.
@lru_cache()
def _production_cost_function(m, g, t, i):
    return m.TimePeriodLengthHours * m.PowerGenerationPiecewiseCostValues[g,t][i]

def _compute_total_production_cost(model):

    ## helper function for PH
    def compute_production_costs_rule(m, g, t, avg_power):
        ## piecewise points for power
        piecewise_points = m.PowerGenerationPiecewisePoints[g,t]
        ## buckets
        piecewise_eval = [0]*(len(piecewise_points)-1)
        ## fill the buckets (skip the first since it's min power)
        for l in range(len(piecewise_eval)):
            ## fill this bucket all the way
            if avg_power >= piecewise_points[l+1]:
                piecewise_eval[l] = piecewise_points[l+1] - piecewise_points[l]
            ## fill the bucket part way and stop
            elif avg_power < piecewise_points[l+1]:
                piecewise_eval[l] = avg_power - piecewise_points[l]
                break
    
                #slope * production
        return sum( (_production_cost_function(m,g,t,l+1) - _production_cost_function(m,g,t,l)) / (piecewise_points[l+1] - piecewise_points[l]) * piecewise_eval[l] for l in range(len(piecewise_eval))) 
    
    model.ComputeProductionCosts = compute_production_costs_rule

def _get_piecewise_production_generators(model):

    # more than two points -> not linear
    def piecewise_generators_time_set(m):
        for g in m.ThermalGenerators:
            for t in m.TimePeriods:
                if len(m.PowerGenerationPiecewisePoints[g,t]) > 2:
                    yield g,t
    model.PiecewiseGeneratorTimeIndexSet = Set(dimen=2, initialize=piecewise_generators_time_set)

    # two points -> linear
    def linear_generators_time_set(m):
        for g in m.ThermalGenerators:
            for t in m.TimePeriods:
                if len(m.PowerGenerationPiecewisePoints[g,t]) == 2:
                    yield g,t
    model.LinearGeneratorTimeIndexSet = Set(dimen=2, initialize=linear_generators_time_set)

    # if there's only 1 or zero points, this has no marginal cost

    # compute the per-generator, per-time period production costs. We'll do this by hand
    def piecewise_production_costs_index_set_generator(m):
        return ((g,t,i) for g,t in m.PiecewiseGeneratorTimeIndexSet for i in range(len(m.PowerGenerationPiecewisePoints[g,t])-1))
    model.PiecewiseProductionCostsIndexSet = Set(initialize=piecewise_production_costs_index_set_generator, dimen=3)


def _basic_production_costs_vars(model):

    _get_piecewise_production_generators(model)

    def piecewise_production_bounds_rule(m, g, t, i):
        return (0, m.PowerGenerationPiecewisePoints[g,t][i+1] - m.PowerGenerationPiecewisePoints[g,t][i])
    
    model.PiecewiseProduction = Var( model.PiecewiseProductionCostsIndexSet, within=NonNegativeReals, bounds = piecewise_production_bounds_rule )
    
    linear_expr = get_linear_expr(model.PowerGeneratedAboveMinimum)

    def piecewise_production_sum_rule(m, g, t):
        linear_vars = list( m.PiecewiseProduction[g,t,i] for i in range(len(m.PowerGenerationPiecewisePoints[g,t])-1))
        linear_coefs = [1.]*len(linear_vars)
        linear_vars.append(m.PowerGeneratedAboveMinimum[g,t])
        linear_coefs.append(-1.)
        return (linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0.)
    model.PiecewiseProductionSum = Constraint( model.PiecewiseGeneratorTimeIndexSet, rule=piecewise_production_sum_rule )

def _basic_production_costs_constr(model):

    model.ProductionCost = Var( model.SingleFuelGenerators, model.TimePeriods, within=Reals )

    linear_expr = get_linear_expr()

    def piecewise_production_costs_rule(m, g, t):
        
        if (g,t) in m.PiecewiseGeneratorTimeIndexSet:
            points = m.PowerGenerationPiecewisePoints[g,t]
            costs = m.PowerGenerationPiecewiseCostValues[g,t]
            time_scale = m.TimePeriodLengthHours
            linear_coefs = [(time_scale*costs[i+1] - time_scale*costs[i])/(points[i+1] - points[i]) \
                        for i in range(len(points)-1)]
            linear_vars = [m.PiecewiseProduction[g,t,i] for i in range(len(points)-1)]
            linear_coefs.append(-1.)
            linear_vars.append(m.ProductionCost[g,t])
            return (linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0.)
        elif (g,t) in m.LinearGeneratorTimeIndexSet:
            i = 0
            points = m.PowerGenerationPiecewisePoints[g,t]
            costs = m.PowerGenerationPiecewiseCostValues[g,t]
            time_scale = m.TimePeriodLengthHours
            slope =(time_scale*costs[i+1] - time_scale*costs[i])/(points[i+1] - points[i])
            linear_vars, linear_coefs = m._get_power_generated_above_minimum_lists(m,g,t)
            linear_coefs = [slope*coef for coef in linear_coefs]
            linear_vars.append(m.ProductionCost[g,t])
            linear_coefs.append(-1.)
            return (linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), 0.)
        else:
            return (m.ProductionCost[g,t], 0.)
    
    model.ProductionCostConstr = Constraint( model.SingleFuelGenerators, model.TimePeriods, rule=piecewise_production_costs_rule )

    _compute_total_production_cost(model)


# BEGIN production cost models
# BEGIN proudction cost models which use auxilarily variables for each
#       piecewise segment
@add_model_attr(component_name, requires = {'data_loader': None,
                                            'power_vars': None,
                                            })
def CA_production_costs(model):

    '''
    This is a production cost model with additional variables for each
    piecewise linear segment, equations (7)--(11) in

    M. Carrion and J. M. Arroyo. A computationally efficient mixed-integer
    linear formulation for the thermal unit commitment problem. IEEE
    Transactions on Power Systems, 21(3):1371–1378, Aug 2006. ISSN 0885-8950.
    doi: 10.1109/TPWRS.2006.876672.
    ''' 
    
    _basic_production_costs_vars(model)
    _basic_production_costs_constr(model)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'power_vars': None,
                                            })
def wu_production_costs(model):

    '''
    This is the (ideal) production cost model introducted by:
    
    Wu, Lei, 2016. Accelerating NCUC Via Binary Variable-Based
    Locally Ideal Formulation and Dynamic Global Cuts,
    IEEE Trans. on Power Sys., Vol 31, No 5.

    equation (14)
    ''' 
    
    _basic_production_costs_vars(model)

    def piecewise_production_limits_rule(m, g, t, i):
        return m.PiecewiseProduction[g,t,i] <= (m.PowerGenerationPiecewisePoints[g,t][i+1] - m.PowerGenerationPiecewisePoints[g,t][i])*m.UnitOn[g,t]
    model.PiecewiseProductionLimits = Constraint( model.PiecewiseProductionCostsIndexSet, rule=piecewise_production_limits_rule )
    
    _basic_production_costs_constr(model) 


# function for helping lay down piecewise curves right
@lru_cache()
def _step_coeff(upper, lower, susd):
    if lower < susd < upper:
        return upper - susd
    elif susd <= lower:
        return upper - lower
    elif upper <= susd:
        return 0
    print("Something went wrong, step_coeff is returning None...")
    print("lower: ", lower)
    print("susd: ", susd)
    print("upper: ", upper)
    return None

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_3bin_relaxed_stop_vars', 'garver_2bin_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            })
def KOW_production_costs_super_tight(model):
    '''
    production costs which take into account the ramping trajectories
    as noted, but not formulated, in text
    '''
    _basic_production_costs_vars(model)

    def piecewise_production_limits_from_start_stops_rule(m, g, t, i):
        RU = m.ScaledNominalRampUpLimit
        RD = m.ScaledNominalRampDownLimit
        if value(RU[g,t]) == 0. or value(RD[g,t]) == 0.:
            ## again, if we can't ramp, there will be no power above minimum
            return m.PowerGenerationPiecewisePoints[g,t][i] <= 0.
        SU = m.ScaledStartupRampLimit
        SD = m.ScaledShutdownRampLimit

        UT = value(m.ScaledMinimumUpTime[g])
        SU_time_limit = _get_look_back_periods(m,g,t,UT-1)
        SD_time_limit = _get_look_forward_periods(m,g,t,UT-2-SU_time_limit)
        full_range = (UT >= max(SU_time_limit,0) + max(SD_time_limit,0) + 2)
        ### these can always be tightened based on SU/SD, regardless of the ramping/aggregation
        ### since PowerGenerationPiecewisePoints are scaled to MinimumPowerOutput, we need to scale Startup/Shutdown ramps to it as well
        upper = value(m.PowerGenerationPiecewisePoints[g,t][i+1])
        lower = value(m.PowerGenerationPiecewisePoints[g,t][i])
        minP = m.MinimumPowerOutput
        su_step = {}
        sd_step = {}
        for j in range(0,SU_time_limit+1):
            su_step[j] = _step_coeff(upper, lower, value(SU[g,t-j]-minP[g,t-j]) + value(sum(RU[g,t-i] for i in range(1,j+1))))
        ## grab one more if we don't cover the full range
        for j in range(0,SD_time_limit+1+(0 if full_range else 1)):
            sd_step[j] = _step_coeff(upper, lower, value(SD[g,t+j]-minP[g,t+j])+value(sum(RD[g,t+i] for i in range(1,j+1))))

        expr = (m.PowerGenerationPiecewisePoints[g,t][i+1] - m.PowerGenerationPiecewisePoints[g,t][i])*m.UnitOn[g,t] \
                                        - sum(su_step[j]*m.UnitStart[g,t-j] for j in range(0,SU_time_limit+1)) \
                                        - sum(sd_step[j]*m.UnitStop[g,t+1+j] for j in range(0,SD_time_limit+1))

        if not full_range: 
            j = SD_time_limit+1
            if t+1+j <= value(m.NumTimePeriods):
                expr -= max(sd_step[j]-su_step[SU_time_limit],0)*m.UnitStop[g,t+1+j]
        return m.PiecewiseProduction[g,t,i] <= expr      
    model.PiecewiseProductionLimits = Constraint( model.PiecewiseProductionCostsIndexSet, rule=piecewise_production_limits_from_start_stops_rule )
    
    def piecewise_production_limits_from_stops_start_rule(m, g, t, i):
        ### these can always be tightened based on SU/SD, regardless of the ramping/aggregation
        ### since PowerGenerationPiecewisePoints are scaled to MinimumPowerOutput, we need to scale Startup/Shutdown ramps to it as well
        RU = m.ScaledNominalRampUpLimit
        RD = m.ScaledNominalRampDownLimit
        if value(RU[g,t]) == 0. or value(RD[g,t]) == 0.:
            ## again, if we can't ramp, there will be no power above minimum
            return m.PowerGenerationPiecewisePoints[g,t][i] <= 0.
        SU = m.ScaledStartupRampLimit
        SD = m.ScaledShutdownRampLimit

        UT = value(m.ScaledMinimumUpTime[g])
        SD_time_limit = _get_look_forward_periods(m,g,t,UT-1)
        SU_time_limit = _get_look_back_periods(m,g,t,UT-2-SD_time_limit)

        # in this case, this will give the same constraint as the rule above
        if (UT >= max(SU_time_limit,0) + max(SD_time_limit,0) + 2) or (SD_time_limit < 0):
            return Constraint.Skip
        ### these can always be tightened based on SU/SD, regardless of the ramping/aggregation
        ### since PowerGenerationPiecewisePoints are scaled to MinimumPowerOutput, we need to scale Startup/Shutdown ramps to it as well
        upper = value(m.PowerGenerationPiecewisePoints[g,t][i+1])
        lower = value(m.PowerGenerationPiecewisePoints[g,t][i])
        minP = m.MinimumPowerOutput
        su_step = {}
        sd_step = {}
        for j in range(0,SU_time_limit+2):
            su_step[j] = _step_coeff(upper, lower, value(SU[g,t-j]-minP[g,t-j]) + value(sum(RU[g,t-i] for i in range(1,j+1))))
        for j in range(0,SD_time_limit+1):
            sd_step[j] = _step_coeff(upper, lower, value(SD[g,t+j]-minP[g,t+j])+value(sum(RD[g,t+i] for i in range(1,j+1))))
        expr = (m.PowerGenerationPiecewisePoints[g,t][i+1] - m.PowerGenerationPiecewisePoints[g,t][i])*m.UnitOn[g,t] \
                                        - sum(su_step[j]*m.UnitStart[g,t-j] for j in range(0,SU_time_limit+1)) \
                                        - sum(sd_step[j]*m.UnitStop[g,t+1+j] for j in range(0,SD_time_limit+1))

        j = SU_time_limit+1
        if (t-j) >= value(m.InitialTime):
            expr -= max(su_step[j]-sd_step[SD_time_limit],0)*m.UnitStart[g,t-j]
        return m.PiecewiseProduction[g,t,i] <= expr
    model.PiecewiseProductionLimits2 = Constraint( model.PiecewiseProductionCostsIndexSet, rule=piecewise_production_limits_from_stops_start_rule )

    _basic_production_costs_constr(model)


def _KOW_production_costs(model, tightened = False):
    '''
    Base for similarities between tightend and not KOW production costs
    '''
    _basic_production_costs_vars(model)

    linear_expr = get_linear_expr(model.UnitOn, model.UnitStart, model.UnitStop)

    def piecewise_production_limits_rule(m, g, t, i):
        ### these can always be tightened based on SU/SD, regardless of the ramping/aggregation
        ### since PowerGenerationPiecewisePoints are scaled to MinimumPowerOutput, we need to scale Startup/Shutdown ramps to it as well
        upper = value(m.PowerGenerationPiecewisePoints[g,t][i+1])
        lower = value(m.PowerGenerationPiecewisePoints[g,t][i])
        SU = value(m.ScaledStartupRampLimit[g,t])
        minP = value(m.MinimumPowerOutput[g,t])

        su_step = _step_coeff(upper, lower, SU-minP)
        if t < value(m.NumTimePeriods):
            SD = value(m.ScaledShutdownRampLimit[g,t])
            UT = value(m.ScaledMinimumUpTime[g])
            sd_step = _step_coeff(upper, lower, SD-minP)
            if UT > 1:
                linear_vars = [m.PiecewiseProduction[g,t,i], m.UnitOn[g,t], m.UnitStart[g,t], m.UnitStop[g,t+1]]
                linear_coefs = [-1., (upper-lower), -su_step, -sd_step]
                return (0, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), None)
            else: ## MinimumUpTime[g] <= 1
                linear_vars = [m.PiecewiseProduction[g,t,i], m.UnitOn[g,t], m.UnitStart[g,t],]
                linear_coefs = [-1., (upper-lower), -su_step,]
                if tightened:
                    coef = -max(sd_step-su_step,0)
                    if coef != 0:
                        linear_vars.append(m.UnitStop[g,t+1])
                        linear_coefs.append(coef)
                return (0, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), None)

        else: ## t >= value(m.NumTimePeriods)
            linear_vars = [m.PiecewiseProduction[g,t,i], m.UnitOn[g,t], m.UnitStart[g,t],]
            linear_coefs = [-1., (upper-lower), -su_step,]
            return (0, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), None)

    model.PiecewiseProductionLimits = Constraint( model.PiecewiseProductionCostsIndexSet, rule=piecewise_production_limits_rule )
    
    def piecewise_production_limits_rule2(m, g, t, i):
        ### these can always be tightened based on SU/SD, regardless of the ramping/aggregation
        ### since PowerGenerationPiecewisePoints are scaled to MinimumPowerOutput, we need to scale Startup/Shutdown ramps to it as well
        UT = value(m.ScaledMinimumUpTime[g])
        if UT <= 1 and t < value(m.NumTimePeriods):
            upper = value(m.PowerGenerationPiecewisePoints[g,t][i+1])
            lower = value(m.PowerGenerationPiecewisePoints[g,t][i])
            SD = value(m.ScaledShutdownRampLimit[g,t])
            minP = value(m.MinimumPowerOutput[g,t])

            sd_step = _step_coeff(upper, lower, SD-minP)
            linear_vars = [m.PiecewiseProduction[g,t,i], m.UnitOn[g,t], m.UnitStop[g,t+1],]
            linear_coefs = [-1., (upper-lower), -sd_step,]
            if tightened:
                SU = value(m.ScaledStartupRampLimit[g,t])
                su_step = _step_coeff(upper, lower, SU-minP)
                coef = -max(su_step - sd_step, 0)
                if coef != 0:
                    linear_vars.append(m.UnitStart[g,t])
                    linear_coefs.append(coef)
            return (0, linear_expr(linear_vars=linear_vars, linear_coefs=linear_coefs), None)
        else: ## MinimumUpTime[g] > 1 or we added it in the t == value(m.NumTimePeriods) clause above
            return Constraint.Skip

        
    model.PiecewiseProductionLimits2 = Constraint( model.PiecewiseProductionCostsIndexSet, rule=piecewise_production_limits_rule2 )
    
    _basic_production_costs_constr(model)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_3bin_relaxed_stop_vars', 'garver_2bin_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            })
def KOW_production_costs(model):

    '''
    this is the (more ideal) production cost model introducted by:

    Ben Knueven, Jim Ostrowski, and Jean-Paul Watson. Exploiting identical
    generators in unit commitment. IEEE Transactions on Power Systems,
    33(4), 2018.

    equations (19d)--(19h)
    '''
    _KOW_production_costs(model, False)
    


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars','garver_3bin_relaxed_stop_vars', 'garver_2bin_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            })
def KOW_production_costs_tightened(model):

    '''
    this is the (more ideal) production cost model introducted by:

    Ben Knueven, Jim Ostrowski, and Jean-Paul Watson. Exploiting identical
    generators in unit commitment. IEEE Transactions on Power Systems,
    33(4), 2018.

    equations (19d)--(19h) with some tightening for when SU != SD, as mentioned in text
    '''
    _KOW_production_costs(model, True)


def _CW_production_costs_garver(model):
    '''
    base for CW_production_costs_garver and CW_production_costs_garver_tightened
    '''
    ## TODO: does this just work with regular production variables, i.e., MinP*u \leq p \leq MaxP*u?

    # compute the per-generator, per-time period production costs. We'll do this by hand
    def piecewise_production_costs_index_set_generator(m):
        return ((g,t,i) for g in m.ThermalGenerators for t in m.TimePeriods for i in range(1,len(m.PowerGenerationPiecewisePoints[g,t])))
    model.PiecewiseProductionCostsIndexSet = Set(initialize=piecewise_production_costs_index_set_generator, dimen=3)

    model.PiecewiseProductionFrac = Var( model.PiecewiseProductionCostsIndexSet, within=UnitInterval ) ## UnitInterval == [0,1]

    def piecewise_production_sum_rule(m, g, t):
        return sum( m.PowerGenerationPiecewisePoints[g,t][i]*m.PiecewiseProductionFrac[g,t,i] for i in range(1,len(m.PowerGenerationPiecewisePoints[g,t]))) == m.PowerGeneratedAboveMinimum[g,t]
    model.PiecewiseProductionSum = Constraint( model.ThermalGenerators, model.TimePeriods, rule=piecewise_production_sum_rule )

    def piecewise_production_frac_limits_rule(m, g, t):
        return sum( m.PiecewiseProductionFrac[g,t,i] for i in range(1,len(m.PowerGenerationPiecewisePoints[g,t]))) <= m.UnitOn[g,t]
    model.PiecewiseProductionFracLimits = Constraint( model.ThermalGenerators, model.TimePeriods, rule=piecewise_production_frac_limits_rule )

    model.ProductionCost = Var( model.SingleFuelGenerators, model.TimePeriods, within=Reals )

    def piecewise_production_costs_rule(m, g, t):
        return m.ProductionCost[g,t] == sum( (_production_cost_function(m, g, t, i))*m.PiecewiseProductionFrac[g,t,i] for i in range(1, len(m.PowerGenerationPiecewisePoints[g,t])))

    model.ProductionCostConstr = Constraint( model.SingleFuelGenerators, model.TimePeriods, rule=piecewise_production_costs_rule )

    _compute_total_production_cost(model)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'power_vars': ['garver_power_vars',],
                                            })
def CW_production_costs_garver(model):

    '''
    This is the ideal SOS2-type model for production costs proposed in:

    Chen, Y. and Wang, F (2017). MIP formulation improvedment for large scale
    security constrained unit commitment with configuration based combined
    cycle modeling. Electric Power Systems Research 148 (2017) pp. 147-154

    which is modified to work with garver (power above minimum) production variables.
    '''
    _CW_production_costs_garver(model)


def _SLL_production_costs(model, ideal=True):
    '''
    Base SSL_production_costs_garver
    '''

    ## TODO: does this just work with regular production variables, i.e., MinP*u \leq p \leq MaxP*u?


    # compute the per-generator, per-time period production costs. We'll do this by hand
    def piecewise_production_costs_index_set_generator(m):
        return ((g,t,i) for g in m.ThermalGenerators for t in m.TimePeriods for i in range(len(m.PowerGenerationPiecewisePoints[g,t])))
    model.PiecewiseProductionCostsIndexSet = Set(initialize=piecewise_production_costs_index_set_generator, dimen=3)

    model.PiecewiseProductionFrac = Var( model.PiecewiseProductionCostsIndexSet, within=UnitInterval ) ## UnitInterval == [0,1]

    def piecewise_production_sum_rule(m, g, t):
        return sum( m.PowerGenerationPiecewisePoints[g,t][i]*m.PiecewiseProductionFrac[g,t,i] for i in range(len(m.PowerGenerationPiecewisePoints[g,t]))) == m.PowerGeneratedAboveMinimum[g,t]
    model.PiecewiseProductionSum = Constraint( model.ThermalGenerators, model.TimePeriods, rule=piecewise_production_sum_rule )

    def piecewise_production_frac_limits_rule(m, g, t):
        return sum( m.PiecewiseProductionFrac[g,t,i] for i in range(len(m.PowerGenerationPiecewisePoints[g,t]))) \
                    == (m.UnitOn[g,t] if ideal else 1)
    model.PiecewiseProductionFracLimits = Constraint( model.ThermalGenerators, model.TimePeriods, rule=piecewise_production_frac_limits_rule )

    model.ProductionCost = Var( model.SingleFuelGenerators, model.TimePeriods, within=Reals )

    def piecewise_production_costs_rule(m, g, t):
        return m.ProductionCost[g,t] == sum( (_production_cost_function(m, g, t, i))*m.PiecewiseProductionFrac[g,t,i] for i in range(len(m.PowerGenerationPiecewisePoints[g,t])))

    model.ProductionCostConstr = Constraint( model.SingleFuelGenerators, model.TimePeriods, rule=piecewise_production_costs_rule )

    _compute_total_production_cost(model)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'power_vars': None,
                                            })
def SOS2_production_costs(model):
    '''
    Based on SOS2 model
    '''
    _SLL_production_costs(model, False)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'power_vars': None,
                                            })
def SLL_production_costs(model):
    '''
    Based on SOS2 model S_2 from

    Sridhar, S. Linderoth, J., and Luedtke, J. Locally ideal formulations with
    indicator variables. Operation Research Letters 41 (2013) 627-632
    '''
    _SLL_production_costs(model)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars', 'garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            })
def SLL_production_costs_tightened(model):
    '''
    Based on SOS2 model S_2 from

    Sridhar, S. Linderoth, J., and Luedtke, J. Locally ideal formulations with
    indicator variables. Operation Research Letters 41 (2013) 627-632
    
    with additional tightening based on SU/SD
    '''

    _SLL_production_costs(model)

    @lru_cache()
    def _get_susd_upper_index(piecewise_points, su_sd):

        if su_sd == piecewise_points[0]:
            return 0
        for i in range(len(piecewise_points)-1):
            if piecewise_points[i] < su_sd <= piecewise_points[i+1]:
                return i+1
        return None ## indicates we ran off the end...which shouldn't happen

    def piecewise_production_frac_limits_startup_rule0(m, g, t):
        upper_limit = _get_susd_upper_index(tuple(m.PowerGenerationPiecewisePoints[g,t]), value(m.ScaledStartupRampLimit[g,t]-m.MinimumPowerOutput[g,t]))
        return sum( m.PiecewiseProductionFrac[g,t,i] for i in range(upper_limit+1) ) >= m.UnitStart[g,t]
    model.PiecewiseProductionSumStartup0 = Constraint( model.ThermalGenerators, model.TimePeriods, rule=piecewise_production_frac_limits_startup_rule0 )

    def piecewise_production_frac_limits_startup_rule1(m, g, t):
        upper_limit = _get_susd_upper_index(tuple(m.PowerGenerationPiecewisePoints[g,t]), value(m.ScaledStartupRampLimit[g,t]-m.MinimumPowerOutput[g,t]))
        if m.ScaledStartupRampLimit[g,t]-m.MinimumPowerOutput[g,t] == m.PowerGenerationPiecewisePoints[g,t][upper_limit]:
            return Constraint.Skip
        frac_from_before = 1 - (m.ScaledStartupRampLimit[g,t]-m.MinimumPowerOutput[g,t] - m.PowerGenerationPiecewisePoints[g,t][upper_limit-1])/(m.PowerGenerationPiecewisePoints[g,t][upper_limit] - m.PowerGenerationPiecewisePoints[g,t][upper_limit-1])
        return sum( m.PiecewiseProductionFrac[g,t,i] for i in range(upper_limit) ) >= frac_from_before* m.UnitStart[g,t]
    model.PiecewiseProductionSumStartup1 = Constraint( model.ThermalGenerators, model.TimePeriods, rule=piecewise_production_frac_limits_startup_rule1 )

    def piecewise_production_frac_limits_shutdown_rule0(m, g, t):
        if t < value(m.NumTimePeriods):
            upper_limit = _get_susd_upper_index(tuple(m.PowerGenerationPiecewisePoints[g,t]), value(m.ScaledShutdownRampLimit[g,t]-m.MinimumPowerOutput[g,t]))
            return sum( m.PiecewiseProductionFrac[g,t,i] for i in range(upper_limit+1) ) >= m.UnitStop[g,t+1]
        else:
            return Constraint.Skip
    model.PiecewiseProductionSumShutdown = Constraint( model.ThermalGenerators, model.TimePeriods, rule=piecewise_production_frac_limits_shutdown_rule0 )

    def piecewise_production_frac_limits_shutdown_rule1(m, g, t):
        if t >= value(m.NumTimePeriods):
            return Constraint.Skip
        upper_limit = _get_susd_upper_index(tuple(m.PowerGenerationPiecewisePoints[g,t]), value(m.ScaledShutdownRampLimit[g,t]-m.MinimumPowerOutput[g,t]))
        if m.ScaledShutdownRampLimit[g,t]-m.MinimumPowerOutput[g,t] == m.PowerGenerationPiecewisePoints[g,t][upper_limit]:
            return Constraint.Skip
        frac_from_before = 1 - (m.ScaledShutdownRampLimit[g,t]-m.MinimumPowerOutput[g,t] - m.PowerGenerationPiecewisePoints[g,t][upper_limit-1])/(m.PowerGenerationPiecewisePoints[g,t][upper_limit] - m.PowerGenerationPiecewisePoints[g,t][upper_limit-1])
        return sum( m.PiecewiseProductionFrac[g,t,i] for i in range(upper_limit) ) >= frac_from_before* m.UnitStop[g,t+1]
    model.PiecewiseProductionSumShutdown1 = Constraint( model.ThermalGenerators, model.TimePeriods, rule=piecewise_production_frac_limits_shutdown_rule1 )

# END production models which use an auxilary variable for each piecewise segment

# BEGIN production models which use one auxilary variable per time step to 
#       represent the cost and construct the lower convex envelope

def _basic_production_costs_envelope_garver(model):
    '''
    nearly everything that's needed for the basic 
    convex envelope production costs
    '''

    # compute the per-generator, per-time period production costs. We'll do this by hand
    def piecewise_production_costs_index_set_generator(m):
        return ((g,t,i) for g in m.ThermalGenerators for t in m.TimePeriods for i in range(len(m.PowerGenerationPiecewisePoints[g,t])-1))
    model.PiecewiseProductionCostsIndexSet = Set(initialize=piecewise_production_costs_index_set_generator, dimen=3)

    # variable on which we'll build the convex envelopes
    # NOTE: this assumes the cost to operate the generator at minimum power has been moved to the
    #       coeffecient on the UnitOn variable, such that 0,0 is the starting point for this curve,
    #       and hence >= 0 is a valid inequality
    model.ProductionCost = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'power_vars': None,
                                            })
def basic_production_costs_envelope(model):

    '''
    this is a basic production cost model with an additional variable
    which is the convex envelope of the production costs
    ''' 
    
    _basic_production_costs_envelope_garver(model)

    def piecewise_production_cost_rule(m, g, t, i):

        x0 = m.PowerGenerationPiecewisePoints[g,t][i]
        x1 = m.PowerGenerationPiecewisePoints[g,t][i+1]
        y0 = _production_cost_function(m,g,t,i)
        y1 = _production_cost_function(m,g,t,i+1)
        slope = (y1 - y0)/ (x1 - x0) 
        intercept = -slope*x0 + y0
        return m.ProductionCost[g,t] >= slope*m.PowerGeneratedAboveMinimum[g,t] + intercept

    model.ProductionCostConstr = Constraint(model.PiecewiseProductionCostsIndexSet, rule=piecewise_production_cost_rule)

    _compute_total_production_cost(model)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': None,
                                            'power_vars': None,
                                            })
def HB_production_costs(model):

    '''
    This is the locally ideal formulation using a convex envelope from

    Hua, B. and Baldick, R. (2017) A Convex Primal Formulation for Convex Hull Pricing
    IEEE Transactions on Power Systems, Vol. 32, No. 5, Sept. 2017.
    ''' 
    
    #TODO: can this be strengthened with the addition of the UnitStart, UnitStop variables? 
    _basic_production_costs_envelope_garver(model)


    def piecewise_production_cost_rule(m, g, t, i):

        x0 = m.PowerGenerationPiecewisePoints[g,t][i]
        x1 = m.PowerGenerationPiecewisePoints[g,t][i+1]
        y0 = _production_cost_function(m,g,t,i)
        y1 = _production_cost_function(m,g,t,i+1)
        slope = (y1 - y0)/ (x1 - x0) 
        intercept = -slope*x0 + y0

        # this will be good regardless
        return m.ProductionCost[g,t] >= slope*m.PowerGeneratedAboveMinimum[g,t] + intercept*m.UnitOn[g,t]

    model.ProductionCostConstr = Constraint(model.PiecewiseProductionCostsIndexSet, rule=piecewise_production_cost_rule)

    _compute_total_production_cost(model)
