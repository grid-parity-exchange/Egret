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

from .uc_utils import add_model_attr, get_linear_expr
from .reserve_vars import check_reserve_requirement
component_name = 'objective'

def _compute_1bin_shutdown_costs(model):
    #############################################################
    # compute the per-generator, per-time period shutdown costs #
    #############################################################
    
    def compute_shutdown_costs_rule(m, g, t):
        if t == m.InitialTime:
          return m.ShutdownCost[g, t] >= m.ShutdownFixedCost[g] * (m.UnitOnT0[g] - m.UnitOn[g, t])
        else:
          return m.ShutdownCost[g, t] >= m.ShutdownFixedCost[g] * (m.UnitOn[g, t-1] - m.UnitOn[g, t])
    
    model.ComputeShutdownCosts = Constraint(model.ThermalGenerators, model.TimePeriods, rule=compute_shutdown_costs_rule)

def _1bin_shutdown_costs(model, add_shutdown_cost_var=True):

    if add_shutdown_cost_var:
        model.ShutdownCost = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)

    _compute_1bin_shutdown_costs(model)

def _3bin_shutdown_costs(model, add_shutdown_cost_var=True):

    #############################################################
    # compute the per-generator, per-time period shutdown costs #
    #############################################################
    ## BK -- replaced with UnitStop

    if add_shutdown_cost_var:
        model.ShutdownCost = Var(model.ThermalGenerators, model.TimePeriods, within=Reals)
    
    linear_expr = get_linear_expr(model.UnitStop)

    def compute_shutdown_costs_rule(m, g, t):
        return (linear_expr(
                    linear_vars=[m.ShutdownCost[g,t], m.UnitStop[g,t]],
                    linear_coefs=[-1., m.ShutdownFixedCost[g]]),
                0.)

    model.ComputeShutdownCosts = Constraint(model.ThermalGenerators, model.TimePeriods, rule=compute_shutdown_costs_rule)

def _add_shutdown_costs(model, add_shutdown_cost_var=True):
    #NOTE: we handle shutdown costs in this manner because it's not a 
    #      common point of contention in the literature, and they're 
    #      often zero as is.
    if model.status_vars in ['garver_3bin_vars','garver_3bin_relaxed_stop_vars','garver_2bin_vars', 'ALS_state_transition_vars']:
        _3bin_shutdown_costs(model, add_shutdown_cost_var)
    elif model.status_vars in ['CA_1bin_vars',]:
        _1bin_shutdown_costs(model, add_shutdown_cost_var)
    else:
        raise Exception("Problem adding shutdown costs, cannot identify status_vars for this model")
    

#TODO: this doesn't check if regulation_service is added first, 
#      but this will only happen when there are regulation_services!
@add_model_attr(component_name, requires = {'data_loader': None,
                                            'status_vars': ['garver_3bin_vars', 'CA_1bin_vars', 'garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                            'power_vars': None,
                                            'startup_costs': None,
                                            'production_costs': None,
                                            'power_balance': None,
                                            'reserve_requirement': None,
                                            'storage_service': None,
                                            'ancillary_service': None,
                                            })
def basic_objective(model):
    '''
    adds the objective and shutdown cost formulation to the model
    '''
    
    #############################################
    # constraints for computing cost components #
    #############################################
    
    def compute_no_load_cost_rule(m,g,t):
        return m.MinimumProductionCost[g,t]*m.UnitOn[g,t]*m.TimePeriodLengthHours
    
    model.NoLoadCost = Expression(model.SingleFuelGenerators, model.TimePeriods, rule=compute_no_load_cost_rule)
    
    _add_shutdown_costs(model)

    # compute the total production costs, across all generators and time periods.
    def compute_total_production_cost_rule(m, t):
        return sum(m.ProductionCost[g, t] for g in m.SingleFuelGenerators) + \
                sum(m.DualFuelProductionCost[g,t] for g in m.DualFuelGenerators)
    
    model.TotalProductionCost = Expression(model.TimePeriods, rule=compute_total_production_cost_rule)

    # 
    # Cost computations
    #

    def commitment_stage_cost_expression_rule(m, st):
        cc = sum(sum(m.NoLoadCost[g,t] + m.StartupCost[g,t] for g in m.SingleFuelGenerators) + \
                 sum(m.DualFuelCommitmentCost[g,t] for g in m.DualFuelGenerators) + \
                 sum(m.ShutdownCost[g,t] for g in m.ThermalGenerators) 
             for t in m.CommitmentTimeInStage[st])
        if m.regulation_service:
            cc += sum(m.RegulationCostCommitment[g,t] for g in m.AGC_Generators for t in m.CommitmentTimeInStage[st])
        return cc
    
    model.CommitmentStageCost = Expression(model.StageSet, rule=commitment_stage_cost_expression_rule)

    def compute_reserve_shortfall_cost_rule(m, t):
        return m.ReserveShortfallPenalty*m.TimePeriodLengthHours*m.ReserveShortfall[t]
    model.ReserveShortfallCost = Expression(model.TimePeriods, rule=compute_reserve_shortfall_cost_rule)
    
    def generation_stage_cost_expression_rule(m, st):
        cc = sum(sum(m.ProductionCost[g, t] for g in m.SingleFuelGenerators) + \
                 sum(m.DualFuelProductionCost[g,t] for g in m.DualFuelGenerators)
                for t in m.GenerationTimeInStage[st]) + \
              sum(m.LoadMismatchCost[t] for t in m.GenerationTimeInStage[st]) + \
              sum(m.ReserveShortfallCost[t] for t in m.GenerationTimeInStage[st]) + \
              sum(m.BranchViolationCost[t] for t in m.GenerationTimeInStage[st]) + \
              sum(m.InterfaceViolationCost[t] for t in m.GenerationTimeInStage[st]) + \
              sum(m.ContingencyViolationCost[t] for t in m.GenerationTimeInStage[st]) + \
              sum(m.StorageCost[s,t] for s in m.Storage for t in m.GenerationTimeInStage[st])
        if m.reactive_power:
            cc += sum(m.LoadMismatchCostReactive[t] for t in m.GenerationTimeInStage[st])
        if m.security_constraints:
            cc += sum(m.SecurityConstraintViolationCost[t] for t in m.GenerationTimeInStage[st])
        if m.regulation_service:
            cc += sum(m.RegulationCostGeneration[g,t] for g in m.AGC_Generators for t in m.GenerationTimeInStage[st]) \
                + sum(m.RegulationCostPenalty[t] for t in m.GenerationTimeInStage[st])
        if m.spinning_reserve:
            cc += sum(m.SpinningReserveCostGeneration[g,t] for g in m.ThermalGenerators for t in m.GenerationTimeInStage[st]) \
                + sum(m.SpinningReserveCostPenalty[t] for t in m.GenerationTimeInStage[st])
        if m.non_spinning_reserve:
            cc += sum(m.NonSpinningReserveCostGeneration[g,t] for g in m.NonSpinGenerators for t in m.GenerationTimeInStage[st]) \
                + sum(m.NonSpinningReserveCostPenalty[t] for t in m.GenerationTimeInStage[st])
        if m.supplemental_reserve:
            cc += sum(m.SupplementalReserveCostGeneration[g,t] for g in m.ThermalGenerators for t in m.GenerationTimeInStage[st]) \
                + sum(m.SupplementalReserveCostPenalty[t] for t in m.GenerationTimeInStage[st])
        if m.flexible_ramping:
            cc += sum(m.FlexibleRampingCostPenalty[t] for t in m.GenerationTimeInStage[st])
        return cc
    model.GenerationStageCost = Expression(model.StageSet, rule=generation_stage_cost_expression_rule)
    
    def stage_cost_expression_rule(m, st):
        return m.GenerationStageCost[st] + m.CommitmentStageCost[st]
    model.StageCost = Expression(model.StageSet, rule=stage_cost_expression_rule)
    
    #
    # Objectives
    #
    
    def total_cost_objective_rule(m):
       return sum(m.StageCost[st] for st in m.StageSet)

    model.TotalCostObjective = Objective(rule=total_cost_objective_rule, sense=minimize)

    return
