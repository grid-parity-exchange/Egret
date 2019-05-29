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

from .uc_utils import add_model_attr 
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
    
    def compute_shutdown_costs_rule(m, g, t):
        return m.ShutdownCost[g,t] ==  m.ShutdownFixedCost[g] * (m.UnitStop[g, t])
    
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
    

#TODO: this doesn't check if regulation_services is added first, 
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
    
    def compute_total_no_load_cost_rule(m,t):
        return sum(m.MinimumProductionCost[g]*m.UnitOn[g,t]*m.TimePeriodLengthHours for g in m.ThermalGenerators)
    
    model.TotalNoLoadCost = Expression(model.TimePeriods, rule=compute_total_no_load_cost_rule)
    
    _add_shutdown_costs(model)

    # 
    # Cost computations
    #

    # gather ancillary services added
    regulation = False
    spin = False
    nspin = False
    supp = False
    flex = False
    if hasattr(model, 'regulation_service'):
        regulation = True
    if hasattr(model, 'spinning_reserve'):
        spin = True
    if hasattr(model, 'non_spinning_reserve'):
        nspin = True
    if hasattr(model, 'supplemental_reserve'):
        supp = True
    if hasattr(model, 'flexible_ramping'):
        flex = True
    
    def commitment_stage_cost_expression_rule(m, st):
        cc = sum(m.StartupCost[g,t] + m.ShutdownCost[g,t] for g in m.ThermalGenerators for t in m.CommitmentTimeInStage[st]) + \
             sum(m.TotalNoLoadCost[t] for t in m.CommitmentTimeInStage[st])
        if regulation:
            cc += sum(m.RegulationCostCommitment[g,t] for g in m.AGC_Generators for t in m.CommitmentTimeInStage[st])
        return cc
    
    model.CommitmentStageCost = Expression(model.StageSet, rule=commitment_stage_cost_expression_rule)

    def compute_load_mismatch_cost_rule(m, t):
        return m.LoadMismatchPenalty*m.TimePeriodLengthHours*sum(m.posLoadGenerateMismatch[b, t] + m.negLoadGenerateMismatch[b, t] for b in m.Buses) 
    model.LoadMismatchCost = Expression(model.TimePeriods, rule=compute_load_mismatch_cost_rule)

    if model.reactive_power:
        def compute_q_load_mismatch_cost_rule(m, t):
            return m.LoadMismatchPenaltyReactive*m.TimePeriodLengthHours*sum(
                        m.posLoadGenerateMismatchReactive[b, t] + m.negLoadGenerateMismatchReactive[b, t] for b in m.Buses) 
        model.LoadMismatchCostReactive = Expression(model.TimePeriods, rule=compute_q_load_mismatch_cost_rule)

    def compute_reserve_shortfall_cost_rule(m, t):
        return m.ReserveShortfallPenalty*m.TimePeriodLengthHours*m.ReserveShortfall[t]
    model.ReserveShortfallCost = Expression(model.TimePeriods, rule=compute_reserve_shortfall_cost_rule)
    
    def generation_stage_cost_expression_rule(m, st):
        cc = sum(m.ProductionCost[g, t] for g in m.ThermalGenerators for t in m.GenerationTimeInStage[st]) + \
              sum(m.LoadMismatchCost[t] for t in m.GenerationTimeInStage[st]) + \
              sum(m.ReserveShortfallCost[t] for t in m.GenerationTimeInStage[st])
        cc += sum(m.StorageCost[s,t] for s in m.Storage for t in m.GenerationTimeInStage[st])
        if m.reactive_power:
            cc += sum(m.LoadMismatchCostReactive[t] for t in m.GenerationTimeInStage[st])
        if regulation:
            cc += sum(m.RegulationCostGeneration[g,t] for g in m.AGC_Generators for t in m.GenerationTimeInStage[st]) \
                + sum(m.RegulationCostPenalty[t] for t in m.GenerationTimeInStage[st])
        if spin:
            cc += sum(m.SpinningReserveCostGeneration[g,t] for g in m.ThermalGenerators for t in m.GenerationTimeInStage[st]) \
                + sum(m.SpinningReserveCostPenalty[t] for t in m.GenerationTimeInStage[st])
        if nspin:
            cc += sum(m.NonSpinningReserveCostGeneration[g,t] for g in m.NonSpinGenerators for t in m.GenerationTimeInStage[st]) \
                + sum(m.NonSpinningReserveCostPenalty[t] for t in m.GenerationTimeInStage[st])
        if supp:
            cc += sum(m.SupplementalReserveCostGeneration[g,t] for g in m.ThermalGenerators for t in m.GenerationTimeInStage[st]) \
                + sum(m.SupplementalReserveCostPenalty[t] for t in m.GenerationTimeInStage[st])
        if flex:
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
