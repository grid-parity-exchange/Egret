#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## system variables and constraints
from pyomo.environ import *
import math

from .uc_utils import add_model_attr
from .power_vars import _add_reactive_power_vars
from .generation_limits import _add_reactive_limits

import numpy as np
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen
import egret.data.data_utils as data_utils
import egret.common.lazy_ptdf_utils as lpu

from egret.model_library.defn import BasePointType, CoordinateType, ApproximationType
from math import pi

component_name = 'power_balance'

def _copperplate_relax_network_model(block,tm):
    m = block.model()

    ## this is not the "real" gens by bus, but the
    gens_by_bus = block.gens_by_bus

    ### declare (and fix) the loads at the buses
    bus_p_loads = {b: value(m.Demand[b,tm]) for b in m.Buses}

    ## index of net injections from the UC model
    libbus.declare_var_pl(block, m.Buses, initialize=bus_p_loads)
    block.pl.fix()

    bus_gs_fixed_shunts = m._bus_gs_fixed_shunts

    ### declare the p balance
    libbus.declare_eq_p_balance_ed(model=block,
                                   index_set=m.Buses,
                                   bus_p_loads=bus_p_loads,
                                   gens_by_bus=gens_by_bus,
                                   bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                   relax_balance = True,
                                   )


def _copperplate_approx_network_model(block,tm):
    m = block.model()

    ## this is not the "real" gens by bus, but the
    gens_by_bus = block.gens_by_bus

    ### declare (and fix) the loads at the buses
    bus_p_loads = {b: value(m.Demand[b,tm]) for b in m.Buses}

    ## index of net injections from the UC model
    libbus.declare_var_pl(block, m.Buses, initialize=bus_p_loads)
    block.pl.fix()

    bus_gs_fixed_shunts = m._bus_gs_fixed_shunts

    ### declare the p balance
    libbus.declare_eq_p_balance_ed(model=block,
                                   index_set=m.Buses,
                                   bus_p_loads=bus_p_loads,
                                   gens_by_bus=gens_by_bus,
                                   bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                   )

def _ptdf_dcopf_network_model(block,tm):
    m = block.model()

    buses = m._buses
    branches = m._branches
    ptdf_options = m._ptdf_options

    branches_in_service = tuple(l for l in m.TransmissionLines if not value(m.LineOutOfService[l,tm]))
    ## this will serve as a key into our dict of PTDF matricies,
    ## so that we can avoid recalculating them each time step
    ## with the same network topology
    branches_out_service = tuple(l for l in m.TransmissionLines if value(m.LineOutOfService[l,tm]))

    gens_by_bus = block.gens_by_bus

    ### declare (and fix) the loads at the buses
    bus_p_loads = {b: value(m.Demand[b,tm]) for b in m.Buses}

    ## this is not the "real" gens by bus, but the
    ## index of net injections from the UC model
    libbus.declare_var_pl(block, m.Buses, initialize=bus_p_loads)
    block.pl.fix()

    libbus.declare_var_p_nw(block, m.Buses)

    ### get the fixed shunts at the buses
    bus_gs_fixed_shunts = m._bus_gs_fixed_shunts

    ### declare net withdraw expression for use in PTDF power flows
    libbus.declare_eq_p_net_withdraw_at_bus(model=block,
                                            index_set=m.Buses,
                                            bus_p_loads=bus_p_loads,
                                            gens_by_bus=gens_by_bus,
                                            bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                            )

    ### declare the p balance
    libbus.declare_eq_p_balance_ed(model=block,
                                   index_set=m.Buses,
                                   bus_p_loads=bus_p_loads,
                                   gens_by_bus=gens_by_bus,
                                   bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                   )

    ### add "blank" power flow expressions
    libbranch.declare_expr_pf(model=block,
                             index_set=branches_in_service,
                             )

    ### Get the PTDF matrix from cache, from file, or create a new one
    if branches_out_service not in m._PTDFs:
        buses_idx = tuple(buses.keys())

        reference_bus = value(m.ReferenceBus)

        PTDF = data_utils.get_ptdf_potentially_from_file(ptdf_options, branches_in_service, buses_idx)
        
        ## NOTE: For now, just use a flat-start for unit commitment
        if PTDF is None:
            PTDF = data_utils.PTDFMatrix(branches, buses, reference_bus, BasePointType.FLATSTART, branches_keys=branches_in_service, buses_keys=buses_idx)

        m._PTDFs[branches_out_service] = PTDF

    else:
        PTDF = m._PTDFs[branches_out_service]

    ### attach the current PTDF object to this block
    block._PTDF = PTDF

    if ptdf_options['lazy']:
        ### add "blank" real power flow limits
        libbranch.declare_ineq_p_branch_thermal_lbub(model=block,
                                                     index_set=branches_in_service,
                                                     branches=branches,
                                                     p_thermal_limits=None,
                                                     approximation_type=None,
                                                     )
        ### add helpers for tracking monitored branches
        lpu.add_monitored_branch_tracker(block)
        
    else: ### add all the dense constraints
        rel_ptdf_tol = m._ptdf_options['rel_ptdf_tol']
        abs_ptdf_tol = m._ptdf_options['abs_ptdf_tol']

        p_max = {k: branches[k]['rating_long_term'] for k in branches_in_service}

        ### declare the branch power flow approximation constraints
        libbranch.declare_eq_branch_power_ptdf_approx(model=block,
                                                      index_set=branches_in_service,
                                                      PTDF=PTDF,
                                                      abs_ptdf_tol=abs_ptdf_tol,
                                                      rel_ptdf_tol=rel_ptdf_tol
                                                      )
        ### declare the real power flow limits
        libbranch.declare_ineq_p_branch_thermal_lbub(model=block,
                                                     index_set=branches_in_service,
                                                     branches=branches,
                                                     p_thermal_limits=p_max,
                                                     approximation_type=ApproximationType.PTDF
                                                     )


def _btheta_dcopf_network_model(block,tm):
    m = block.model()

    buses = m._buses
    branches = m._branches

    branches_in_service = tuple(l for l in m.TransmissionLines if not value(m.LineOutOfService[l,tm]))
    ## this will serve as a key into our dict of PTDF matricies,
    ## so that we can avoid recalculating them each time step
    ## with the same network topology
    branches_out_service = tuple(l for l in m.TransmissionLines if value(m.LineOutOfService[l,tm]))

    ## need the inlet/outlet relationship given some lines may be out
    inlet_branches_by_bus = dict()
    outlet_branches_by_bus = dict()
    for b in m.Buses:
        inlet_branches_by_bus[b] = list()
        for l in m.LinesTo[b]:
            if l not in branches_out_service:
                inlet_branches_by_bus[b].append(l)
        outlet_branches_by_bus[b] = list()
        for l in m.LinesFrom[b]:
            if l not in branches_out_service:
                outlet_branches_by_bus[b].append(l)

    ## this is not the "real" gens by bus, but the
    ## index of net injections from the UC model
    gens_by_bus = block.gens_by_bus

    ### declare (and fix) the loads at the buses
    bus_p_loads = {b: value(m.Demand[b,tm]) for b in m.Buses}

    libbus.declare_var_pl(block, buses.keys(), initialize=bus_p_loads)
    block.pl.fix()

    ### get the fixed shunts at the buses
    bus_gs_fixed_shunts = m._bus_gs_fixed_shunts

    ### declare the polar voltages
    va_bounds = {k: (-pi, pi) for k in buses.keys()}
    libbus.declare_var_va(block, buses.keys(), initialize=None,
                          bounds=va_bounds
                          )

    ### fix the reference bus
    ref_bus = value(m.ReferenceBus)
    ref_angle = value(m.ReferenceBusAngle)
    block.va[ref_bus].fix(math.radians(ref_angle))

    p_max = {k: branches[k]['rating_long_term'] for k in branches_in_service}
    p_lbub = {k: (-p_max[k],p_max[k]) for k in branches_in_service}
    pf_bounds = p_lbub

    libbranch.declare_var_pf(model=block,
                             index_set=branches_in_service,
                             initialize=None,
                             bounds=pf_bounds
                             )

    ### declare the branch power flow approximation constraints
    libbranch.declare_eq_branch_power_btheta_approx(model=block,
                                                    index_set=branches_in_service,
                                                    branches=branches
                                                    )

    ### declare the p balance
    libbus.declare_eq_p_balance_dc_approx(model=block,
                                          index_set=buses.keys(),
                                          bus_p_loads=bus_p_loads,
                                          gens_by_bus=gens_by_bus,
                                          bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                          inlet_branches_by_bus=inlet_branches_by_bus,
                                          outlet_branches_by_bus=outlet_branches_by_bus,
                                          approximation_type=ApproximationType.BTHETA
                                          )

    ### declare the real power flow limits
    libbranch.declare_ineq_p_branch_thermal_lbub(model=block,
                                                 index_set=branches_in_service,
                                                 branches=branches,
                                                 p_thermal_limits=p_max,
                                                 approximation_type=ApproximationType.BTHETA
                                                 )

    ### declare angle difference limits on interconnected buses
    libbranch.declare_ineq_angle_diff_branch_lbub(model=block,
                                                  index_set=branches_in_service,
                                                  branches=branches,
                                                  coordinate_type=CoordinateType.POLAR
                                                  )

    return block

def _add_system_load_mismatch(model):

    #####################################################
    # load "shedding" can be both positive and negative #
    #####################################################
    model.posLoadGenerateMismatch = Var(model.TimePeriods, within=NonNegativeReals) # load shedding
    model.negLoadGenerateMismatch = Var(model.TimePeriods, within=NonNegativeReals) # over generation

    ## for interfacing with the rest of the model code
    def define_pos_neg_load_generate_mismatch_rule(m, b, t):
        if b == value(m.ReferenceBus):
            return m.posLoadGenerateMismatch[t] - m.negLoadGenerateMismatch[t]
        else:
            return 0
    model.LoadGenerateMismatch = Expression(model.Buses, model.TimePeriods, rule = define_pos_neg_load_generate_mismatch_rule )

    # the following constraints are necessarily, at least in the case of CPLEX 12.4, to prevent
    # the appearance of load generation mismatch component values in the range of *negative* e-5.
    # what these small negative values do is to cause the optimal objective to be a very large negative,
    # due to obviously large penalty values for under or over-generation. JPW would call this a heuristic
    # at this point, but it does seem to work broadly. we tried a single global constraint, across all
    # buses, but that failed to correct the problem, and caused the solve times to explode.
    
    def pos_load_generate_mismatch_tolerance_rule(m):
       return sum((m.posLoadGenerateMismatch[t] for t in m.TimePeriods)) >= 0.0
    model.PosLoadGenerateMismatchTolerance = Constraint(rule=pos_load_generate_mismatch_tolerance_rule)
    
    def neg_load_generate_mismatch_tolerance_rule(m):
       return sum((m.negLoadGenerateMismatch[t] for t in m.TimePeriods)) >= 0.0
    model.NegLoadGenerateMismatchTolerance = Constraint(rule=neg_load_generate_mismatch_tolerance_rule)

    def compute_load_mismatch_cost_rule(m, t):
        return m.LoadMismatchPenalty*m.TimePeriodLengthHours*(m.posLoadGenerateMismatch[t] + m.negLoadGenerateMismatch[t]) 
    model.LoadMismatchCost = Expression(model.TimePeriods, rule=compute_load_mismatch_cost_rule)

def _add_load_mismatch(model):

    def load_shedding_bus_times(m):
        for b in m.Buses:
            for t in m.TimePeriods:
                if value(m.Demand[b,t]) > 0:
                    yield b,t
    model.LoadSheddingBusTimes = Set(dimen=2, initialize=load_shedding_bus_times)

    load_shedding_times_per_bus = {b: list() for b in model.Buses}
    for b,t in model.LoadSheddingBusTimes:
        load_shedding_times_per_bus[b].append(t)

    def load_shedding_bounds_times(m,b,t):
        ## NOTE: depends on above which ensures no engative demand
        return (0, value(m.Demand[b,t]))
    model.LoadShedding = Var(model.LoadSheddingBusTimes, within=NonNegativeReals, bounds=load_shedding_bounds_times) # load shedding

    over_gen_maxes = {}
    over_gen_times_per_bus = {b: list() for b in model.Buses}
    for b in model.Buses:
        gen = sum(value(model.MaximumPowerOutput[g]) for g in model.ThermalGeneratorsAtBus[b])
        for t in model.TimePeriods:
            total_gen = gen + sum(value(model.MinNondispatchablePower[n,t]) for n in model.NondispatchableGeneratorsAtBus[b])
            total_gen -= value(model.Demand[b,t])
            if total_gen > 0:
                over_gen_maxes[b,t] = total_gen
                over_gen_times_per_bus[b].append(t)

    model.OverGenerationBusTimes = Set(dimen=2, initialize=over_gen_maxes.keys())

    def get_over_gen_bounds(m,b,t):
        return (0,over_gen_maxes[b,t])
    model.OverGeneration = Var(model.OverGenerationBusTimes, within=NonNegativeReals, bounds=get_over_gen_bounds) # over generation
    
    # the following constraints are necessarily, at least in the case of CPLEX 12.4, to prevent
    # the appearance of load generation mismatch component values in the range of *negative* e-5.
    # what these small negative values do is to cause the optimal objective to be a very large negative,
    # due to obviously large penalty values for under or over-generation. JPW would call this a heuristic
    # at this point, but it does seem to work broadly. we tried a single global constraint, across all
    # buses, but that failed to correct the problem, and caused the solve times to explode.


    def pos_load_generate_mismatch_tolerance_rule(m, b):
        if load_shedding_times_per_bus[b]:
            return sum(m.LoadShedding[b,t] for t in load_shedding_times_per_bus[b]) >= 0
        else:
            return Constraint.Feasible
    model.PosLoadGenerateMismatchTolerance = Constraint(model.Buses, 
                                                                rule=pos_load_generate_mismatch_tolerance_rule)
    
    def neg_load_generate_mismatch_tolerance_rule(m, b):
        if over_gen_times_per_bus[b]:
            return sum(m.OverGeneration[b,t] for t in m.TimePeriods if (b,t) in m.OverGenerationBusTimes) >= 0
        else:
            return Constraint.Feasible
    model.NegLoadGenerateMismatchTolerance = Constraint(model.Buses,
                                                                rule=neg_load_generate_mismatch_tolerance_rule)
    
    #####################################################
    # load "shedding" can be both positive and negative #
    #####################################################
    model.LoadGenerateMismatch = Expression(model.Buses, model.TimePeriods)
    for b in model.Buses:
        for t in model.TimePeriods:
            model.LoadGenerateMismatch[b,t].expr = 0
    for b,t in model.LoadSheddingBusTimes:
        model.LoadGenerateMismatch[b,t].expr += model.LoadShedding[b,t]
    for b,t in model.OverGenerationBusTimes:
        model.LoadGenerateMismatch[b,t].expr -= model.OverGeneration[b,t]

    model.LoadMismatchCost = Expression(model.TimePeriods)
    for t in model.TimePeriods:
        model.LoadMismatchCost[t].expr = 0
    for b,t in model.LoadSheddingBusTimes:
        model.LoadMismatchCost[t].expr += model.LoadMismatchPenalty*model.TimePeriodLengthHours*model.LoadShedding[b,t]
    for b,t in model.OverGenerationBusTimes:
        model.LoadMismatchCost[t].expr += model.LoadMismatchPenalty*model.TimePeriodLengthHours*model.OverGeneration[b,t]


def _add_q_load_mismatch(model):

    #####################################################
    # load "shedding" can be both positive and negative #
    #####################################################
    model.LoadGenerateMismatchReactive = Var(model.Buses, model.TimePeriods, within=Reals)
    model.posLoadGenerateMismatchReactive = Var(model.Buses, model.TimePeriods, within=NonNegativeReals) # load shedding
    model.negLoadGenerateMismatchReactive = Var(model.Buses, model.TimePeriods, within=NonNegativeReals) # over generation
    
    def define_pos_neg_load_generate_mismatch_rule_reactive(m, b, t):
        return m.posLoadGenerateMismatchReactive[b, t] - m.negLoadGenerateMismatchReactive[b, t] \
                == m.LoadGenerateMismatchReactive[b, t]
    model.DefinePosNegLoadGenerateMismatchReactive = Constraint(model.Buses, model.TimePeriods, rule = define_pos_neg_load_generate_mismatch_rule_reactive)

    # the following constraints are necessarily, at least in the case of CPLEX 12.4, to prevent
    # the appearance of load generation mismatch component values in the range of *negative* e-5.
    # what these small negative values do is to cause the optimal objective to be a very large negative,
    # due to obviously large penalty values for under or over-generation. JPW would call this a heuristic
    # at this point, but it does seem to work broadly. we tried a single global constraint, across all
    # buses, but that failed to correct the problem, and caused the solve times to explode.
    
    def pos_load_generate_mismatch_tolerance_rule_reactive(m, b):
       return sum((m.posLoadGenerateMismatchReactive[b,t] for t in m.TimePeriods)) >= 0.0
    model.PosLoadGenerateMismatchToleranceReactive = Constraint(model.Buses, 
                                                                rule=pos_load_generate_mismatch_tolerance_rule_reactive)
    
    def neg_load_generate_mismatch_tolerance_rule_reactive(m, b):
       return sum((m.negLoadGenerateMismatchReactive[b,t] for t in m.TimePeriods)) >= 0.0
    model.NegLoadGenerateMismatchToleranceReactive = Constraint(model.Buses,
                                                                rule=neg_load_generate_mismatch_tolerance_rule_reactive)

    def compute_q_load_mismatch_cost_rule(m, t):
        return m.LoadMismatchPenaltyReactive*m.TimePeriodLengthHours*sum(
                    m.posLoadGenerateMismatchReactive[b, t] + m.negLoadGenerateMismatchReactive[b, t] for b in m.Buses) 
    model.LoadMismatchCostReactive = Expression(model.TimePeriods, rule=compute_q_load_mismatch_cost_rule)

def _add_blank_load_mismatch(model):
    model.LoadGenerateMismatch = Param(model.Buses, model.TimePeriods, default=0.)
    model.posLoadGenerateMismatch = Param(model.Buses, model.TimePeriods, default=0.)
    model.negLoadGenerateMismatch = Param(model.Buses, model.TimePeriods, default=0.)
    model.LoadMismatchCost = Param(model.TimePeriods, default=0.)

def _add_blank_q_load_mismatch(model):
    model.LoadGenerateMismatchReactive = Param(model.Buses, model.TimePeriods, default=0.)
    model.posLoadGenerateMismatchReactive = Param(model.Buses, model.TimePeriods, default=0.)
    model.negLoadGenerateMismatchReactive = Param(model.Buses, model.TimePeriods, default=0.)
    model.LoadMismatchCostReactive = Param(model.TimePeriods, default=0.)

## helper defining real power injection at a bus
def _get_pg_expr_rule(t):
    def pg_expr_rule(block,b):
        m = block.model()
        # bus b, time t (S)
        return sum(m.PowerGenerated[g, t] for g in m.ThermalGeneratorsAtBus[b]) \
                + sum(m.PowerOutputStorage[s, t] for s in m.StorageAtBus[b])\
                - sum(m.PowerInputStorage[s, t] for s in m.StorageAtBus[b])\
                + sum(m.NondispatchablePowerUsed[g, t] for g in m.NondispatchableGeneratorsAtBus[b]) \
                + m.LoadGenerateMismatch[b,t]
    return pg_expr_rule

## helper defining reacative power injection at a bus
def _get_qg_expr_rule(t):
    def qg_expr_rule(block,b):
        m = block.model()
        # bus b, time t (S)
        return sum(m.ReactivePowerGenerated[g, t] for g in m.ThermalGeneratorsAtBus[b]) \
            + m.LoadGenerateMismatchReactive[b,t]
    return qg_expr_rule

## Defines generic interface for egret tramsmission models
def _add_egret_power_flow(model, network_model_builder, reactive_power=False, slacks=True):

    ## save flag for objective
    model.reactive_power = reactive_power

    system_load_mismatch = (network_model_builder == _copperplate_approx_network_model)

    if slacks:
        if system_load_mismatch:
            _add_system_load_mismatch(model)
        else:
            _add_load_mismatch(model)
    else:
        if system_load_mismatch:
            _add_blank_system_load_mismatch(model)
        else:
            _add_blank_load_mismatch(model)

    if reactive_power:
        if system_load_mismatch:
            raise Exception("Need to implement system mismatch for reactive power")
        _add_reactive_power_vars(model)
        _add_reactive_limits(model)
        if slacks:
            _add_q_load_mismatch(model)
        else:
            _add_blank_q_load_mistmatch(model)

    # for transmission network
    model.TransmissionBlock = Block(model.TimePeriods, concrete=True)

    for tm in model.TimePeriods:
        b = model.TransmissionBlock[tm]
        ## this creates a fake bus generator for all the
        ## appropriate injection/withdraws from the unit commitment
        ## model
        b.pg = Expression(model.Buses, rule=_get_pg_expr_rule(tm))
        if reactive_power:
            b.qg = Expression(model.Buses, rule=_get_qg_expr_rule(tm))
        b.gens_by_bus = {bus : [bus] for bus in model.Buses}
        network_model_builder(b,tm)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'power_vars': None,
                                            'non_dispatchable_vars': None,
                                            'storage_service': None,
                                            })
def copperplate_power_flow(model, slacks=True):
    _add_egret_power_flow(model, _copperplate_approx_network_model, reactive_power=False, slacks=slacks)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'power_vars': None,
                                            'non_dispatchable_vars': None,
                                            'storage_service': None,
                                            })
def copperplate_relaxed_power_flow(model, slacks=True):
    _add_egret_power_flow(model, _copperplate_relax_network_model, reactive_power=False, slacks=slacks)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'power_vars': None,
                                            'non_dispatchable_vars': None,
                                            'storage_service': None,
                                            })
def ptdf_power_flow(model, slacks=True):
    model._PTDFs = dict()
    _add_egret_power_flow(model, _ptdf_dcopf_network_model, reactive_power=False, slacks=slacks)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'power_vars': None,
                                            'non_dispatchable_vars': None,
                                            'storage_service': None,
                                            })
def btheta_power_flow(model, slacks=True):
    _add_egret_power_flow(model, _btheta_dcopf_network_model, reactive_power=False, slacks=slacks)


@add_model_attr(component_name, requires = {'data_loader': None,
                                            'power_vars': None,
                                            'non_dispatchable_vars': None,
                                            'storage_service': None,
                                            })
def power_balance_constraints(model, slacks=True):
    '''
    adds the demand and network constraints to the model
    '''
    model.reactive_power = False

    # system variables
    # amount of power flowing along each line, at each time period
    def line_power_bounds_rule(m, l, t):
        return (-m.ThermalLimit[l], m.ThermalLimit[l])
    model.LinePower = Var(model.TransmissionLines, model.TimePeriods, bounds=line_power_bounds_rule)
    
    # voltage angles at the buses (S) (lock the first bus at 0) in radians
    model.Angle = Var(model.Buses, model.TimePeriods, within=Reals, bounds=(-3.14159265,3.14159265))
    
    def fix_first_angle_rule(m,t):
        first_bus = list(sorted(m.Buses))[0]
        return m.Angle[first_bus,t] == 0.0
    model.FixFirstAngle = Constraint(model.TimePeriods, rule=fix_first_angle_rule)

    def line_power_rule(m, l, t):
        if value(m.LineOutOfService[l,t]):
            return m.LinePower[l,t] == 0.0
        else:
            return m.LinePower[l,t] == (m.Angle[m.BusFrom[l], t] - m.Angle[m.BusTo[l], t]) / m.Impedence[l]
    model.CalculateLinePower = Constraint(model.TransmissionLines, model.TimePeriods, rule=line_power_rule)
    
    def interface_from_limit_rule(m,i,t):
        return sum(m.LinePower[l,t] for l in m.InterfaceLines[i]) <= m.InterfaceFromLimit[i]
    model.InterfaceFromLimitConstr = Constraint(model.Interfaces, model.TimePeriods, rule=interface_from_limit_rule)

    def interface_to_limit_rule(m,i,t):
        return sum(m.LinePower[l,t] for l in m.InterfaceLines[i]) >= -m.InterfaceToLimit[i]
    model.InterfaceToLimitConstr = Constraint(model.Interfaces, model.TimePeriods, rule=interface_to_limit_rule)
    
    if slacks:
        _add_load_mismatch(model)
    else:
        _add_blank_load_mismatch(model)
    
    # Power balance at each node (S)
    def power_balance(m, b, t):
        # bus b, time t (S)
        return sum(m.PowerGenerated[g, t] for g in m.ThermalGeneratorsAtBus[b]) \
                + sum(m.PowerOutputStorage[s, t] for s in m.StorageAtBus[b])\
                - sum(m.PowerInputStorage[s, t] for s in m.StorageAtBus[b])\
                + sum(m.NondispatchablePowerUsed[g, t] for g in m.NondispatchableGeneratorsAtBus[b]) \
                + sum(m.LinePower[l,t] for l in m.LinesTo[b]) \
                - sum(m.LinePower[l,t] for l in m.LinesFrom[b]) \
                + m.LoadGenerateMismatch[b,t] \
                == m.Demand[b, t] 
    model.PowerBalance = Constraint(model.Buses, model.TimePeriods, rule=power_balance)

    return
## end power_balance_constraints
