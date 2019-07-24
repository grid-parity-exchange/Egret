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

from egret.model_library.defn import BasePointType, CoordinateType, ApproximationType
from math import pi

component_name = 'power_balance'

def _copperplate_approx_network_model(md,block,tm,td):
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

def _lazy_ptdf_dcopf_network_model(md,block,tm,td):
    m = block.model()

    buses = m._buses
    branches = m._branches

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

    ### get the fixed shunts at the buses
    bus_gs_fixed_shunts = m._bus_gs_fixed_shunts

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

    ### add "blank" real power flow limits
    libbranch.declare_ineq_p_branch_thermal_lbub(model=block,
                                                 index_set=branches_in_service,
                                                 branches=branches,
                                                 p_thermal_limits=None,
                                                 approximation_type=None,
                                                 )
    ### end initial set-up
    if branches_out_service not in m._PTDFs_dict:
        ## make a new PTDF matrix for this topology
        _PTDF_dict = dict()
        ## to keep things in order
        buses_idx = tuple(buses.keys())
        branches_idx = branches_in_service

        reference_bus = md.data['system']['reference_bus']

        ## calculate PTDFs, but don't do anything with them (for now..)
        #from pyutilib.misc.timing import TicTocTimer
        #timer = TicTocTimer()
        #timer.tic('starting PTDF calculation')
        PTDFM = tx_calc.calculate_ptdf(branches,buses,branches_idx,buses_idx,reference_bus,BasePointType.FLATSTART)
        #timer.toc('done')

        ## store some information we'll need when iterating on the model object
        _PTDF_dict['PTDFM'] = PTDFM
        _PTDF_dict['buses_idx'] = buses_idx
        _PTDF_dict['branches_idx'] = branches_idx
        _PTDF_dict['branch_limits'] = np.array([ branches[branch]['rating_long_term'] for branch in branches_idx ])

        m._PTDFs_dict[branches_out_service] = _PTDF_dict

    ## create a pointer on this block to all this PTDF data,
    ## for easy iteration within the solve loop
    block._PTDF_dict = m._PTDFs_dict[branches_out_service]

    ## this expression is specific to each block
    block._PTDF_bus_nw_exprs = \
      [ block.pl[bus] + bus_gs_fixed_shunts[bus] - \
      sum(block.pg[g] for g in gens_by_bus[bus]) for bus in buses]
    block._PTDF_bus_p_loads = bus_p_loads



def _ptdf_dcopf_network_model(md,block,tm,td):

    m = block.model()
    rel_ptdf_tol = m._ptdf_options_dict['rel_ptdf_tol']
    abs_ptdf_tol = m._ptdf_options_dict['abs_ptdf_tol']

    buses = m._buses
    branches = m._branches
    shunts = m._shunts

    branches_in_service = tuple(l for l in m.TransmissionLines if not value(m.LineOutOfService[l,tm]))

    ## this will serve as a key into our dict of PTDF matricies,
    ## so that we can avoid recalculating them each time step
    ## with the same network topology
    branches_out_service = tuple(l for l in m.TransmissionLines if value(m.LineOutOfService[l,tm]))

    ## this is not the "real" gens by bus, but the
    ## index of net injections from the UC model
    gens_by_bus = block.gens_by_bus

    ### declare (and fix) the loads at the buses
    bus_p_loads = {b: value(m.Demand[b,tm]) for b in m.Buses}

    libbus.declare_var_pl(block, buses.keys(), initialize=bus_p_loads)
    block.pl.fix()

    ### get the fixed shunts at the buses
    bus_gs_fixed_shunts = m._bus_gs_fixed_shunts


    if branches_out_service not in m._PTDFs_dict:
        ## make a new PTDF matrix for this topology

        _PTDF_dict = dict()
        ## to keep things in order
        buses_idx = tuple(buses.keys())
        branches_idx = branches_in_service

        reference_bus = md.data['system']['reference_bus']

        ## calculate PTDFs, but don't do anything with them (for now..)
        #from pyutilib.misc.timing import TicTocTimer
        #timer = TicTocTimer()
        #timer.tic('starting PTDF calculation')
        PTDFM = tx_calc.calculate_ptdf(branches,buses,branches_idx,buses_idx,reference_bus,BasePointType.FLATSTART)
        #timer.toc('done')

        ## store some information we'll need when iterating on the model object
        _PTDF_dict['PTDFM'] = PTDFM
        _PTDF_dict['buses_idx'] = buses_idx
        _PTDF_dict['branches_idx'] = branches_idx
        m._PTDFs_dict[branches_out_service] = _PTDF_dict

    ## create a pointer on this block to all this PTDF data,
    ## for easy iteration within the solve loop
    block._PTDF_dict = m._PTDFs_dict[branches_out_service]

    ## this expression is specific to each block
    block._PTDF_bus_nw_exprs = \
      [ block.pl[bus] + bus_gs_fixed_shunts[bus] - \
      sum(block.pg[g] for g in gens_by_bus[bus]) for bus in buses]
    block._PTDF_bus_p_loads = bus_p_loads

    ### get the sets for this block, based on the logic above
    PTDF_dict = block._PTDF_dict

    PTDFM = PTDF_dict['PTDFM']
    buses_idx = PTDF_dict['buses_idx']
    branches_idx = PTDF_dict['branches_idx']

    for i,branch_name in enumerate(branches_idx):
        branch = branches[branch_name]
        ptdf_row = {bus : PTDFM[i,j] for j, bus in enumerate(buses_idx)}
        branch['ptdf'] = ptdf_row

    p_max = {k: branches[k]['rating_long_term'] for k in branches_idx}

    ### declare the power flows and their limits
    libbranch.declare_expr_pf(model=block,
                             index_set=branches_in_service,
                             )

    ### declare the branch power flow approximation constraints
    libbranch.declare_eq_branch_power_ptdf_approx(model=block,
                                                  index_set=branches_in_service,
                                                  branches=branches,
                                                  bus_p_loads=bus_p_loads,
                                                  gens_by_bus=gens_by_bus,
                                                  bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                                  abs_ptdf_tol=abs_ptdf_tol,
                                                  rel_ptdf_tol=rel_ptdf_tol
                                                  )
    ### declare the p balance
    libbus.declare_eq_p_balance_ed(model=block,
                                   index_set=buses.keys(),
                                   bus_p_loads=bus_p_loads,
                                   gens_by_bus=gens_by_bus,
                                   bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                   )

    ### declare the real power flow limits
    libbranch.declare_ineq_p_branch_thermal_lbub(model=block,
                                                 index_set=branches_in_service,
                                                 branches=branches,
                                                 p_thermal_limits=p_max,
                                                 approximation_type=ApproximationType.PTDF
                                                 )

def _btheta_dcopf_network_model(md,block,tm,td):
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
    ref_bus = md.data['system']['reference_bus']
    block.va[ref_bus].fix(0.0)

    ref_angle = md.data['system']['reference_bus_angle']
    if ref_angle != 0.0:
        raise ValueError('The BTHETA DCOPF formulation currently only supports'
                         ' a reference bus angle of 0 degrees, but an angle'
                         ' of {} degrees was found.'.format(ref_angle))

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

    md = model.model_data

    if 'reference_bus' in md.data['system'] and md.data['system']['reference_bus'] in model.Buses:
        reference_bus = md.data['system']['reference_bus']
    else:
        reference_bus = list(sorted(m.Buses))[0]
    
    ## for interfacing with the rest of the model code
    def define_pos_neg_load_generate_mismatch_rule(m, b, t):
        if b == reference_bus:
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

    #####################################################
    # load "shedding" can be both positive and negative #
    #####################################################
    model.LoadGenerateMismatch = Var(model.Buses, model.TimePeriods, within=Reals)
    model.posLoadGenerateMismatch = Var(model.Buses, model.TimePeriods, within=NonNegativeReals) # load shedding
    model.negLoadGenerateMismatch = Var(model.Buses, model.TimePeriods, within=NonNegativeReals) # over generation
    
    def define_pos_neg_load_generate_mismatch_rule(m, b, t):
        return m.posLoadGenerateMismatch[b, t] - m.negLoadGenerateMismatch[b, t] == m.LoadGenerateMismatch[b, t]
    model.DefinePosNegLoadGenerateMismatch = Constraint(model.Buses, model.TimePeriods, rule = define_pos_neg_load_generate_mismatch_rule)

    # the following constraints are necessarily, at least in the case of CPLEX 12.4, to prevent
    # the appearance of load generation mismatch component values in the range of *negative* e-5.
    # what these small negative values do is to cause the optimal objective to be a very large negative,
    # due to obviously large penalty values for under or over-generation. JPW would call this a heuristic
    # at this point, but it does seem to work broadly. we tried a single global constraint, across all
    # buses, but that failed to correct the problem, and caused the solve times to explode.
    
    def pos_load_generate_mismatch_tolerance_rule(m, b):
       return sum((m.posLoadGenerateMismatch[b,t] for t in m.TimePeriods)) >= 0.0
    model.PosLoadGenerateMismatchTolerance = Constraint(model.Buses, rule=pos_load_generate_mismatch_tolerance_rule)
    
    def neg_load_generate_mismatch_tolerance_rule(m, b):
       return sum((m.negLoadGenerateMismatch[b,t] for t in m.TimePeriods)) >= 0.0
    model.NegLoadGenerateMismatchTolerance = Constraint(model.Buses, rule=neg_load_generate_mismatch_tolerance_rule)

    def compute_load_mismatch_cost_rule(m, t):
        return m.LoadMismatchPenalty*m.TimePeriodLengthHours*sum(m.posLoadGenerateMismatch[b, t] + m.negLoadGenerateMismatch[b, t] for b in m.Buses) 
    model.LoadMismatchCost = Expression(model.TimePeriods, rule=compute_load_mismatch_cost_rule)

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

    md = model.model_data

    # for transmission network
    model.TransmissionBlock = Block(model.TimePeriods, concrete=True)

    for tm, td in zip(model.TimePeriods, md.data['system']['time_indices']):
        b = model.TransmissionBlock[tm]
        ## this creates a fake bus generator for all the
        ## appropriate injection/withdraws from the unit commitment
        ## model
        b.pg = Expression(model.Buses, rule=_get_pg_expr_rule(tm))
        if reactive_power:
            b.qg = Expression(model.Buses, rule=_get_qg_expr_rule(tm))
        b.gens_by_bus = {bus : [bus] for bus in model.Buses}
        md_t = md.clone_at_timestamp(td)
        network_model_builder(md_t,b,tm,td)

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
def ptdf_power_flow(model, slacks=True):
    model._PTDFs_dict = dict()
    _add_egret_power_flow(model, _ptdf_dcopf_network_model, reactive_power=False, slacks=slacks)

@add_model_attr(component_name, requires = {'data_loader': None,
                                            'power_vars': None,
                                            'non_dispatchable_vars': None,
                                            'storage_service': None,
                                            })
def lazy_ptdf_power_flow(model, slacks=True):
    model._PTDFs_dict = dict()
    _add_egret_power_flow(model, _lazy_ptdf_dcopf_network_model, reactive_power=False, slacks=slacks)

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
