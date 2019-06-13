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

import pyomo.environ as pe
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen

from egret.model_library.defn import CoordinateType, ApproximationType
from math import pi

component_name = 'power_balance'

def _btheta_dcopf_network_model(md,block):
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    loads = dict(md.elements(element_type='load'))
    shunts = dict(md.elements(element_type='shunt'))

    bus_attrs = md.attributes(element_type='bus')
    branch_attrs = md.attributes(element_type='branch')

    inlet_branches_by_bus, outlet_branches_by_bus = \
        tx_utils.inlet_outlet_branches_by_bus(branches, buses)

    ## this is not the "real" gens by bus, but the
    ## index of net injections from the UC model
    gens_by_bus = block.gens_by_bus

    ### declare (and fix) the loads at the buses
    bus_p_loads, _ = tx_utils.dict_of_bus_loads(buses, loads)

    libbus.declare_var_pl(block, bus_attrs['names'], initialize=bus_p_loads)
    block.pl.fix()

    ### declare the fixed shunts at the buses
    _, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    ### declare the polar voltages
    va_bounds = {k: (-pi, pi) for k in bus_attrs['names']}
    libbus.declare_var_va(block, bus_attrs['names'], initialize=None,
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

    p_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
    p_lbub = {k: (-p_max[k],p_max[k]) for k in branches.keys()}
    pf_bounds = p_lbub

    libbranch.declare_var_pf(model=block,
                             index_set=branch_attrs['names'],
                             initialize=None,
                             bounds=pf_bounds
                             )

    ### declare the branch power flow approximation constraints
    libbranch.declare_eq_branch_power_btheta_approx(model=block,
                                                    index_set=branch_attrs['names'],
                                                    branches=branches
                                                    )

    ### declare the p balance
    libbus.declare_eq_p_balance_dc_approx(model=block,
                                          index_set=bus_attrs['names'],
                                          bus_p_loads=bus_p_loads,
                                          gens_by_bus=gens_by_bus,
                                          bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                          inlet_branches_by_bus=inlet_branches_by_bus,
                                          outlet_branches_by_bus=outlet_branches_by_bus,
                                          approximation_type=ApproximationType.BTHETA
                                          )

    ### declare the real power flow limits
    libbranch.declare_ineq_p_branch_thermal_lbub(model=block,
                                                 index_set=branch_attrs['names'],
                                                 branches=branches,
                                                 p_thermal_limits=p_max,
                                                 approximation_type=ApproximationType.BTHETA
                                                 )

    ### declare angle difference limits on interconnected buses
    libbranch.declare_ineq_angle_diff_branch_lbub(model=block,
                                                  index_set=branch_attrs['names'],
                                                  branches=branches,
                                                  coordinate_type=CoordinateType.POLAR
                                                  )

    return block


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

def _add_blank_load_mismatch(model):
    model.LoadGenerateMismatch = Param(model.Buses, model.TimePeriods, default=0.)
    model.posLoadGenerateMismatch = Param(model.Buses, model.TimePeriods, default=0.)
    model.negLoadGenerateMismatch = Param(model.Buses, model.TimePeriods, default=0.)

def _add_blank_q_load_mismatch(model):
    model.LoadGenerateMismatchReactive = Param(model.Buses, model.TimePeriods, default=0.)
    model.posLoadGenerateMismatchReactive = Param(model.Buses, model.TimePeriods, default=0.)
    model.negLoadGenerateMismatchReactive = Param(model.Buses, model.TimePeriods, default=0.)

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

    if slacks:
        _add_load_mismatch(model)
    else:
        _add_blank_load_mismatch(model)

    if reactive_power:
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
        network_model_builder(md_t,b)


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
