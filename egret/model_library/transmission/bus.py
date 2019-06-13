#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module contains the declarations for the modeling components
typically used for buses (including loads and shunts)
"""
import pyomo.environ as pe
import egret.model_library.decl as decl
from egret.model_library.defn import FlowType, CoordinateType, ApproximationType


def declare_var_vr(model, index_set, **kwargs):
    """
    Create variable for the real component of the voltage at a bus
    """
    decl.declare_var('vr', model=model, index_set=index_set, **kwargs)


def declare_var_vj(model, index_set, **kwargs):
    """
    Create variable for the imaginary component of the voltage at a bus
    """
    decl.declare_var('vj', model=model, index_set=index_set, **kwargs)


def declare_var_vm(model, index_set, **kwargs):
    """
    Create variable for the voltage magnitude of the voltage at a bus
    """
    decl.declare_var('vm', model=model, index_set=index_set, **kwargs)


def declare_var_va(model, index_set, **kwargs):
    """
    Create variable for the phase angle of the voltage at a bus
    """
    decl.declare_var('va', model=model, index_set=index_set, **kwargs)


def declare_var_ir_aggregation_at_bus(model, index_set, **kwargs):
    """
    Create a variable for the aggregated real current at a bus
    """
    decl.declare_var('ir_aggregation_at_bus', model=model, index_set=index_set, **kwargs)


def declare_var_ij_aggregation_at_bus(model, index_set, **kwargs):
    """
    Create a variable for the aggregated imaginary current at a bus
    """
    decl.declare_var('ij_aggregation_at_bus', model=model, index_set=index_set, **kwargs)


def declare_var_pl(model, index_set, **kwargs):
    """
    Create variable for the real power load at a bus
    """
    decl.declare_var('pl', model=model, index_set=index_set, **kwargs)


def declare_var_ql(model, index_set, **kwargs):
    """
    Create variable for the reactive power load at a bus
    """
    decl.declare_var('ql', model=model, index_set=index_set, **kwargs)


def declare_expr_shunt_power_at_bus(model, index_set, shunt_attrs,
                                    coordinate_type=CoordinateType.POLAR):
    """
    Create the expression for the shunt power at the bus
    """
    m = model
    expr_set = decl.declare_set('_expr_shunt_at_bus_set', model, index_set)

    m.shunt_p = pe.Expression(expr_set, initialize=0.0)
    m.shunt_q = pe.Expression(expr_set, initialize=0.0)

    if coordinate_type == CoordinateType.POLAR:
        for bus_name in expr_set:
            if bus_name in shunt_attrs['bus']:
                vmsq = m.vm[bus_name]**2
                m.shunt_p[bus_name] = shunt_attrs['gs'][bus_name]*vmsq
                m.shunt_q[bus_name] = -shunt_attrs['bs'][bus_name]*vmsq
    elif coordinate_type == CoordinateType.RECTANGULAR:
        for bus_name in expr_set:
            if bus_name in shunt_attrs['bus']:
                vmsq = m.vr[bus_name]**2 + m.vj[bus_name]**2
                m.shunt_p[bus_name] = shunt_attrs['gs'][bus_name]*vmsq
                m.shunt_q[bus_name] = -shunt_attrs['bs'][bus_name]*vmsq


def declare_eq_i_aggregation_at_bus(model, index_set,
                                    bus_bs_fixed_shunts, bus_gs_fixed_shunts,
                                    inlet_branches_by_bus, outlet_branches_by_bus):
    """
    Create the equality constraints for the aggregated real and imaginary
    currents at the bus
    """
    m = model
    con_set = decl.declare_set('_con_eq_i_aggregation_at_bus_set', model, index_set)

    m.eq_ir_aggregation_at_bus = pe.Constraint(con_set)
    m.eq_ij_aggregation_at_bus = pe.Constraint(con_set)

    for bus_name in con_set:
        ir_expr = sum([m.ifr[branch_name] for branch_name in outlet_branches_by_bus[bus_name]])
        ir_expr += sum([m.itr[branch_name] for branch_name in inlet_branches_by_bus[bus_name]])
        ij_expr = sum([m.ifj[branch_name] for branch_name in outlet_branches_by_bus[bus_name]])
        ij_expr += sum([m.itj[branch_name] for branch_name in inlet_branches_by_bus[bus_name]])

        if bus_bs_fixed_shunts[bus_name] != 0.0:
            ir_expr -= bus_bs_fixed_shunts[bus_name] * m.vj[bus_name]
            ij_expr += bus_bs_fixed_shunts[bus_name] * m.vr[bus_name]
        if bus_gs_fixed_shunts[bus_name] != 0.0:
            ir_expr += bus_gs_fixed_shunts[bus_name] * m.vr[bus_name]
            ij_expr += bus_gs_fixed_shunts[bus_name] * m.vj[bus_name]

        ir_expr -= m.ir_aggregation_at_bus[bus_name]
        ij_expr -= m.ij_aggregation_at_bus[bus_name]

        m.eq_ir_aggregation_at_bus[bus_name] = ir_expr == 0
        m.eq_ij_aggregation_at_bus[bus_name] = ij_expr == 0


def declare_eq_p_balance_ed(model, index_set, bus_p_loads, gens_by_bus, bus_gs_fixed_shunts, **rhs_kwargs):
    """
    Create the equality constraints for the real power balance
    at a bus using the variables for real power flows, respectively.

    NOTE: Equation build orientates constants to the RHS in order to compute the correct dual variable sign
    """
    m = model

    p_expr = sum(m.pg[gen_name] for bus_name in index_set for gen_name in gens_by_bus[bus_name])
    p_expr -= sum(m.pl[bus_name] for bus_name in index_set if bus_p_loads[bus_name] is not None)
    p_expr -= sum(bus_gs_fixed_shunts[bus_name] for bus_name in index_set if bus_gs_fixed_shunts[bus_name] != 0.0)

    if rhs_kwargs:
        for idx,val in rhs_kwargs.items():
            if idx == 'include_feasibility_slack_pos':
                p_expr -= eval("m." + val)
            if idx == 'include_feasibility_slack_neg':
                p_expr += eval("m." + val)
            if idx == 'include_losses':
                p_expr -= sum(m.pfl[branch_name] for branch_name in val)

    m.eq_p_balance = pe.Constraint(expr = p_expr == 0.0)


def declare_eq_p_balance_dc_approx(model, index_set,
                                   bus_p_loads,
                                   gens_by_bus,
                                   bus_gs_fixed_shunts,
                                   inlet_branches_by_bus, outlet_branches_by_bus,
                                   approximation_type=ApproximationType.BTHETA,
                                   **rhs_kwargs):
    """
    Create the equality constraints for the real power balance
    at a bus using the variables for real power flows, respectively.

    NOTE: Equation build orientates constants to the RHS in order to compute the correct dual variable sign
    """
    m = model
    con_set = decl.declare_set('_con_eq_p_balance', model, index_set)

    m.eq_p_balance = pe.Constraint(con_set)

    for bus_name in con_set:
        if approximation_type == ApproximationType.BTHETA:
            p_expr = -sum([m.pf[branch_name] for branch_name in outlet_branches_by_bus[bus_name]])
            p_expr += sum([m.pf[branch_name] for branch_name in inlet_branches_by_bus[bus_name]])
        elif approximation_type == ApproximationType.BTHETA_LOSSES:
            p_expr = -0.5*sum([m.pfl[branch_name] for branch_name in inlet_branches_by_bus[bus_name]])
            p_expr -= 0.5*sum([m.pfl[branch_name] for branch_name in outlet_branches_by_bus[bus_name]])
            p_expr -= sum([m.pf[branch_name] for branch_name in outlet_branches_by_bus[bus_name]])
            p_expr += sum([m.pf[branch_name] for branch_name in inlet_branches_by_bus[bus_name]])

        if bus_gs_fixed_shunts[bus_name] != 0.0:
            p_expr -= bus_gs_fixed_shunts[bus_name]

        if bus_p_loads[bus_name] != 0.0: # only applies to fixed loads, otherwise may cause an error
            p_expr -= m.pl[bus_name]

        if rhs_kwargs:
            for idx, val in rhs_kwargs.items():
                if idx == 'include_feasibility_slack_pos':
                    p_expr -= eval("m." + val)[bus_name]
                if idx == 'include_feasibility_slack_neg':
                    p_expr += eval("m." + val)[bus_name]

        for gen_name in gens_by_bus[bus_name]:
            p_expr += m.pg[gen_name]

        m.eq_p_balance[bus_name] = \
            p_expr == 0.0


def declare_eq_p_balance(model, index_set,
                         bus_p_loads,
                         gens_by_bus,
                         bus_gs_fixed_shunts,
                         inlet_branches_by_bus, outlet_branches_by_bus,
                         coordinate_type=CoordinateType.RECTANGULAR,
                         **rhs_kwargs):
    """
    Create the equality constraints for the real power balance
    at a bus using the variables for real power flows, respectively.

    NOTE: Equation build orientates constants to the RHS in order to compute the correct dual variable sign
    """

    m = model
    con_set = decl.declare_set('_con_eq_p_balance', model, index_set)

    m.eq_p_balance = pe.Constraint(con_set)

    for bus_name in con_set:
        p_expr = -sum([m.pf[branch_name] for branch_name in outlet_branches_by_bus[bus_name]])
        p_expr -= sum([m.pt[branch_name] for branch_name in inlet_branches_by_bus[bus_name]])

        if bus_gs_fixed_shunts[bus_name] != 0.0:
            if coordinate_type == CoordinateType.RECTANGULAR:
                vmsq = m.vr[bus_name] ** 2 + m.vj[bus_name] ** 2
            elif coordinate_type == CoordinateType.POLAR:
                vmsq = m.vm[bus_name] ** 2
            p_expr -= bus_gs_fixed_shunts[bus_name] * vmsq

        if bus_p_loads[bus_name] != 0.0: # only applies to fixed loads, otherwise may cause an error
            p_expr -= m.pl[bus_name]

        if rhs_kwargs:
            for idx, val in rhs_kwargs.items():
                if idx == 'include_feasibility_slack_pos':
                    p_expr -= eval("m." + val)[bus_name]
                if idx == 'include_feasibility_slack_neg':
                    p_expr += eval("m." + val)[bus_name]

        for gen_name in gens_by_bus[bus_name]:
            p_expr += m.pg[gen_name]

        m.eq_p_balance[bus_name] = \
            p_expr == 0.0


def declare_eq_p_balance_with_i_aggregation(model, index_set,
                                            bus_p_loads,
                                            gens_by_bus,
                                            **rhs_kwargs):
    """
    Create the equality constraints for the real power balance
    at a bus using the variables for real power flows, respectively.

    NOTE: Equation build orientates constants to the RHS in order to compute the correct dual variable sign
    """
    m = model
    con_set = decl.declare_set('_con_eq_p_balance', model, index_set)

    m.eq_p_balance = pe.Constraint(con_set)

    for bus_name in con_set:
        p_expr = -m.vr[bus_name] * m.ir_aggregation_at_bus[bus_name] + \
                 -m.vj[bus_name] * m.ij_aggregation_at_bus[bus_name]

        if bus_p_loads[bus_name] != 0.0: # only applies to fixed loads, otherwise may cause an error
            p_expr -= m.pl[bus_name]

        if rhs_kwargs:
            for idx, val in rhs_kwargs.items():
                if idx == 'include_feasibility_slack_pos':
                    p_expr -= eval("m." + val)[bus_name]
                if idx == 'include_feasibility_slack_neg':
                    p_expr += eval("m." + val)[bus_name]

        for gen_name in gens_by_bus[bus_name]:
            p_expr += m.pg[gen_name]

        m.eq_p_balance[bus_name] = \
            p_expr == 0.0


def declare_eq_q_balance(model, index_set,
                         bus_q_loads,
                         gens_by_bus,
                         bus_bs_fixed_shunts,
                         inlet_branches_by_bus, outlet_branches_by_bus,
                         coordinate_type=CoordinateType.POLAR,
                         **rhs_kwargs):
    """
    Create the equality constraints for the reactive power balance
    at a bus using the variables for reactive power flows, respectively.

    NOTE: Equation build orientates constants to the RHS in order to compute the correct dual variable sign
    """
    m = model
    con_set = decl.declare_set('_con_eq_q_balance', model, index_set)

    m.eq_q_balance = pe.Constraint(con_set)

    for bus_name in con_set:
        q_expr = -sum([m.qf[branch_name] for branch_name in outlet_branches_by_bus[bus_name]])
        q_expr -= sum([m.qt[branch_name] for branch_name in inlet_branches_by_bus[bus_name]])

        if bus_bs_fixed_shunts[bus_name] != 0.0:
            if coordinate_type == CoordinateType.RECTANGULAR:
                vmsq = m.vr[bus_name] ** 2 + m.vj[bus_name] ** 2
            elif coordinate_type == CoordinateType.POLAR:
                vmsq = m.vm[bus_name] ** 2
            q_expr += bus_bs_fixed_shunts[bus_name] * vmsq

        if bus_q_loads[bus_name] != 0.0: # only applies to fixed loads, otherwise may cause an error
            q_expr -= m.ql[bus_name]

        if rhs_kwargs:
            for idx, val in rhs_kwargs.items():
                if idx == 'include_feasibility_slack_pos':
                    q_expr -= eval("m." + val)[bus_name]
                if idx == 'include_feasibility_slack_neg':
                    q_expr += eval("m." + val)[bus_name]

        for gen_name in gens_by_bus[bus_name]:
            q_expr += m.qg[gen_name]

        m.eq_q_balance[bus_name] = \
            q_expr == 0.0


def declare_eq_q_balance_with_i_aggregation(model, index_set,
                                            bus_q_loads,
                                            gens_by_bus,
                                            **rhs_kwargs):
    """
    Create the equality constraints for the reactive power balance
    at a bus using the variables for reactive power flows, respectively.

    NOTE: Equation build orientates constants to the RHS in order to compute the correct dual variable sign
    """
    m = model
    con_set = decl.declare_set('_con_eq_q_balance', model, index_set)

    m.eq_q_balance = pe.Constraint(con_set)

    for bus_name in con_set:
        q_expr = m.vr[bus_name] * m.ij_aggregation_at_bus[bus_name] + \
                 -m.vj[bus_name] * m.ir_aggregation_at_bus[bus_name]

        if bus_q_loads[bus_name] != 0.0: # only applies to fixed loads, otherwise may cause an error
            q_expr -= m.ql[bus_name]

        if rhs_kwargs:
            for idx, val in rhs_kwargs.items():
                if idx == 'include_feasibility_slack_pos':
                    q_expr -= eval("m." + val)[bus_name]
                if idx == 'include_feasibility_slack_neg':
                    q_expr += eval("m." + val)[bus_name]

        for gen_name in gens_by_bus[bus_name]:
            q_expr += m.qg[gen_name]

        m.eq_q_balance[bus_name] = \
            q_expr == 0.0


def declare_ineq_vm_bus_lbub(model, index_set, buses, coordinate_type=CoordinateType.POLAR):
    """
    Create the inequalities for the voltage magnitudes from the
    voltage variables
    """
    m = model
    con_set = decl.declare_set('_con_ineq_vm_bus_lbub',
                               model=model, index_set=index_set)

    m.ineq_vm_bus_lb = pe.Constraint(con_set)
    m.ineq_vm_bus_ub = pe.Constraint(con_set)

    if coordinate_type == CoordinateType.POLAR:
        for bus_name in con_set:
            m.ineq_vm_bus_lb[bus_name] = \
                buses[bus_name]['v_min'] <= m.vm[bus_name]
            m.ineq_vm_bus_ub[bus_name] = \
                m.vm[bus_name] <= buses[bus_name]['v_max']
    elif coordinate_type == CoordinateType.RECTANGULAR:
        for bus_name in con_set:
            m.ineq_vm_bus_lb[bus_name] = \
                buses[bus_name]['v_min']**2 <= m.vr[bus_name]**2 + m.vj[bus_name]**2
            m.ineq_vm_bus_ub[bus_name] = \
                m.vr[bus_name]**2 + m.vj[bus_name]**2 <= buses[bus_name]['v_max']**2
