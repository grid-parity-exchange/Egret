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
typically used for transmission lines
"""
import math
import pyomo.environ as pe
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.decl as decl
from egret.model_library.defn import FlowType, CoordinateType, ApproximationType, RelaxationType
from egret.data.data_utils import zip_items
from pyomo.core.util import quicksum
from pyomo.core.expr.numeric_expr import LinearExpression
from collections import OrderedDict
from pyomo.contrib.fbbt.fbbt import fbbt
try:
    import coramin
    coramin_available = True
except ImportError:
    coramin_available = False


def declare_var_dva(model, index_set, **kwargs):
    """
    Create variable or the angle difference between interconnected bus pairs
    """
    decl.declare_var('dva', model=model, index_set=index_set, **kwargs)


def declare_var_pfl(model, index_set, **kwargs):
    """
    Create variable for the real part of the power loss in the transmission
    line
    """
    decl.declare_var('pfl', model=model, index_set=index_set, **kwargs)


def declare_var_pf(model, index_set, **kwargs):
    """
    Create variable for the real part of the power flow in the "from"
    end of the transmission line
    """
    decl.declare_var('pf', model=model, index_set=index_set, **kwargs)


def declare_expr_pf(model, index_set, **kwargs):
    """
    Create expression for the real part of the power flow in the "from"
    end of the transmission line
    """
    decl.declare_expr('pf', model=model, index_set=index_set, **kwargs)


def declare_var_pf_slack_pos(model, index_set, **kwargs):
    """
    Create the positive slack variable for the real part of power flow
    in the "from" end of the transmission line
    """
    decl.declare_var('pf_slack_pos', model=model, index_set=index_set, **kwargs)


def declare_var_pf_slack_neg(model, index_set, **kwargs):
    """
    Create the negative slack variable for the real part of power flow
    in the "from" end of the transmission line
    """
    decl.declare_var('pf_slack_neg', model=model, index_set=index_set, **kwargs)


def declare_var_pfi(model, index_set, **kwargs):
    """
    Create variable for the real part of the power flow through an interface
    """
    decl.declare_var('pfi', model=model, index_set=index_set, **kwargs)


def declare_expr_pfi(model, index_set, **kwargs):
    """
    Create expression for the real part of the power flow through an interface
    """
    decl.declare_expr('pfi', model=model, index_set=index_set, **kwargs)


def declare_var_pfi_slack_pos(model, index_set, **kwargs):
    """
    Create the positive slack variable for the real part of power flow
    through an interface
    """
    decl.declare_var('pfi_slack_pos', model=model, index_set=index_set, **kwargs)


def declare_var_pfi_slack_neg(model, index_set, **kwargs):
    """
    Create the negative slack variable for the real part of power flow
    through an interface
    """
    decl.declare_var('pfi_slack_neg', model=model, index_set=index_set, **kwargs)


def declare_var_qf(model, index_set, **kwargs):
    """
    Create variable for the imaginary part of the power flow in the "from"
    end of the transmission line
    """
    decl.declare_var('qf', model=model, index_set=index_set, **kwargs)


def declare_var_pt(model, index_set, **kwargs):
    """
    Create variable for the real part of the power flow in the "to"
    end of the transmission line
    """
    decl.declare_var('pt', model=model, index_set=index_set, **kwargs)


def declare_var_qt(model, index_set, **kwargs):
    """
    Create variable for the imaginary part of the power flow in the "to"
    end of the transmission line
    """
    decl.declare_var('qt', model=model, index_set=index_set, **kwargs)


def declare_var_ifr(model, index_set, **kwargs):
    """
    Create variable for the real part of the current flow in the "from"
    end of the transmission line
    """
    decl.declare_var('ifr', model=model, index_set=index_set, **kwargs)


def declare_var_ifj(model, index_set, **kwargs):
    """
    Create variable for the imaginary part of the current flow in the "from"
    end of the transmission line
    """
    decl.declare_var('ifj', model=model, index_set=index_set, **kwargs)


def declare_var_itr(model, index_set, **kwargs):
    """
    Create variable for the real part of the current flow in the "to"
    end of the transmission line
    """
    decl.declare_var('itr', model=model, index_set=index_set, **kwargs)


def declare_var_itj(model, index_set, **kwargs):
    """
    Create variable for the imaginary part of the current flow in the "to"
    end of the transmission line
    """
    decl.declare_var('itj', model=model, index_set=index_set, **kwargs)


def declare_eq_branch_dva(model, index_set, branches):
    """
    Create the equality constraints for the angle difference
    in the branch

    dva = va[from_bus] - va[to_bus] + transformer_phase_shift
    """
    m = model

    con_set = decl.declare_set("_con_eq_branch_dva_set", model, index_set)

    m.eq_dva_branch = pe.Constraint(con_set)
    for branch_name in con_set:
        branch = branches[branch_name]

        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        shift = 0.0
        if branch['branch_type'] == 'transformer':
            shift = math.radians(branch['transformer_phase_shift'])

        m.eq_dva_branch[branch_name] = \
            m.dva[branch_name] == \
            m.va[from_bus] - m.va[to_bus] + shift


def declare_eq_delta_va(model, index_set):
    """
    Create the equality constraints for the angle difference
    in the branch

    dva = va[from_bus] - va[to-bus]
    """
    m = model
    con_set = decl.declare_set("_con_eq_delta_va_set", model, index_set)
    m.eq_delta_va = pe.Constraint(con_set)

    for from_bus, to_bus in con_set:
        m.eq_delta_va[(from_bus, to_bus)] = m.dva[(from_bus, to_bus)] == m.va[from_bus] - m.va[to_bus]


def declare_expr_c(model, index_set, coordinate_type=CoordinateType.POLAR):
    """
    Create expression for the nonlinear, nonconvex term based on cosine
    of the phase angle difference (polar) or bilinear voltages (rectangular)
    """
    m = model
    expr_set = decl.declare_set('_expr_c', model, index_set)
    m.c = pe.Expression(expr_set)

    if coordinate_type == CoordinateType.RECTANGULAR:
        for from_bus, to_bus in expr_set:
            m.c[(from_bus,to_bus)] = m.vr[from_bus]*m.vr[to_bus] + m.vj[from_bus]*m.vj[to_bus]
    elif coordinate_type == CoordinateType.POLAR:
        for from_bus, to_bus in expr_set:
            m.c[(from_bus,to_bus)] = m.vm[from_bus]*m.vm[to_bus]*pe.cos(m.va[from_bus]-m.va[to_bus])


def declare_expr_s(model, index_set, coordinate_type=CoordinateType.POLAR):
    """
    Create expression for the nonlinear, nonconvex term based on cosine
    of the phase angle difference (polar) or bilinear voltages (rectangular)
    """
    m = model
    expr_set = decl.declare_set('_expr_s', model, index_set)
    m.s = pe.Expression(expr_set)

    if coordinate_type == CoordinateType.RECTANGULAR:
        for from_bus, to_bus in expr_set:
            m.s[(from_bus,to_bus)] = m.vj[from_bus]*m.vr[to_bus] - m.vr[from_bus]*m.vj[to_bus]
    elif coordinate_type == CoordinateType.POLAR:
        for from_bus, to_bus in expr_set:
            m.s[(from_bus,to_bus)] = m.vm[from_bus]*m.vm[to_bus]*pe.sin(m.va[from_bus]-m.va[to_bus])


def declare_var_c(model, index_set, **kwargs):
    """
    Create an auxiliary variable for vf * vt * cos(theta_f - theta_t)
    """
    decl.declare_var('c', model=model, index_set=index_set, **kwargs)


def declare_var_s(model, index_set, **kwargs):
    """
    Create an auxiliary variable for vf * vt * sin(theta_f - theta_t)
    """
    decl.declare_var('s', model=model, index_set=index_set, **kwargs)


def declare_eq_c(model, index_set, coordinate_type=CoordinateType.POLAR):
    """
    Create a constraint relating c to the voltages
    """
    m = model
    con_set = decl.declare_set('_con_eq_c', model, index_set)
    m.eq_c = pe.Constraint(con_set)

    if coordinate_type == CoordinateType.POLAR:
        for from_bus, to_bus in con_set:
            m.eq_c[(from_bus, to_bus)] = (m.c[(from_bus, to_bus)] ==
                                          m.vm[from_bus] * m.vm[to_bus] * pe.cos(m.dva[(from_bus, to_bus)]))
    elif coordinate_type == CoordinateType.RECTANGULAR:
        for from_bus, to_bus in con_set:
            m.eq_c[(from_bus, to_bus)] = (m.c[(from_bus, to_bus)] ==
                                          m.vr[from_bus] * m.vr[to_bus] + m.vj[from_bus] * m.vj[to_bus])
    else:
        raise ValueError('unexpected coordinate_type: {0}'.format(str(coordinate_type)))


def declare_eq_s(model, index_set, coordinate_type=CoordinateType.POLAR):
    """
    Create a constraint relating s to the voltages
    """
    m = model
    con_set = decl.declare_set('_con_eq_s', model, index_set)
    m.eq_s = pe.Constraint(con_set)

    if coordinate_type == CoordinateType.POLAR:
        for from_bus, to_bus in con_set:
            m.eq_s[(from_bus, to_bus)] = (m.s[(from_bus, to_bus)] ==
                                          m.vm[from_bus] * m.vm[to_bus] * pe.sin(m.dva[(from_bus, to_bus)]))
    elif coordinate_type == CoordinateType.RECTANGULAR:
        for from_bus, to_bus in con_set:
            m.eq_s[(from_bus, to_bus)] = (m.s[(from_bus, to_bus)] ==
                                          m.vj[from_bus] * m.vr[to_bus] - m.vr[from_bus] * m.vj[to_bus])
    else:
        raise ValueError('unexpected coordinate_type: {0}'.format(str(coordinate_type)))


def declare_eq_branch_current(model, index_set, branches, coordinate_type=CoordinateType.RECTANGULAR):
    """
    Create the equality constraints for the real and imaginary current
    in the branch
    """
    assert(coordinate_type != CoordinateType.POLAR
           and "Branch current in polar coordinates not implemented.")

    m = model
    con_set = decl.declare_set("_con_eq_branch_current_set", model, index_set)

    m.eq_ifr_branch = pe.Constraint(con_set)
    m.eq_ifj_branch = pe.Constraint(con_set)
    m.eq_itr_branch = pe.Constraint(con_set)
    m.eq_itj_branch = pe.Constraint(con_set)
    for branch_name in con_set:
        branch = branches[branch_name]

        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        g = tx_calc.calculate_conductance(branch)
        b = tx_calc.calculate_susceptance(branch)
        bc = branch['charging_susceptance']
        tau = 1.0
        shift = 0.0

        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])

        g11 = g / tau**2
        g12 = (g * math.cos(shift) - b * math.sin(shift)) / tau
        g21 = (g * math.cos(shift) + b * math.sin(shift)) / tau
        g22 = g

        b11 = (b + bc / 2) / tau**2
        b12 = (b * math.cos(shift) + g*math.sin(shift)) / tau
        b21 = (b * math.cos(shift) - g*math.sin(shift)) / tau
        b22 = b + bc / 2

        m.eq_ifr_branch[branch_name] = \
            m.ifr[branch_name] == \
            g11 * m.vr[from_bus] - g12 * m.vr[to_bus] - (b11 * m.vj[from_bus] - b12 * m.vj[to_bus])

        m.eq_ifj_branch[branch_name] = \
            m.ifj[branch_name] == \
            g11 * m.vj[from_bus] - g12 * m.vj[to_bus] + (b11 * m.vr[from_bus] - b12 * m.vr[to_bus])

        m.eq_itr_branch[branch_name] = \
            m.itr[branch_name] == \
            -(g21 * m.vr[from_bus] - g22 * m.vr[to_bus] - (b21 * m.vj[from_bus] - b22 * m.vj[to_bus]))

        m.eq_itj_branch[branch_name] = \
            m.itj[branch_name] == \
            -(g21 * m.vj[from_bus] - g22 * m.vj[to_bus] + (b21 * m.vr[from_bus] - b22 * m.vr[to_bus]))


def declare_eq_branch_power(model, index_set, branches):
    """
    Create the equality constraints for the real and reactive power
    in the branch
    """
    m = model
    con_set = decl.declare_set("_con_eq_branch_power_set", model, index_set)

    m.eq_pf_branch = pe.Constraint(con_set)
    m.eq_pt_branch = pe.Constraint(con_set)
    m.eq_qf_branch = pe.Constraint(con_set)
    m.eq_qt_branch = pe.Constraint(con_set)
    for branch_name in con_set:
        branch = branches[branch_name]

        from_bus = branch['from_bus']
        to_bus = branch['to_bus']
        vmsq_from_bus = m.vmsq[from_bus]
        vmsq_to_bus = m.vmsq[to_bus]

        g = tx_calc.calculate_conductance(branch)
        b = tx_calc.calculate_susceptance(branch)
        bc = branch['charging_susceptance']
        tau = 1.0
        shift = 0.0

        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])

        g11 = g / tau ** 2
        g12 = g * math.cos(shift) / tau
        g21 = g * math.sin(shift) / tau
        g22 = g

        b11 = (b + bc / 2) / tau ** 2
        b12 = b * math.cos(shift) / tau
        b21 = b * math.sin(shift) / tau
        b22 = b + bc / 2

        m.eq_pf_branch[branch_name] = \
            m.pf[branch_name] == \
            g11 * vmsq_from_bus - \
            (g12 * m.c[(from_bus,to_bus)] +
             g21 * m.s[(from_bus,to_bus)] +
             b12 * m.s[(from_bus,to_bus)] -
             b21 * m.c[(from_bus,to_bus)])

        m.eq_pt_branch[branch_name] = \
            m.pt[branch_name] == \
            g22 * vmsq_to_bus - \
            (g12 * m.c[(from_bus,to_bus)] +
             g21 * m.s[(from_bus,to_bus)] -
             b12 * m.s[(from_bus,to_bus)] +
             b21 * m.c[(from_bus,to_bus)])

        m.eq_qf_branch[branch_name] = \
            m.qf[branch_name] == \
            -b11 * vmsq_from_bus + \
            (b12 * m.c[(from_bus,to_bus)] +
             b21 * m.s[(from_bus,to_bus)] -
             g12 * m.s[(from_bus,to_bus)] +
             g21 * m.c[(from_bus,to_bus)])

        m.eq_qt_branch[branch_name] = \
            m.qt[branch_name] == \
            -b22 * vmsq_to_bus + \
            (b12 * m.c[(from_bus,to_bus)] +
             b21 * m.s[(from_bus,to_bus)] +
             g12 * m.s[(from_bus,to_bus)] -
             g21 * m.c[(from_bus,to_bus)])


def declare_ineq_soc(model, index_set, use_outer_approximation=False):
    """
    create the constraint for the second order cone
    """
    m = model
    if not use_outer_approximation:
        con_set = decl.declare_set("_con_ineq_soc", model, index_set)
        m.ineq_soc = pe.Constraint(con_set)
        for from_bus, to_bus in con_set:
            m.ineq_soc[(from_bus, to_bus)] = m.c[from_bus, to_bus] ** 2 + m.s[from_bus, to_bus] ** 2 <= m.vmsq[from_bus] * m.vmsq[to_bus]
    else:
        if not coramin_available:
            raise ImportError('Cannot create SOC relaxation with outer approximation unless coramin is available.')
        """
        in order to use outer approximation, we have to reformulate 
        
        c**2 + s**2 <= vmsq[from_bus] * vmsq[to_bus]
        
        to
        
        (c**2 + s**2 + z1**2) ** 0.5 <= z2
        z1 = 0.5 * (vmsq[from_bus] - vmsq[to_bus])
        z2 = 0.5 * (vmsq[from_bus] + vmsq[to_bus]) 
        """
        con_set = decl.declare_set("_con_ineq_soc", model, index_set)
        decl.declare_var('_z1', model=model, index_set=con_set)
        decl.declare_var('_z2', model=model, index_set=con_set)
        m._eq_z1 = pe.Constraint(con_set)
        m._eq_z2 = pe.Constraint(con_set)
        m.ineq_soc_OA = coramin.relaxations.MultivariateRelaxation(con_set)
        for from_bus, to_bus in con_set:
            m._eq_z1[from_bus, to_bus] = m._z1[from_bus, to_bus] == 0.5 * (m.vmsq[from_bus] - m.vmsq[to_bus])
            m._eq_z2[from_bus, to_bus] = m._z2[from_bus, to_bus] == 0.5 * (m.vmsq[from_bus] + m.vmsq[to_bus])
            fbbt(m._eq_z1[from_bus, to_bus])
            fbbt(m._eq_z2[from_bus, to_bus])
            m.ineq_soc_OA[from_bus, to_bus].build(aux_var=m._z2[from_bus, to_bus],
                                                  shape=coramin.utils.FunctionShape.CONVEX,
                                                  f_x_expr=(m.c[from_bus, to_bus]**2 +
                                                            m.s[from_bus, to_bus]**2 +
                                                            m._z1[from_bus, to_bus]**2)**0.5)


def declare_eq_branch_power_btheta_approx(model, index_set, branches, approximation_type=ApproximationType.BTHETA):
    """
    Create the equality constraints for power (from BTHETA approximation)
    in the branch
    """
    m = model

    con_set = decl.declare_set("_con_eq_branch_power_btheta_approx_set", model, index_set)


    m.eq_pf_branch = pe.Constraint(con_set)
    for branch_name in con_set:
        branch = branches[branch_name]

        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        tau = 1.0
        shift = 0.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])

        if approximation_type == ApproximationType.BTHETA:
            x = branch['reactance']
            b = -1/(tau*x)
        elif approximation_type == ApproximationType.BTHETA_LOSSES:
            b = tx_calc.calculate_susceptance(branch)/tau

        m.eq_pf_branch[branch_name] = \
            m.pf[branch_name] == \
            b * (m.va[from_bus] - m.va[to_bus] + shift)


def declare_eq_branch_loss_btheta_approx(model, index_set, branches, relaxation_type = RelaxationType.NONE):
    """
    Create the equality constraints for losses (from BTHETA approximation)
    in the branch
    """
    m = model

    con_set = decl.declare_set("_con_eq_branch_loss_btheta_approx_set", model, index_set)

    m.eq_pfl_branch = pe.Constraint(con_set)
    for branch_name in con_set:
        branch = branches[branch_name]

        tau = 1.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
        g = tx_calc.calculate_conductance(branch)/tau

        if relaxation_type == RelaxationType.NONE:
            m.eq_pfl_branch[branch_name] = \
                m.pfl[branch_name] == \
                g * (m.dva[branch_name])**2
        elif relaxation_type == RelaxationType.SOC:
            m.eq_pfl_branch[branch_name] = \
                m.pfl[branch_name] >= \
                g * (m.dva[branch_name])**2


def declare_eq_interface_power_btheta_approx(model, index_set, interfaces):
    """
    Create the equality constraints for interface real power flow
    """
    m = model
    con_set = decl.declare_set("_con_eq_interface_power_btheta_approx", model, index_set)

    m.eq_pf_interface = pe.Constraint(con_set)
    for interface_name in con_set:
        interface = interfaces[interface_name]

        expr = 0.
        for line, orientation in zip(interface['lines'], interface['line_orientation']):
            ### the later case could happen
            ### if the line is out of service
            if orientation == 0 or line not in m.pf:
                continue
            elif orientation == 1:
                expr += m.pf[line]
            elif orientation == -1:
                expr -= m.pf[line]
            else:
                raise Exception("line_orientation must be in [-1,0,1], found "\
                        "line_orientation {} for line {} in interface {}".format(
                            orientation, line, interface_name))
        m.eq_pf_interface[interface_name] = m.pfi[interface_name] == expr


def get_power_flow_expr_ptdf_approx(model, branch_name, PTDF, rel_ptdf_tol=None, abs_ptdf_tol=None):
    """
    Create a pyomo power flow expression from PTDF matrix
    """

    if rel_ptdf_tol is None:
        rel_ptdf_tol = 0.
    if abs_ptdf_tol is None:
        abs_ptdf_tol = 0.

    const = PTDF.get_branch_phase_shift(branch_name) + PTDF.get_branch_phi_adj(branch_name)

    max_coef = PTDF.get_branch_ptdf_abs_max(branch_name)

    ptdf_tol = max(abs_ptdf_tol, rel_ptdf_tol*max_coef) 
    ## NOTE: It would be easy to hold on to the 'ptdf' dictionary here,
    ##       if we wanted to
    m_p_nw = model.p_nw
    ## if model.p_nw is Var, we can use LinearExpression
    ## to build these dense constraints much faster
    if isinstance(m_p_nw, pe.Var):
        coef_list = list()
        var_list = list()
        for bus_name, coef in PTDF.get_branch_ptdf_iterator(branch_name):
            if abs(coef) >= ptdf_tol:
                coef_list.append(coef)
                var_list.append(m_p_nw[bus_name])

        lin_expr_list = [const] + coef_list + var_list 
        expr = LinearExpression(lin_expr_list)
    else:
        expr = quicksum( (coef*m_p_nw[bus_name] for bus_name, coef in PTDF.get_branch_ptdf_iterator(branch_name) if abs(coef) >= ptdf_tol), start=const, linear=True)

    return expr


def declare_eq_branch_power_ptdf_approx(model, index_set, PTDF, rel_ptdf_tol=None, abs_ptdf_tol=None):
    """
    Create the equality constraints or expressions for power (from PTDF 
    approximation) in the branch
    """

    m = model

    con_set = decl.declare_set("_con_eq_branch_power_ptdf_approx_set", model, index_set)

    pf_is_var = isinstance(m.pf, pe.Var)

    if pf_is_var:
        m.eq_pf_branch = pe.Constraint(con_set)
    else:
        if not isinstance(m.pf, pe.Expression):
            raise Exception("Unrecognized type for m.pf", m.pf.pprint())

    for branch_name in con_set:
        expr = \
            get_power_flow_expr_ptdf_approx(m, branch_name, PTDF, rel_ptdf_tol=rel_ptdf_tol, abs_ptdf_tol=abs_ptdf_tol)

        if pf_is_var:
            m.eq_pf_branch[branch_name] = \
                m.pf[branch_name] == expr
        else:
            m.pf[branch_name] = expr


def get_branch_loss_expr_ptdf_approx(model, branch_name, PTDF, rel_ptdf_tol=None, abs_ptdf_tol=None): 
    """
    Create a pyomo power flow loss expression from PTDF matrix
    """
    if rel_ptdf_tol is None:
        rel_ptdf_tol = 0.
    if abs_ptdf_tol is None:
        abs_ptdf_tol = 0.

    const = PTDF.get_branch_losses_phase_shift(branch_name)
    const += PTDF.get_branch_ldf_c(branch_name)
    const += PTDF.get_branch_phi_losses_adj(branch_name)

    max_coef = PTDF.get_branch_ldf_abs_max(branch_name)

    ptdf_tol = max(abs_ptdf_tol, rel_ptdf_tol*max_coef) 
    ## NOTE: It would be easy to hold on to the 'ptdf' dictionary here,
    ##       if we wanted to
    m_p_nw = model.p_nw
    ## if model.p_nw is Var, we can use LinearExpression
    ## to build these dense constraints much faster
    if isinstance(m_p_nw, pe.Var):
        coef_list = list()
        var_list = list()
        for bus_name, coef in PTDF.get_branch_ldf_iterator(branch_name):
            if abs(coef) >= ptdf_tol:
                coef_list.append(coef)
                var_list.append(m_p_nw[bus_name])

        lin_expr_list = [const] + coef_list + var_list 
        expr = LinearExpression(lin_expr_list)
    else:
        expr = quicksum( (coef*m_p_nw[bus_name] for bus_name, coef in PTDF.get_branch_ldf_iterator(branch_name) if abs(coef) >= ptdf_tol), start=const, linear=True)

    return expr


def declare_eq_branch_loss_ptdf_approx(model, index_set, PTDF, rel_ptdf_tol=None, abs_ptdf_tol=None):
    """
    Create the equality constraints or expressions for losses (from PTDF 
    approximation) in the branch
    """
    m = model

    con_set = decl.declare_set("_con_eq_branch_loss_ptdf_approx_set", model, index_set)
    pfl_is_var = isinstance(m.pfl, pe.Var)
    if pfl_is_var:
        m.eq_pfl_branch = pe.Constraint(con_set)
    else:
        if not isinstance(m.pfl, pe.Expression):
            raise Exception("Unrecognized type for m.pfl", m.pfl.pprint())

    for branch_name in con_set:
        expr = \
            get_branch_loss_expr_ptdf_approx(m, branch_name, PTDF, rel_ptdf_tol=rel_ptdf_tol, abs_ptdf_tol=abs_ptdf_tol)

        if pfl_is_var:
            m.eq_pfl_branch[branch_name] = \
                m.pfl[branch_name] == expr
        else:
            m.pfl[branch_name] = expr


def get_power_flow_interface_expr_ptdf(model, interface_name, PTDF, rel_ptdf_tol=None, abs_ptdf_tol=None):
    """
    Create a pyomo power flow expression from PTDF matrix for an interface
    """
    if rel_ptdf_tol is None:
        rel_ptdf_tol = 0.
    if abs_ptdf_tol is None:
        abs_ptdf_tol = 0.

    const = PTDF.get_interface_const(interface_name)
    max_coef = PTDF.get_interface_ptdf_abs_max(interface_name)

    ptdf_tol = max(abs_ptdf_tol, rel_ptdf_tol*max_coef)

    m_p_nw = model.p_nw

    ## if model.p_nw is Var, we can use LinearExpression
    ## to build these dense constraints much faster
    if isinstance(m_p_nw, pe.Var):
        coef_list = list()
        var_list = list()
        for bus_name, coef in PTDF.get_interface_ptdf_iterator(interface_name):
            if abs(coef) >= ptdf_tol:
                coef_list.append(coef)
                var_list.append(m_p_nw[bus_name])

        lin_expr_list = [const] + coef_list + var_list
        expr = LinearExpression(lin_expr_list)
    else:
        expr = quicksum( (coef*m_p_nw[bus_name] for bus_name, coef in PTDF.get_interface_ptdf_iterator(interface_name) if abs(coef) >= ptdf_tol), start=const, linear=True)

    return expr


def declare_eq_interface_power_ptdf_approx(model, index_set, PTDF, rel_ptdf_tol=None, abs_ptdf_tol=None):
    """
    Create equality constraints or expressions for power (from PTDF
    approximation) across the interface
    """

    m = model
    con_set = decl.declare_set("_con_eq_interface_power_ptdf_approx_set", model, index_set)

    pfi_is_var = isinstance(m.pfi, pe.Var)

    if pfi_is_var:
        m.eq_pf_interface = pe.Constraint(con_set)
    else:
        if not isinstance(m.pfi, pe.Expression):
            raise Exception("Unrecognized type for m.pfi", m.pfi.pprint())

    for interface_name in con_set:
        expr = \
            get_power_flow_interface_expr_ptdf(m, interface_name, PTDF,
                    rel_ptdf_tol=rel_ptdf_tol, abs_ptdf_tol=abs_ptdf_tol)

        if pfi_is_var:
            m.eq_pf_interface[interface_name] = \
                    m.pfi[interface_name] == expr
        else:
            m.pfi[interface_name] = expr


def declare_ineq_s_branch_thermal_limit(model, index_set,
                                        branches, s_thermal_limits,
                                        flow_type=FlowType.POWER):
    """
    Create the inequality constraints for the branch thermal limits
    based on the power variables.
    """
    m = model
    con_set = decl.declare_set('_con_ineq_s_branch_thermal_limit',
                               model=model, index_set=index_set)

    m.ineq_sf_branch_thermal_limit = pe.Constraint(con_set)
    m.ineq_st_branch_thermal_limit = pe.Constraint(con_set)

    if flow_type == FlowType.CURRENT:
        for branch_name in con_set:
            if s_thermal_limits[branch_name] is None:
                continue

            from_bus = branches[branch_name]['from_bus']
            to_bus = branches[branch_name]['to_bus']
            m.ineq_sf_branch_thermal_limit[branch_name] = \
                (m.vr[from_bus] ** 2 + m.vj[from_bus] ** 2) * (m.ifr[branch_name] ** 2 + m.ifj[branch_name] ** 2) \
                <= s_thermal_limits[branch_name] ** 2
            m.ineq_st_branch_thermal_limit[branch_name] = \
                (m.vr[to_bus] ** 2 + m.vj[to_bus] ** 2) * (m.itr[branch_name] ** 2 + m.itj[branch_name] ** 2) \
                <= s_thermal_limits[branch_name] ** 2
    elif flow_type == FlowType.POWER:
        for branch_name in con_set:
            if s_thermal_limits[branch_name] is None:
                continue

            m.ineq_sf_branch_thermal_limit[branch_name] = \
                m.pf[branch_name] ** 2 + m.qf[branch_name] ** 2 \
                <= s_thermal_limits[branch_name] ** 2
            m.ineq_st_branch_thermal_limit[branch_name] = \
                m.pt[branch_name] ** 2 + m.qt[branch_name] ** 2 \
                <= s_thermal_limits[branch_name] ** 2


def declare_ineq_p_branch_thermal_lbub(model, index_set,
                                        branches, p_thermal_limits,
                                        approximation_type=ApproximationType.BTHETA,
                                        slacks=False):
    """
    Create the inequality constraints for the branch thermal limits
    based on the power variables or expressions.
    """
    m = model
    con_set = decl.declare_set('_con_ineq_p_branch_thermal_lbub',
                               model=model, index_set=index_set)

    # flag for if slacks are on the model
    if slacks:
        if not hasattr(model, 'pf_slack_pos'):
            raise Exception('No positive slack branch variables on model, but slacks=True')
        if not hasattr(model, 'pf_slack_neg'):
            raise Exception('No negative slack branch variables on model, but slacks=True')

    m.ineq_pf_branch_thermal_lb = pe.Constraint(con_set)
    m.ineq_pf_branch_thermal_ub = pe.Constraint(con_set)

    if approximation_type == ApproximationType.BTHETA or \
            approximation_type == ApproximationType.PTDF:
        for branch_name in con_set:
            if p_thermal_limits[branch_name] is None:
                continue

            if slacks and branch_name in m.pf_slack_neg:
                pf_bn = m.pf[branch_name]
                if hasattr(pf_bn, 'expr') and isinstance(pf_bn.expr, LinearExpression):
                    ## create a copy
                    old_expr = pf_bn.expr
                    expr = LinearExpression(constant=old_expr.constant,
                                            linear_vars = old_expr.linear_vars[:] + [m.pf_slack_neg[branch_name]],
                                            linear_coefs = old_expr.linear_coefs[:] + [1],
                                            )

                else:
                    expr = m.pf[branch_name] + m.pf_slack_neg[branch_name]
                m.ineq_pf_branch_thermal_lb[branch_name] = \
                    (-p_thermal_limits[branch_name], expr, None)
            else:
                m.ineq_pf_branch_thermal_lb[branch_name] = \
                    (-p_thermal_limits[branch_name], m.pf[branch_name], None)

            if slacks and branch_name in m.pf_slack_pos:
                pf_bn = m.pf[branch_name]
                if hasattr(pf_bn, 'expr') and isinstance(pf_bn.expr, LinearExpression):
                    ## create a copy
                    old_expr = pf_bn.expr
                    expr = LinearExpression(constant=old_expr.constant,
                                            linear_vars = old_expr.linear_vars[:] + [m.pf_slack_pos[branch_name]],
                                            linear_coefs = old_expr.linear_coefs[:] + [-1],
                                            )
                else:
                    expr = m.pf[branch_name] - m.pf_slack_pos[branch_name]
                m.ineq_pf_branch_thermal_lb[branch_name] = \
                    (None, expr, p_thermal_limits[branch_name])
            else:
                m.ineq_pf_branch_thermal_ub[branch_name] = \
                    (None, m.pf[branch_name], p_thermal_limits[branch_name])

def generate_thermal_bounds(pf, llimit, ulimit, neg_slack=None, pos_slack=None):
    """
    Create a constraint for thermal limits on a line given the power flow
    expression or variable pf, a lower limit llimit, a uppder limit ulimit,
    and the negative slack variable, neg_slack, (None if not needed) and
    positive slack variable, pos_slack, (None if not needed) added to this 
    constraint.
    """
    if hasattr(pf, 'expr') and isinstance(pf.expr, LinearExpression):
        ## if necessary, copy again, so that m.pf[bn] **is** the flow
        add_vars = list()
        add_coefs = list()
        if neg_slack is not None:
            add_vars.append(neg_slack)
            add_coefs.append(1)
        if pos_slack is not None:
            add_vars.append(pos_slack)
            add_coefs.append(-1)
        if add_vars:
            ## create a copy
            old_expr = pf.expr
            expr = LinearExpression(constant = old_expr.constant,
                                    linear_vars = old_expr.linear_vars[:] + add_vars,
                                    linear_coefs = old_expr.linear_coefs[:] + add_coefs,
                                   )
        else:
            expr = pf
    else:
        expr = pf
        if neg_slack is not None:
            expr += neg_slack
        if pos_slack is not None:
            expr -= pos_slack
    return (llimit, expr, ulimit)

def declare_ineq_p_branch_thermal_bounds(model, index_set,
                                        branches, p_thermal_limits,
                                        approximation_type=ApproximationType.BTHETA,
                                        slacks=False):
    """
    Create an inequality constraint for the branch thermal limits
    based on the power variables or expressions.
    """
    m = model
    con_set = decl.declare_set('_con_ineq_p_branch_thermal_lbub',
                               model=model, index_set=index_set)
    # flag for if slacks are on the model
    if slacks:
        if not hasattr(model, 'pf_slack_pos'):
            raise Exception('No positive slack branch variables on model, but slacks=True')
        if not hasattr(model, 'pf_slack_neg'):
            raise Exception('No negative slack branch variables on model, but slacks=True')

    m.ineq_pf_branch_thermal_bounds = pe.Constraint(con_set)

    if approximation_type == ApproximationType.BTHETA or \
            approximation_type == ApproximationType.PTDF:
        for branch_name in con_set:
            limit = p_thermal_limits[branch_name]
            if limit is None:
                continue

            if slacks and branch_name in m.pf_slack_neg:
                neg_slack = m.pf_slack_neg[branch_name]
            else:
                neg_slack = None

            if slacks and branch_name in m.pf_slack_pos:
                pos_slack = m.pf_slack_pos[branch_name]
            else:
                pos_slack = None

            m.ineq_pf_branch_thermal_bounds[branch_name] = \
                    generate_thermal_bounds(m.pf[branch_name], -limit, limit, neg_slack, pos_slack)

def declare_ineq_angle_diff_branch_lbub_c_s(model, index_set, branches):
    """
    Create the inequality constraints for the angle difference
    bounds between interconnected buses.
    """
    m = model
    con_set = decl.declare_set('_con_ineq_angle_diff_branch_lbub',
                               model=model, index_set=index_set)

    m.ineq_angle_diff_branch_lb = pe.Constraint(con_set)
    m.ineq_angle_diff_branch_ub = pe.Constraint(con_set)

    for branch_name in con_set:
        from_bus = branches[branch_name]['from_bus']
        to_bus = branches[branch_name]['to_bus']

        m.ineq_angle_diff_branch_lb[branch_name] = (math.tan(math.radians(branches[branch_name]['angle_diff_min'])) *
                                                    m.c[(from_bus, to_bus)] <= m.s[(from_bus, to_bus)])
        m.ineq_angle_diff_branch_ub[branch_name] = (m.s[(from_bus, to_bus)] <=
                                                    math.tan(math.radians(branches[branch_name]['angle_diff_max'])) *
                                                    m.c[(from_bus, to_bus)])


def declare_ineq_angle_diff_branch_lbub(model, index_set, branches, coordinate_type=CoordinateType.POLAR):
    """
    Create the inequality constraints for the angle difference
    bounds between interconnected buses.
    """
    m = model
    con_set = decl.declare_set('_con_ineq_angle_diff_branch_lbub',
                               model=model, index_set=index_set)

    m.ineq_angle_diff_branch_lb = pe.Constraint(con_set)
    m.ineq_angle_diff_branch_ub = pe.Constraint(con_set)

    if coordinate_type == CoordinateType.POLAR:
        for branch_name in con_set:
            from_bus = branches[branch_name]['from_bus']
            to_bus = branches[branch_name]['to_bus']

            m.ineq_angle_diff_branch_lb[branch_name] = \
                math.radians(branches[branch_name]['angle_diff_min']) <= m.va[from_bus] - m.va[to_bus]
            m.ineq_angle_diff_branch_ub[branch_name] = \
                m.va[from_bus] - m.va[to_bus] <= math.radians(branches[branch_name]['angle_diff_max'])
    elif coordinate_type == CoordinateType.RECTANGULAR:
        for branch_name in con_set:
            from_bus = branches[branch_name]['from_bus']
            to_bus = branches[branch_name]['to_bus']

            m.ineq_angle_diff_branch_lb[branch_name] = (math.tan(math.radians(branches[branch_name]['angle_diff_min'])) *
                                                        (m.vr[from_bus] * m.vr[to_bus] + m.vj[from_bus] * m.vj[to_bus]) <=
                                                        m.vj[from_bus] * m.vr[to_bus] - m.vr[from_bus] * m.vj[to_bus])
            m.ineq_angle_diff_branch_ub[branch_name] = (m.vj[from_bus] * m.vr[to_bus] - m.vr[from_bus] * m.vj[to_bus] <=
                                                        math.tan(math.radians(branches[branch_name]['angle_diff_max'])) *
                                                        (m.vr[from_bus] * m.vr[to_bus] + m.vj[from_bus] * m.vj[to_bus]))


def declare_ineq_p_interface_bounds(model, index_set, interfaces,
                                        approximation_type=ApproximationType.BTHETA,
                                        slacks=False):
    """
    Create the inequality constraints for the interface limits
    based on the power variables or expressions.

    p_interface_limits should be (lower, upper) tuple
    """
    m = model
    con_set = decl.declare_set('_con_ineq_p_interface_bounds',
                               model=model, index_set=index_set)

    m.ineq_pf_interface_bounds = pe.Constraint(con_set)

    # flag for if slacks are on the model
    if slacks:
        if not hasattr(model, 'pfi_slack_pos'):
            raise Exception('No positive slack interface variables on model, but slacks=True')
        if not hasattr(model, 'pfi_slack_neg'):
            raise Exception('No negative slack interface variables on model, but slacks=True')

    if approximation_type == ApproximationType.BTHETA or \
            approximation_type == ApproximationType.PTDF:
        for interface_name in con_set:
            interface = interfaces[interface_name]
            if interface['minimum_limit'] is not None or \
                    interface['maximum_limit'] is not None:

                if slacks and interface_name in m.pfi_slack_neg:
                    neg_slack = m.pfi_slack_neg[interface_name]
                else:
                    neg_slack = None

                if slacks and interface_name in m.pfi_slack_pos:
                    pos_slack = m.pfi_slack_pos[interface_name]
                else:
                    pos_slack = None

                m.ineq_pf_interface_bounds[interface_name] = \
                    generate_thermal_bounds(m.pfi[interface_name], interface['minimum_limit'], interface['maximum_limit'],
                                            neg_slack, pos_slack)
