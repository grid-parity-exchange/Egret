#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module contains the modeling components used when modeling generators
in transmission models
"""
import pyomo.environ as pe
import egret.model_library.decl as decl
from egret.model_library.transmission import tx_utils
from pyomo.core.expr.numeric_expr import LinearExpression


def declare_var_pg(model, index_set, **kwargs):
    """
    Create a variable for the real component of the power at a generator
    """
    decl.declare_var('pg', model=model, index_set=index_set, **kwargs)


def declare_var_qg(model, index_set, **kwargs):
    """
    Create a variable for the reactive component of the power at a generator
    """
    decl.declare_var('qg', model=model, index_set=index_set, **kwargs)


def pw_gen_generator(index_set, costs):
    for gen_name in index_set:
        if gen_name not in costs:
            continue
        curve = costs[gen_name]
        assert curve['cost_curve_type'] in {'piecewise', 'polynomial'}
        if curve['cost_curve_type'] != 'piecewise':
            continue
        yield gen_name


def declare_var_delta_pg(model, index_set, p_costs):
    m = model

    m.delta_pg_set = pe.Set(dimen=2)
    m.delta_pg = pe.Var(m.delta_pg_set)
    for gen_name in pw_gen_generator(index_set, p_costs):
        p_min = m.pg[gen_name].lb
        p_max = m.pg[gen_name].ub
        curve = p_costs[gen_name]
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=p_min,
                                                                p_max=p_max,
                                                                gen_name=gen_name)
        for ndx, ((o1, c1), (o2, c2)) in enumerate(zip(cleaned_values, cleaned_values[1:])):
            m.delta_pg_set.add((gen_name, ndx))
            m.delta_pg[gen_name, ndx].setlb(0)
            m.delta_pg[gen_name, ndx].setub(o2 - o1)


def declare_var_delta_qg(model, index_set, q_costs):
    m = model

    m.delta_qg_set = pe.Set(dimen=2)
    m.delta_qg = pe.Var(m.delta_qg_set)
    for gen_name in pw_gen_generator(index_set, q_costs):
        q_min = m.qg[gen_name].lb
        q_max = m.qg[gen_name].ub
        curve = q_costs[gen_name]
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=q_min,
                                                                p_max=q_max,
                                                                gen_name=gen_name)
        for ndx, ((o1, c1), (o2, c2)) in enumerate(zip(cleaned_values, cleaned_values[1:])):
            m.delta_qg_set.add((gen_name, ndx))
            m.delta_qg[gen_name, ndx].setlb(0)
            m.delta_qg[gen_name, ndx].setub(o2 - o1)


def declare_pg_delta_pg_con(model, index_set, p_costs):
    m = model

    m.pg_delta_pg_con_set = pe.Set()
    m.pg_delta_pg_con = pe.Constraint(m.pg_delta_pg_con_set)
    for gen_name in pw_gen_generator(index_set, p_costs):
        p_min = m.pg[gen_name].lb
        p_max = m.pg[gen_name].ub
        curve = p_costs[gen_name]
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=p_min,
                                                                p_max=p_max,
                                                                gen_name=gen_name)
        m.pg_delta_pg_con_set.add(gen_name)
        lin_coefs = [-1]
        lin_vars = [m.pg[gen_name]]
        for ndx, ((o1, c1), (o2, c2)) in enumerate(zip(cleaned_values, cleaned_values[1:])):
            lin_coefs.append(1)
            lin_vars.append(m.delta_pg[gen_name, ndx])
        expr = LinearExpression(constant=cleaned_values[0][0], linear_coefs=lin_coefs, linear_vars=lin_vars)
        m.pg_delta_pg_con[gen_name] = (expr, 0)


def declare_qg_delta_qg_con(model, index_set, q_costs):
    m = model

    m.qg_delta_qg_con_set = pe.Set()
    m.qg_delta_qg_con = pe.Constraint(m.qg_delta_qg_con_set)
    for gen_name in pw_gen_generator(index_set, q_costs):
        q_min = m.qg[gen_name].lb
        q_max = m.qg[gen_name].ub
        curve = q_costs[gen_name]
        cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                curve_type='cost_curve',
                                                                p_min=q_min,
                                                                p_max=q_max,
                                                                gen_name=gen_name)
        m.qg_delta_qg_con_set.add(gen_name)
        lin_coefs = [-1]
        lin_vars = [m.qg[gen_name]]
        for ndx, ((o1, c1), (o2, c2)) in enumerate(zip(cleaned_values, cleaned_values[1:])):
            lin_coefs.append(1)
            lin_vars.append(m.delta_qg[gen_name, ndx])
        expr = LinearExpression(constant=cleaned_values[0][0], linear_coefs=lin_coefs, linear_vars=lin_vars)
        m.qg_delta_qg_con[gen_name] = (expr, 0)


def declare_var_pg_cost(model, index_set, p_costs):
    m = model
    actual_indices = list(pw_gen_generator(index_set, p_costs))
    m.pg_cost_set = pe.Set(initialize=actual_indices)
    m.pg_cost = pe.Var(m.pg_cost_set)


def declare_var_qg_cost(model, index_set, q_costs):
    m = model
    actual_indices = list(pw_gen_generator(index_set, q_costs))
    m.qg_cost_set = pe.Set(initialize=actual_indices)
    m.qg_cost = pe.Var(m.qg_cost_set)


def _pw_cost_helper(cost_dict, cost_var, gen_var, pw_cost_set, gen_name, indexed_pw_cost_con):
    if cost_dict['cost_curve_type'] == 'polynomial':
        pass
    elif cost_dict['cost_curve_type'] == 'piecewise':
        cleaned_values = tx_utils.validate_and_clean_cost_curve(cost_dict,
                                                                curve_type='cost_curve',
                                                                p_min=gen_var.lb,
                                                                p_max=gen_var.ub,
                                                                gen_name=gen_name)
        if len(cleaned_values) > 1:
            for ndx, ((pt1, cost1), (pt2, cost2)) in enumerate(zip(cleaned_values, cleaned_values[1:])):
                slope = (cost2 - cost1) / (pt2 - pt1)
                intercept = cost2 - slope * pt2
                pw_cost_set.add((gen_name, ndx))
                indexed_pw_cost_con[gen_name, ndx] = cost_var >= slope * gen_var + intercept
        else:
            intercept = cleaned_values[0][1]
            pw_cost_set.add((gen_name, 0))
            indexed_pw_cost_con[gen_name, 0] = cost_var == intercept
    else:
        raise ValueError(f"Unrecognized cost_curve_type: {cost_dict['cost_curve_type']}")


def declare_piecewise_pg_cost_cons(model, index_set, p_costs):
    m = model

    m.pg_piecewise_cost_set = pe.Set(dimen=2)
    m.pg_piecewise_cost_cons = pe.Constraint(m.pg_piecewise_cost_set)

    for gen_name in pw_gen_generator(index_set=index_set, costs=p_costs):
        _pw_cost_helper(cost_dict=p_costs[gen_name],
                        cost_var=m.pg_cost[gen_name],
                        gen_var=m.pg[gen_name],
                        pw_cost_set=m.pg_piecewise_cost_set,
                        gen_name=gen_name,
                        indexed_pw_cost_con=m.pg_piecewise_cost_cons)


def declare_piecewise_qg_cost_cons(model, index_set, q_costs):
    m = model

    m.qg_piecewise_cost_set = pe.Set(dimen=2)
    m.qg_piecewise_cost_cons = pe.Constraint(m.qg_piecewise_cost_set)

    for gen_name in pw_gen_generator(index_set=index_set, costs=q_costs):
        _pw_cost_helper(cost_dict=q_costs[gen_name],
                        cost_var=m.qg_cost[gen_name],
                        gen_var=m.qg[gen_name],
                        pw_cost_set=m.qg_piecewise_cost_set,
                        gen_name=gen_name,
                        indexed_pw_cost_con=m.qg_piecewise_cost_cons)


def declare_expression_pg_operating_cost(model, index_set, p_costs, pw_formulation='delta'):
    """
    Create the Expression objects to represent the operating costs
    for the real power of each of the generators.
    """
    m = model
    expr_set = decl.declare_set('_expr_pg_operating_cost',
                                model=model, index_set=index_set)
    m.pg_operating_cost = pe.Expression(expr_set)

    for gen_name in expr_set:
        if gen_name in p_costs:
            if p_costs[gen_name]['cost_curve_type'] == 'polynomial':
                m.pg_operating_cost[gen_name] = sum(v*m.pg[gen_name]**i for i, v in p_costs[gen_name]['values'].items())
            elif p_costs[gen_name]['cost_curve_type'] == 'piecewise':
                if pw_formulation == 'delta':
                    p_min = m.pg[gen_name].lb
                    p_max = m.pg[gen_name].ub
                    curve = p_costs[gen_name]
                    cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                            curve_type='cost_curve',
                                                                            p_min=p_min,
                                                                            p_max=p_max,
                                                                            gen_name=gen_name)
                    expr = cleaned_values[0][1]
                    if len(cleaned_values) > 1:
                        for ndx, ((o1, c1), (o2, c2)) in enumerate(zip(cleaned_values, cleaned_values[1:])):
                            slope = (c2 - c1) / (o2 - o1)
                            expr += slope * m.delta_pg[gen_name, ndx]
                    m.pg_operating_cost[gen_name] = expr
                else:
                    m.pg_operating_cost[gen_name] = m.pg_cost[gen_name]
            else:
                raise ValueError(f"Unrecognized cost_cureve_type: {p_costs[gen_name]['cost_curve_type']}")
        else:
            m.pg_operating_cost[gen_name] = 0


def declare_expression_qg_operating_cost(model, index_set, q_costs, pw_formulation='delta'):
    """
    Create the Expression objects to represent the operating costs
    for the reactive power of each of the generators.
    """
    m = model
    expr_set = decl.declare_set('_expr_qg_operating_cost',
                                model=model, index_set=index_set)
    m.qg_operating_cost = pe.Expression(expr_set)

    for gen_name in expr_set:
        if gen_name in q_costs:
            if q_costs[gen_name]['cost_curve_type'] == 'polynomial':
                m.qg_operating_cost[gen_name] = sum(v*m.qg[gen_name]**i for i, v in q_costs[gen_name]['values'].items())
            elif q_costs[gen_name]['cost_curve_type'] == 'piecewise':
                if pw_formulation == 'delta':
                    q_min = m.qg[gen_name].lb
                    q_max = m.qg[gen_name].ub
                    curve = q_costs[gen_name]
                    cleaned_values = tx_utils.validate_and_clean_cost_curve(curve=curve,
                                                                            curve_type='cost_curve',
                                                                            p_min=q_min,
                                                                            p_max=q_max,
                                                                            gen_name=gen_name)
                    expr = cleaned_values[0][1]
                    for ndx, ((o1, c1), (o2, c2)) in enumerate(zip(cleaned_values, cleaned_values[1:])):
                        slope = (c2 - c1) / (o2 - o1)
                        expr += slope * m.delta_pg[gen_name, ndx]
                    m.qg_operating_cost[gen_name] = expr
                else:
                    m.qg_operating_cost[gen_name] = m.qg_cost[gen_name]
            else:
                raise ValueError(f"Unrecognized cost_cureve_type: {q_costs[gen_name]['cost_curve_type']}")
        else:
            m.qg_operating_cost[gen_name] = 0
