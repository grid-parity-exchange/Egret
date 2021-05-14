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


def declare_var_pg_cost(model, index_set, **kwargs):
    decl.declare_var('pg_cost', model, index_set=index_set, **kwargs)


def declare_var_qg_cost(model, index_set, **kwargs):
    decl.declare_var('qg_cost', model, index_set=index_set, **kwargs)


def _pw_cost_helper(cost_dict, cost_var, gen_var, pw_cost_set, gen_name, indexed_pw_cost_con):
    if cost_dict['cost_curve_type'] == 'polynomial':
        pass
    elif cost_dict['cost_curve_type'] == 'piecewise':
        pt0, cost0 = cost_dict['values'][0]
        if gen_var.lb < pt0:
            raise ValueError('Piecewise costs require that the lower bound on pg be greater than or equal to the first piecewise point.')

        last_slope = None
        for ndx in range(len(cost_dict['values']) - 1):
            pt1, cost1 = cost_dict['values'][ndx]
            pt2, cost2 = cost_dict['values'][ndx + 1]
            slope = (cost2 - cost1) / (pt2 - pt1)
            if last_slope is not None:
                if slope < last_slope:
                    all_slopes = [(cost_dict['values'][i + 1][1] - cost_dict['values'][i][1]) / (cost_dict['values'][i + 1][0] - cost_dict['values'][i][0]) for i in range(len(cost_dict['values']) - 1)]
                    raise ValueError(f'Piecewise costs must be convex; generator: {gen_name}; slopes: {all_slopes}')
            intercept = cost2 - slope * pt2
            pw_cost_set.add((gen_name, ndx))
            indexed_pw_cost_con[gen_name, ndx] = cost_var >= slope * gen_var + intercept
            last_slope = slope
    else:
        raise ValueError(f"Unrecognized cost_cureve_type: {cost_dict['cost_curve_type']}")


def declare_piecewise_cost_cons(model, index_set, p_costs, q_costs=None):
    m = model

    m.pg_piecewise_cost_set = pe.Set(dimen=2)
    m.qg_piecewise_cost_set = pe.Set(dimen=2)
    m.pg_piecewise_cost_cons = pe.Constraint(m.pg_piecewise_cost_set)
    m.qg_piecewise_cost_cons = pe.Constraint(m.qg_piecewise_cost_set)

    for gen_name in index_set:
        if p_costs is not None and gen_name in p_costs:
            _pw_cost_helper(cost_dict=p_costs[gen_name],
                            cost_var=m.pg_cost[gen_name],
                            gen_var=m.pg[gen_name],
                            pw_cost_set=m.pg_piecewise_cost_set,
                            gen_name=gen_name,
                            indexed_pw_cost_con=m.pg_piecewise_cost_cons)

        if q_costs is not None and gen_name in q_costs:
            _pw_cost_helper(cost_dict=q_costs[gen_name],
                            cost_var=m.qg_cost[gen_name],
                            gen_var=m.qg[gen_name],
                            pw_cost_set=m.qg_piecewise_cost_set,
                            gen_name=gen_name,
                            indexed_pw_cost_con=m.qg_piecewise_cost_cons)


def declare_expression_pgqg_operating_cost(model, index_set,
                                           p_costs, q_costs=None):
    """
    Create the Expression objects to represent the operating costs
    for the real and reactive (if present) power of each of the
    generators.
    """
    m = model
    expr_set = decl.declare_set('_expr_g_operating_cost',
                                model=model, index_set=index_set)
    m.pg_operating_cost = pe.Expression(expr_set)
    m.qg_operating_cost = pe.Expression(expr_set)

    for gen_name in expr_set:
        if p_costs is not None and gen_name in p_costs:
            if p_costs[gen_name]['cost_curve_type'] == 'polynomial':
                m.pg_operating_cost[gen_name] = sum(v*m.pg[gen_name]**i for i, v in p_costs[gen_name]['values'].items())
            elif p_costs[gen_name]['cost_curve_type'] == 'piecewise':
                m.pg_operating_cost[gen_name] = m.pg_cost[gen_name]
            else:
                raise ValueError(f"Unrecognized cost_cureve_type: {p_costs[gen_name]['cost_curve_type']}")
        else:
            m.pg_operating_cost[gen_name] = 0

        if q_costs is not None and gen_name in q_costs:
            if q_costs[gen_name]['cost_curve_type'] == 'polynomial':
                m.qg_operating_cost[gen_name] = sum(v*m.qg[gen_name]**i for i, v in q_costs[gen_name]['values'].items())
            elif q_costs[gen_name]['cost_curve_type'] == 'piecewise':
                m.qg_operating_cost[gen_name] = m.qg_cost[gen_name]
            else:
                raise ValueError(f"Unrecognized cost_cureve_type: {q_costs[gen_name]['cost_curve_type']}")
        else:
            m.qg_operating_cost[gen_name] = 0
