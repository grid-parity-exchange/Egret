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

    found_q_costs = False
    for gen_name in expr_set:
        if p_costs is not None and gen_name in p_costs:
            #TODO: We assume that the costs are polynomial type
            assert p_costs[gen_name]['cost_curve_type'] == 'polynomial'
            m.pg_operating_cost[gen_name] = \
                sum(v*m.pg[gen_name]**i for i, v in p_costs[gen_name]['values'].items())

        if q_costs is not None and gen_name in q_costs:
            #TODO: We assume that the costs are polynomial type
            assert q_costs['cost_curve_type'] == 'polynomial'
            found_q_costs = True
            m.qg_operating_cost[gen_name] = \
                sum(v*m.qg[gen_name]**i for i, v in q_costs[gen_name]['values'].items())

    if not found_q_costs:
        del m.qg_operating_cost
