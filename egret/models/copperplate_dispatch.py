#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module provides functions that create the modules for typical copperplate dispatch formulations.

#TODO: document this with examples
"""
import pyomo.environ as pe
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen

from egret.model_library.defn import ApproximationType
from egret.data.data_utils import map_items, zip_items
from math import pi

def _include_system_feasibility_slack(model, bus_p_loads, gen_attrs, p_marginal_slack_penalty):
    import egret.model_library.decl as decl
    load = sum(bus_p_loads.values())

    over_gen_bounds = (0, tx_utils.over_gen_limit(load, gen_attrs['names'], gen_attrs['p_max']))
    decl.declare_var('p_over_generation', model=model, index_set=None,
                     initialize=0., bounds=over_gen_bounds
                     )

    load_shed_bounds = (0, tx_utils.load_shed_limit(load, gen_attrs['names'], gen_attrs['p_min']))
    decl.declare_var('p_load_shed', model=model, index_set=None,
                     initialize=0., bounds=load_shed_bounds
                     )
    p_rhs_kwargs = {'include_feasibility_load_shed':'p_load_shed', 'include_feasibility_over_generation':'p_over_generation'}
    penalty_expr = p_marginal_slack_penalty * (model.p_load_shed + model.p_over_generation)
    return p_rhs_kwargs, penalty_expr

def _validate_and_extract_slack_penalty(model_data):
    assert('load_mismatch_cost' in model_data.data['system'])
    return model_data.data['system']['load_mismatch_cost']

def create_copperplate_dispatch_approx_model(model_data, include_feasibility_slack=False):
    md = model_data.clone_in_service()
    tx_utils.scale_ModelData_to_pu(md, inplace = True)

    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    loads = dict(md.elements(element_type='load'))
    shunts = dict(md.elements(element_type='shunt'))

    gen_attrs = md.attributes(element_type='generator')
    bus_attrs = md.attributes(element_type='bus')

    inlet_branches_by_bus, outlet_branches_by_bus = \
        tx_utils.inlet_outlet_branches_by_bus(branches, buses)
    gens_by_bus = tx_utils.gens_by_bus(buses, gens)

    model = pe.ConcreteModel()

    ### declare (and fix) the loads at the buses
    bus_p_loads, _ = tx_utils.dict_of_bus_loads(buses, loads)

    libbus.declare_var_pl(model, bus_attrs['names'], initialize=bus_p_loads)
    model.pl.fix()

    ### declare the fixed shunts at the buses
    _, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    ### declare the generator real power
    pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model, gen_attrs['names'], initialize=pg_init,
                          bounds=zip_items(gen_attrs['p_min'], gen_attrs['p_max'])
                          )

    ### include the feasibility slack for the system balance
    p_rhs_kwargs = {}
    if include_feasibility_slack:
        p_marginal_slack_penalty = _validate_and_extract_slack_penalty(model_data)                
        p_rhs_kwargs, penalty_expr = _include_system_feasibility_slack(model, bus_p_loads, gen_attrs, p_marginal_slack_penalty)

    ### declare the p balance
    libbus.declare_eq_p_balance_ed(model=model,
                                   index_set=bus_attrs['names'],
                                   bus_p_loads=bus_p_loads,
                                   gens_by_bus=gens_by_bus,
                                   bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                   **p_rhs_kwargs
                                   )

    ### declare the generator cost objective
    libgen.declare_expression_pgqg_operating_cost(model=model,
                                                  index_set=gen_attrs['names'],
                                                  p_costs=gen_attrs['p_cost']
                                                  )

    obj_expr = sum(model.pg_operating_cost[gen_name] for gen_name in model.pg_operating_cost)
    if include_feasibility_slack:
        obj_expr += penalty_expr

    model.obj = pe.Objective(expr=obj_expr)

    return model, md


def solve_copperplate_dispatch(model_data,
                solver,
                timelimit = None,
                solver_tee = True,
                symbolic_solver_labels = False,
                options = None,
                copperplate_dispatch_model_generator = create_copperplate_dispatch_approx_model,
                return_model = False,
                return_results = False,
                **kwargs):
    '''
    Create and solve a new copperplate dispatch model

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    solver : str or pyomo.opt.base.solvers.OptSolver
        Either a string specifying a pyomo solver name, or an instantiated pyomo solver
    timelimit : float (optional)
        Time limit for dcopf run. Default of None results in no time
        limit being set.
    solver_tee : bool (optional)
        Display solver log. Default is True.
    symbolic_solver_labels : bool (optional)
        Use symbolic solver labels. Useful for debugging; default is False.
    options : dict (optional)
        Other options to pass into the solver. Default is dict().
    copperplate_dispatch_model_generator : function (optional)
        Function for generating the copperplate dispatch model. Default is
        egret.models.copperplate_dispatch.create_copperplate_dispatch_approx_model
    return_model : bool (optional)
        If True, returns the pyomo model object
    return_results : bool (optional)
        If True, returns the pyomo results object
    kwargs : dictionary (optional)
        Additional arguments for building model
    '''

    import pyomo.environ as pe
    from pyomo.environ import value
    from egret.common.solver_interface import _solve_model
    from egret.model_library.transmission.tx_utils import \
        scale_ModelData_to_pu, unscale_ModelData_to_pu

    m, md = copperplate_dispatch_model_generator(model_data, **kwargs)

    m.dual = pe.Suffix(direction=pe.Suffix.IMPORT)

    m, results = _solve_model(m,solver,timelimit=timelimit,solver_tee=solver_tee,
                              symbolic_solver_labels=symbolic_solver_labels,solver_options=options)

    md = model_data.clone_in_service()

    # save results data to ModelData object
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))

    md.data['system']['total_cost'] = value(m.obj)
    if hasattr(m, 'p_load_shed'):
        md.data['system']['p_balance_violation'] = value(m.p_load_shed) - value(m.p_over_generation)

    for g,g_dict in gens.items():
        g_dict['pg'] = value(m.pg[g])

    for b,b_dict in buses.items():
        b_dict['pl'] = value(m.pl[b])
        b_dict['lmp'] = value(m.dual[m.eq_p_balance])

    unscale_ModelData_to_pu(md, inplace=True)

    if return_model and return_results:
        return md, m, results
    elif return_model:
        return md, m
    elif return_results:
        return md, results
    return md


# if __name__ == '__main__':
#     import os
#     from egret.parsers.matpower_parser import create_ModelData
#
#     path = os.path.dirname(__file__)
#     filename = 'pglib_opf_case3_lmbd.m'
#     matpower_file = os.path.join(path, '../../download/pglib-opf/', filename)
#     md = create_ModelData(matpower_file)
#     kwargs = {'include_feasibility_slack': 'True'}
#     md = solve_copperplate_dispatch(md, "gurobi", **kwargs)
