#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module provides functions that create the modules for typical DCOPF formulations.

#TODO: document this with examples
"""
import pyomo.environ as pe
import numpy as np
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen
import egret.common.lazy_ptdf_utils as lpu
import egret.data.ptdf_utils as ptdf_utils

from egret.model_library.defn import CoordinateType, ApproximationType, BasePointType
from egret.data.data_utils import map_items, zip_items
from egret.models.copperplate_dispatch import (_include_system_feasibility_slack,
                                               _validate_and_extract_slack_penalty)
from egret.common.log import logger
from math import pi, radians


def create_scopf_model(model_data, include_feasibility_slack=False, base_point=BasePointType.FLATSTART, ptdf_options=None):

    ptdf_options = lpu.populate_default_ptdf_options(ptdf_options)

    baseMVA = model_data.data['system']['baseMVA']
    lpu.check_and_scale_ptdf_options(ptdf_options, baseMVA)
    
    md = model_data.clone_in_service()
    tx_utils.scale_ModelData_to_pu(md, inplace = True)


    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    loads = dict(md.elements(element_type='load'))
    shunts = dict(md.elements(element_type='shunt'))

    dc_branches = dict(md.elements(element_type='dc_branch'))
    contingencies = dict(md.elements(element_type='contingency'))

    gen_attrs = md.attributes(element_type='generator')
    ## to keep things in order
    buses_idx = tuple(buses.keys())
    branches_idx = tuple(branches.keys())

    inlet_branches_by_bus, outlet_branches_by_bus = \
        tx_utils.inlet_outlet_branches_by_bus(branches, buses)
    gens_by_bus = tx_utils.gens_by_bus(buses, gens)

    model = pe.ConcreteModel()

    ### declare (and fix) the loads at the buses
    bus_p_loads, _ = tx_utils.dict_of_bus_loads(buses, loads)

    libbus.declare_var_pl(model, buses_idx, initialize=bus_p_loads)
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
        p_marginal_slack_penalty = _validate_and_extract_slack_penalty(md)                
        p_rhs_kwargs, penalty_expr = _include_system_feasibility_slack(model, bus_p_loads, gen_attrs, p_marginal_slack_penalty)

    if dc_branches:
        dcpf_bounds = dict()
        for k, k_dict in dc_branches.items():
            kp_max = k_dict['rating_long_term']
            if kp_max is None:
                dcpf_bounds[k] = (None, None)
            else:
                dcpf_bounds[k] = (-kp_max, kp_max)
        libbranch.declare_var_dcpf(model=model,
                                   index_set=dc_branches.keys(),
                                   initialize=0.,
                                   bounds=dcpf_bounds,
                                  )
        dc_inlet_branches_by_bus, dc_outlet_branches_by_bus = \
                tx_utils.inlet_outlet_branches_by_bus(dc_branches, buses)
    else:
        dc_inlet_branches_by_bus = None
        dc_outlet_branches_by_bus = None

    ### declare the p balance
    libbus.declare_eq_p_balance_ed(model=model,
                                   index_set=buses_idx,
                                   bus_p_loads=bus_p_loads,
                                   gens_by_bus=gens_by_bus,
                                   bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                   **p_rhs_kwargs
                                   )

    ### declare net withdraw expression for use in PTDF power flows
    libbus.declare_expr_p_net_withdraw_at_bus(model=model,
                                              index_set=buses_idx,
                                              bus_p_loads=bus_p_loads,
                                              gens_by_bus=gens_by_bus,
                                              bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                              dc_inlet_branches_by_bus=dc_inlet_branches_by_bus,
                                              dc_outlet_branches_by_bus=dc_outlet_branches_by_bus,
                                              )
    
    ### add "blank" power flow expressions
    libbranch.declare_expr_pf(model=model,
                              index_set=branches_idx,
                              )

    ### add "blank" power flow expressions
    libbranch.declare_expr_pfc(model=model,
                               index_set=(contingencies.keys(), branches_idx),
                               )

    ## Do and store PTDF calculation
    reference_bus = md.data['system']['reference_bus']

    PTDF = ptdf_utils.VirtualPTDFMatrix(branches, buses, reference_bus, base_point, ptdf_options,\
                                        contingencies=contingencies, branches_keys=branches_idx, buses_keys=buses_idx)

    model._PTDF = PTDF
    model._ptdf_options = ptdf_options

    if not ptdf_options['lazy']:
        raise RuntimeError("scopf only supports lazy constraint generation")

    ### add "blank" real power flow limits
    libbranch.declare_ineq_p_branch_thermal_bounds(model=model,
                                                   index_set=branches_idx,
                                                   branches=branches,
                                                   p_thermal_limits=None,
                                                   approximation_type=None,
                                                   )

    ### add "blank" real power flow limits
    libbranch.declare_ineq_p_contingency_branch_thermal_bounds(model=model,
                                                               index_set=(contingencies.keys(), branches_idx),
                                                               pc_thermal_limits=None,
                                                               approximation_type=None,
                                                               )

    ### add helpers for tracking monitored branches
    lpu.add_monitored_flow_tracker(model)

    ### add initial branches to monitored set
    lpu.add_initial_monitored_branches(model, branches, branches_idx, ptdf_options, PTDF)

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

def _lazy_ptdf_dcopf_model_solve_loop(m, md, solver, solver_tee=True, symbolic_solver_labels=False, iteration_limit=100000):
    '''
    The lazy PTDF DCOPF solver loop. This function iteratively
    adds violated transmission constraints until either the result is
    transmission feasible or we're tracking every violated constraint
    in the model

    Parameters
    ----------
    m : pyomo.environ.ConcreteModel
        An egret DCOPF model with no transmission constraints
    md : egret.data.ModelData
        An egret ModelData object
    solver : pyomo.opt.solver
        A pyomo solver object
    solver_tee : bool (optional)
        For displaying the solver log (default is True)
    symbolic_solver_labels : bool (optional)
        Use symbolic solver labels when writing to the solver (default is False)
    iteration_limit : int (optional)
        Number of iterations before a hard termination (default is 100000)

    Returns
    -------
    egret.common.lazy_ptdf_utils.LazyPTDFTerminationCondition : the termination status
    pyomo.opt.results.SolverResults : The results object from the pyomo solver
    int : The number of iterations before termination

    '''
    from pyomo.solvers.plugins.solvers.persistent_solver import PersistentSolver

    PTDF = m._PTDF

    ptdf_options = m._ptdf_options

    persistent_solver = isinstance(solver, PersistentSolver)

    for i in range(iteration_limit):

        flows, viol_num, mon_viol_num, viol_lazy \
                = lpu.check_violations(m, md, PTDF, ptdf_options['max_violations_per_iteration'])

        iter_status_str = "iteration {0}, found {1} violation(s)".format(i,viol_num)
        if mon_viol_num:
            iter_status_str += ", {} of which are already monitored".format(mon_viol_num)

        logger.info(iter_status_str)

        if viol_num <= 0:
            ## in this case, there are no violations!
            ## load the duals now too, if we're using a persistent solver
            if persistent_solver:
                solver.load_duals()
            return lpu.LazyPTDFTerminationCondition.NORMAL

        elif viol_num == mon_viol_num:
            logger.warning('WARNING: Terminating with monitored violations! Result is not transmission feasible.')
            if persistent_solver:
                solver.load_duals()
            return lpu.LazyPTDFTerminationCondition.FLOW_VIOLATION

        lpu.add_violations(viol_lazy, flows, m, md, solver, ptdf_options, PTDF)
        total_flow_constr_added = len(viol_lazy)
        logger.info( "iteration {0}, added {1} flow constraint(s)".format(i,total_flow_constr_added))

        if persistent_solver:
            solver.solve(m, tee=solver_tee, load_solutions=False, save_results=False)
            solver.load_vars()
        else:
            solver.solve(m, tee=solver_tee, symbolic_solver_labels=symbolic_solver_labels)

    else: # we hit the iteration limit
        logger.warning('WARNING: Exiting on maximum iterations for lazy PTDF model. Result is not transmission feasible.')
        if persistent_solver:
            solver.load_duals()
        return lpu.LazyPTDFTerminationCondition.ITERATION_LIMIT


def solve_scopf(model_data,
                solver,
                timelimit = None,
                solver_tee = True,
                symbolic_solver_labels = False,
                options = None,
                dcopf_model_generator = create_btheta_dcopf_model,
                return_model = False,
                return_results = False,
                **kwargs):
    '''
    Create and solve a new dcopf model

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
    dcopf_model_generator : function (optional)
        Function for generating the dcopf model. Default is
        egret.models.dcopf.create_btheta_dcopf_model
    return_model : bool (optional)
        If True, returns the pyomo model object
    return_results : bool (optional)
        If True, returns the pyomo results object
    kwargs : dictionary (optional)
        Additional arguments for building model
    '''

    import pyomo.environ as pe
    import pyomo.opt as po
    from pyomo.environ import value
    from egret.common.solver_interface import _solve_model
    from egret.model_library.transmission.tx_utils import \
        scale_ModelData_to_pu, unscale_ModelData_to_pu

    m, md = dcopf_model_generator(model_data, **kwargs)

    m.dual = pe.Suffix(direction=pe.Suffix.IMPORT)

    m, results, solver = _solve_model(m,solver,timelimit=timelimit,solver_tee=solver_tee,
                              symbolic_solver_labels=symbolic_solver_labels,solver_options=options, return_solver=True)

    if dcopf_model_generator == create_ptdf_dcopf_model and m._ptdf_options['lazy']:
        iter_limit = m._ptdf_options['iteration_limit']
        term_cond = _lazy_ptdf_dcopf_model_solve_loop(m, md, solver, solver_tee=solver_tee, symbolic_solver_labels=symbolic_solver_labels,iteration_limit=iter_limit)

    # save results data to ModelData object
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))

    dc_branches = dict(md.elements(element_type='dc_branch'))

    md.data['system']['total_cost'] = value(m.obj)

    for g,g_dict in gens.items():
        g_dict['pg'] = value(m.pg[g])

    ## calculate the power flows from our PTDF matrix for maximum precision
    ## calculate the LMPC (LMP congestion) using numpy
    if dcopf_model_generator == create_ptdf_dcopf_model:
        PTDF = m._PTDF

        PFV, _, VA = PTDF.calculate_PFV(m)

        branches_idx = PTDF.branches_keys
        for i,bn in enumerate(branches_idx):
            branches[bn]['pf'] = PFV[i]
    else:
        for k, k_dict in branches.items():
            k_dict['pf'] = value(m.pf[k])

    if dcopf_model_generator == create_ptdf_dcopf_model:
        if hasattr(m, 'p_load_shed'):
            md.data['system']['p_balance_violation'] = value(m.p_load_shed) - value(m.p_over_generation)
        buses_idx = PTDF.buses_keys
        LMP = PTDF.calculate_LMP(m, m.dual, m.eq_p_balance)
        for i,b in enumerate(buses_idx):
            b_dict = buses[b]
            b_dict['lmp'] = LMP[i]
            b_dict['pl'] = value(m.pl[b])
            b_dict['va'] = VA[i]
    else:
        for b,b_dict in buses.items():
            if hasattr(m, 'p_load_shed'):
                b_dict['p_balance_violation'] = value(m.p_load_shed[b]) - value(m.p_over_generation[b])
            b_dict['pl'] = value(m.pl[b])
            if dcopf_model_generator == create_btheta_dcopf_model:
                b_dict['lmp'] = value(m.dual[m.eq_p_balance[b]])
                b_dict['va'] = value(m.va[b])
            else:
                raise Exception("Unrecognized dcopf_model_generator {}".format(dcopf_model_generator))

    for k, k_dict in dc_branches.items():
        k_dict['pf'] = value(m.dcpf[k])

    unscale_ModelData_to_pu(md, inplace=True)

    if return_model and return_results:
        return md, m, results
    elif return_model:
        return md, m
    elif return_results:
        return md, results
    return md
