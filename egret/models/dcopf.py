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
                                               _validate_and_extract_slack_penalty,
                                               create_copperplate_dispatch_approx_model)
from egret.common.log import logger
from math import pi, radians


def _include_feasibility_slack(model, bus_names, bus_p_loads, gens_by_bus, gen_attrs, p_marginal_slack_penalty):
    import egret.model_library.decl as decl

    load_shed_bounds  = {k: (0, tx_utils.load_shed_limit(bus_p_loads[k], gens_by_bus[k], gen_attrs['p_min'])) for k in bus_names}
    decl.declare_var('p_load_shed', model=model, index_set=bus_names,
                     initialize=0., bounds=load_shed_bounds
                     )
    over_gen_bounds = {k: (0, tx_utils.over_gen_limit(bus_p_loads[k], gens_by_bus[k], gen_attrs['p_max'])) for k in bus_names}
    decl.declare_var('p_over_generation', model=model, index_set=bus_names,
                     initialize=0., bounds=over_gen_bounds
                     )

    p_rhs_kwargs = {'include_feasibility_load_shed':'p_load_shed', 'include_feasibility_over_generation':'p_over_generation'}

    penalty_expr = sum(p_marginal_slack_penalty * (model.p_load_shed[bus_name] + model.p_over_generation[bus_name])
                    for bus_name in bus_names)
    return p_rhs_kwargs, penalty_expr

def create_btheta_dcopf_model(model_data, include_angle_diff_limits=False, include_feasibility_slack=False):
    md = model_data.clone_in_service()
    tx_utils.scale_ModelData_to_pu(md, inplace = True)

    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    loads = dict(md.elements(element_type='load'))
    shunts = dict(md.elements(element_type='shunt'))

    dc_branches = dict(md.elements(element_type='dc_branch'))

    gen_attrs = md.attributes(element_type='generator')
    bus_attrs = md.attributes(element_type='bus')
    branch_attrs = md.attributes(element_type='branch')

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

    ### declare the polar voltages
    va_bounds = {k: (-pi, pi) for k in bus_attrs['va']}
    libbus.declare_var_va(model, bus_attrs['names'], initialize=bus_attrs['va'],
                          bounds=va_bounds
                          )

    ### include the feasibility slack for the bus balances
    p_rhs_kwargs = {}
    penalty_expr = None
    if include_feasibility_slack:
        p_marginal_slack_penalty = _validate_and_extract_slack_penalty(md)        
        p_rhs_kwargs, penalty_expr = _include_feasibility_slack(model, bus_attrs['names'], bus_p_loads,
                                                                gens_by_bus, gen_attrs, p_marginal_slack_penalty)

    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    ref_angle = md.data['system']['reference_bus_angle']
    model.va[ref_bus].fix(radians(ref_angle))

    ### declare the generator real power
    pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model, gen_attrs['names'], initialize=pg_init,
                          bounds=zip_items(gen_attrs['p_min'], gen_attrs['p_max'])
                          )

    ### declare the current flows in the branches
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    vj_init = {k: bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    p_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
    p_lbub = dict()
    for k in branches.keys():
        k_pmax = p_max[k]
        if k_pmax is None:
            p_lbub[k] = (None, None)
        else:
            p_lbub[k] = (-k_pmax,k_pmax)
    pf_bounds = p_lbub
    pf_init = dict()
    for branch_name, branch in branches.items():
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']
        y_matrix = tx_calc.calculate_y_matrix_from_branch(branch)
        ifr_init = tx_calc.calculate_ifr(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        ifj_init = tx_calc.calculate_ifj(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        pf_init[branch_name] = tx_calc.calculate_p(ifr_init, ifj_init, vr_init[from_bus], vj_init[from_bus])

    libbranch.declare_var_pf(model=model,
                             index_set=branch_attrs['names'],
                             initialize=pf_init,
                             bounds=pf_bounds
                             )

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


    ### declare the branch power flow approximation constraints
    libbranch.declare_eq_branch_power_btheta_approx(model=model,
                                                    index_set=branch_attrs['names'],
                                                    branches=branches
                                                    )

    ### declare the p balance
    libbus.declare_eq_p_balance_dc_approx(model=model,
                                          index_set=bus_attrs['names'],
                                          bus_p_loads=bus_p_loads,
                                          gens_by_bus=gens_by_bus,
                                          bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                          inlet_branches_by_bus=inlet_branches_by_bus,
                                          outlet_branches_by_bus=outlet_branches_by_bus,
                                          approximation_type=ApproximationType.BTHETA,
                                          dc_inlet_branches_by_bus=dc_inlet_branches_by_bus,
                                          dc_outlet_branches_by_bus=dc_outlet_branches_by_bus,
                                          **p_rhs_kwargs
                                          )

    ### declare the real power flow limits
    libbranch.declare_ineq_p_branch_thermal_lbub(model=model,
                                                 index_set=branch_attrs['names'],
                                                 branches=branches,
                                                 p_thermal_limits=p_max,
                                                 approximation_type=ApproximationType.BTHETA
                                                 )

    ### declare angle difference limits on interconnected buses
    if include_angle_diff_limits:
        libbranch.declare_ineq_angle_diff_branch_lbub(model=model,
                                                      index_set=branch_attrs['names'],
                                                      branches=branches,
                                                      coordinate_type=CoordinateType.POLAR
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

def create_ptdf_dcopf_model(model_data, include_feasibility_slack=False, base_point=BasePointType.FLATSTART, ptdf_options=None):

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

    ## Do and store PTDF calculation
    reference_bus = md.data['system']['reference_bus']

    PTDF = ptdf_utils.VirtualPTDFMatrix(branches, buses, reference_bus, base_point, ptdf_options, branches_keys=branches_idx, buses_keys=buses_idx)

    model._PTDF = PTDF
    model._ptdf_options = ptdf_options

    if ptdf_options['lazy']:

        ### add "blank" real power flow limits
        libbranch.declare_ineq_p_branch_thermal_bounds(model=model,
                                                       index_set=branches_idx,
                                                       branches=branches,
                                                       p_thermal_limits=None,
                                                       approximation_type=None,
                                                       )

        ### add helpers for tracking monitored branches
        lpu.add_monitored_flow_tracker(model)

        ### add initial branches to monitored set
        lpu.add_initial_monitored_branches(model, branches, branches_idx, ptdf_options, PTDF)

    else:
        p_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
        ## add all the constraints
        ### declare the branch power flow approximation constraints
        libbranch.declare_eq_branch_power_ptdf_approx(model=model,
                                                      index_set=branches_idx,
                                                      PTDF=PTDF,
                                                      abs_ptdf_tol=ptdf_options['abs_ptdf_tol'],
                                                      rel_ptdf_tol=ptdf_options['rel_ptdf_tol'],
                                                      )

        ### add all the limits
        libbranch.declare_ineq_p_branch_thermal_lbub(model=model,
                                                     index_set=branches_idx,
                                                     branches=branches,
                                                     p_thermal_limits=p_max,
                                                     approximation_type=ApproximationType.PTDF,
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
    from egret.common.solver_interface import _solve_model

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
            m, results, solver = _solve_model(m, solver, solver_tee=solver_tee, return_solver=True, vars_to_load=[], set_instance=False)
            solver.load_vars()
        else:
            m, results, solver = _solve_model(m, solver, solver_tee=solver_tee, return_solver=True)

    else: # we hit the iteration limit
        logger.warning('WARNING: Exiting on maximum iterations for lazy PTDF model. Result is not transmission feasible.')
        if persistent_solver:
            solver.load_duals()
        return lpu.LazyPTDFTerminationCondition.ITERATION_LIMIT


def solve_dcopf(model_data,
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

def check_instance_feasibility(instance, tolerance, active_only=True):
    infeasibilities = list()

    for con in instance.component_data_objects(pe.Constraint, descend_into=True, sort=True):
        if active_only == False or con.active == True:
            resid = compute_constraint_resid(con)
            if (resid > tolerance):
                infeasibilities.append(constraint_resid_to_string(con.getname(True), con, resid))

    for var in instance.component_data_objects(pe.Var, descend_into=True, sort=True):
        lb = var.lb
        ub = var.ub

        if (ub is not None and lb is not None and ub < lb):
            infeasibility_found = True
            infeasibilities.append('Var: {0} has an upper bound ({1}) that is smaller than its lower bound ({2})'.format(
                var.getname(True), ub, lb))
        if (ub is not None and value(var) > ub):
            infeasibility_found = True
            infeasibilities.append('Var: {0} has an value ({1} that is greater than its upper bound ({2})'.format(
                var.getname(True), value(var), ub))
        if (lb and value(var) < lb):
            infeasibility_found = True
            infeasibilities.append('Var: {0} has an value ({1}) that is less than its lower bound ({2})'.format(
                var.getname(True), value(var), lb))

    if len(infeasibilities) > 0:
        print("*** Infeasibilities found in check_instance_feasibility")
        for s in infeasibilities:
            print(s)
        print("***")

    return len(infeasibilities) == 0

def compute_constraint_resid(con):
    bodyval = value(con.body)
    upper_resid = 0
    if con.upper is not None:
        upper_resid = max(0, bodyval - value(con.upper))
    lower_resid = 0
    if con.lower is not None:
        lower_resid = max(0, value(con.lower) - bodyval)
    return  max(upper_resid, lower_resid)

def constraint_resid_to_string(name, con, resid):
    if con.lower is None and con.upper is None:
        return '{0:10.4g} | {2:10s} <= {3:10.4g} <= {4:10s} : {1}'.format(resid, name, '-', value(con.body), '-')
    elif con.lower is None:
        return '{0:10.4g} | {2:10s} <= {3:10.4g} <= {4:10.4g} : {1}'.format(resid, name, '-', value(con.body), value(con.upper))
    elif con.upper is None:
        return '{0:10.4g} | {2:10.4} <= {3:10.4g} <= {4:10s} : {1}'.format(resid, name, value(con.lower), value(con.body), '-')
    else:
        return '{0:10.4g} | {2:10.4} <= {3:10.4g} <= {4:10.4g} : {1}'.format(resid, name, value(con.lower), value(con.body), value(con.upper))


# if __name__ == '__main__':
#     import os
#     from egret.parsers.matpower_parser import create_ModelData
#
#     path = os.path.dirname(__file__)
#     print(path)
#     filename = 'pglib_opf_case30_ieee.m'
#     test_case = os.path.join(path, '../../download/pglib-opf-master/', filename)
#     md_dict = create_ModelData(test_case)
#
#     dcopf_model = create_ptdf_dcopf_model
#
#     kwargs = {'ptdf_options': {'save_to': test_case + '.pickle'}}
#     md_serialization, results = solve_dcopf(md_dict, "ipopt", dcopf_model_generator=dcopf_model, solver_tee=False,
#                                             return_results=True, **kwargs)
#
#     kwargs = {'ptdf_options': {'load_from': test_case + '.pickle'}}
#     md_deserialization, results = solve_dcopf(md_dict, "ipopt", dcopf_model_generator=dcopf_model, solver_tee=False,
#                                               return_results=True, **kwargs)

