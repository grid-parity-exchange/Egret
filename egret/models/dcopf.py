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
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen

import egret.data.data_utils as data_utils
from egret.model_library.defn import CoordinateType, ApproximationType, BasePointType
from egret.data.model_data import map_items, zip_items
from egret.models.copperplate_dispatch import _include_system_feasibility_slack, create_copperplate_dispatch_approx_model
from math import pi

import numpy as np

def _include_feasibility_slack(model, bus_attrs, gen_attrs, bus_p_loads, penalty=1000):
    import egret.model_library.decl as decl
    slack_init = {k: 0 for k in bus_attrs['names']}
    slack_bounds = {k: (0, sum(bus_p_loads.values())) for k in bus_attrs['names']}
    decl.declare_var('p_slack_pos', model=model, index_set=bus_attrs['names'],
                     initialize=slack_init, bounds=slack_bounds
                     )
    decl.declare_var('p_slack_neg', model=model, index_set=bus_attrs['names'],
                     initialize=slack_init, bounds=slack_bounds
                     )
    p_rhs_kwargs = {'include_feasibility_slack_pos':'p_slack_pos','include_feasibility_slack_neg':'p_slack_neg'}

    p_penalty = penalty * (max([gen_attrs['p_cost'][k]['values'][1] for k in gen_attrs['names']]) + 1)

    penalty_expr = sum(p_penalty * (model.p_slack_pos[bus_name] + model.p_slack_neg[bus_name])
                    for bus_name in bus_attrs['names'])
    return p_rhs_kwargs, penalty_expr


def create_btheta_dcopf_model(model_data, include_angle_diff_limits=False, include_feasibility_slack=False):
    md = model_data.clone_in_service()
    tx_utils.scale_ModelData_to_pu(md, inplace = True)

    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    loads = dict(md.elements(element_type='load'))
    shunts = dict(md.elements(element_type='shunt'))

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
        p_rhs_kwargs, penalty_expr = _include_feasibility_slack(model, bus_attrs, gen_attrs, bus_p_loads)

    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    model.va[ref_bus].fix(0.0)

    ref_angle = md.data['system']['reference_bus_angle']
    if ref_angle != 0.0:
        raise ValueError('The BTHETA DCOPF formulation currently only supports'
                         ' a reference bus angle of 0 degrees, but an angle'
                         ' of {} degrees was found.'.format(ref_angle))

    ### declare the generator real power
    pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model, gen_attrs['names'], initialize=pg_init,
                          bounds=zip_items(gen_attrs['p_min'], gen_attrs['p_max'])
                          )

    ### declare the current flows in the branches
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    vj_init = {k: bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    p_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
    p_lbub = {k: (-p_max[k],p_max[k]) for k in branches.keys()}
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


def create_ptdf_dcopf_model(model_data, include_feasibility_slack=False):
    md = model_data.clone_in_service()
    tx_utils.scale_ModelData_to_pu(md, inplace = True)

    data_utils.create_dicts_of_ptdf(md)

    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    loads = dict(md.elements(element_type='load'))
    shunts = dict(md.elements(element_type='shunt'))

    gen_attrs = md.attributes(element_type='generator')
    bus_attrs = md.attributes(element_type='bus')
    branch_attrs = md.attributes(element_type='branch')
    load_attrs = md.attributes(element_type='load')
    shunt_attrs = md.attributes(element_type='shunt')

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
        p_rhs_kwargs, penalty_expr = _include_system_feasibility_slack(model, gen_attrs, bus_p_loads)

    ### declare the current flows in the branches
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    vj_init = {k: bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    p_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
    p_lbub = {k: (-p_max[k],p_max[k]) for k in branches.keys()}
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

    ### declare the branch power flow approximation constraints
    libbranch.declare_eq_branch_power_ptdf_approx(model=model,
                                                  index_set=branch_attrs['names'],
                                                  branches=branches,
                                                  bus_p_loads=bus_p_loads,
                                                  gens_by_bus=gens_by_bus,
                                                  bus_gs_fixed_shunts=bus_gs_fixed_shunts
                                                  )

    ### declare the p balance
    libbus.declare_eq_p_balance_ed(model=model,
                                   index_set=bus_attrs['names'],
                                   bus_p_loads=bus_p_loads,
                                   gens_by_bus=gens_by_bus,
                                   bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                   **p_rhs_kwargs
                                   )

    ### declare the real power flow limits
    libbranch.declare_ineq_p_branch_thermal_lbub(model=model,
                                                 index_set=branch_attrs['names'],
                                                 branches=branches,
                                                 p_thermal_limits=p_max,
                                                 approximation_type=ApproximationType.PTDF
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

def create_lazy_ptdf_dcopf_model(model_data, include_feasibility_slack=False):
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
        p_rhs_kwargs, penalty_expr = _include_system_feasibility_slack(model, gen_attrs, bus_p_loads)

    ### declare the p balance
    libbus.declare_eq_p_balance_ed(model=model,
                                   index_set=bus_attrs['names'],
                                   bus_p_loads=bus_p_loads,
                                   gens_by_bus=gens_by_bus,
                                   bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                   **p_rhs_kwargs
                                   )

    ### add "blank" power flow expressions
    libbranch.declare_expr_pf(model=model,
                             index_set=branches.keys(),
                             )
    

    ### add "blank" real power flow limits
    libbranch.declare_ineq_p_branch_thermal_lbub(model=model,
                                                 index_set=branches.keys(),
                                                 branches=branches,
                                                 p_thermal_limits=None,
                                                 approximation_type=None,
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

    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))

    ## to keep things in order
    buses_idx = list(buses.keys())
    branches_idx = list(branches.keys())

    reference_bus = md.data['system']['reference_bus']

    ## calculate PTDFs, but don't do anything with them (for now..)
    #from pyutilib.misc.timing import TicTocTimer
    #timer = TicTocTimer()
    #timer.tic('starting PTDF calculation')
    PTDFM = tx_calc.calculate_ptdf(branches,buses,branches_idx,buses_idx,reference_bus,BasePointType.FLATSTART)
    #timer.toc('done')

    ## store some information we'll need when iterating on the model object
    model._PTDFM = PTDFM
    model._PTDF_bus_idx = buses_idx
    model._PTDF_branch_idx = branches_idx

    model._PTDF_bus_nw_exprs = [ model.pl[bus] + bus_gs_fixed_shunts[bus] - \
                     sum(model.pg[g] for g in gens_by_bus[bus]) for bus in buses_idx]

    model._PTDF_branch_limits = np.array([ branches[branch]['rating_long_term'] for branch in branches_idx ])

    model._PTDF_branches = branches
    model._PTDF_bus_p_loads = bus_p_loads
    model._PTDF_gens_by_bus = gens_by_bus
    model._PTDF_bus_gs_fixed_shunts = bus_gs_fixed_shunts

    return model, md

def _lazy_ptdf_dcopf_model_solve_loop(m, md, solver, timelimit, solver_tee=True, symbolic_solver_labels=True,flow_vio_tol=1.e-5,iteration_limit=100000):
    from pyomo.solvers.plugins.solvers.persistent_solver import PersistentSolver

    buses_idx = m._PTDF_bus_idx
    branches_idx = m._PTDF_branch_idx
    bus_nw_exprs = m._PTDF_bus_nw_exprs
    branch_limits = m._PTDF_branch_limits
    PTDFM = m._PTDFM
    
    branches = m._PTDF_branches
    bus_p_loads = m._PTDF_bus_p_loads
    gens_by_bus = m._PTDF_gens_by_bus
    bus_gs_fixed_shunts = m._PTDF_bus_gs_fixed_shunts

    def _iter_over_viol_set(viol_set):
        for i in viol_set:
            bn = branches_idx[i]
            branch = branches[bn]
            if 'ptdf' in branch:
                ## This case should be rare, but it could happen that
                ## we've already added the lower (upper) constraint
                ## also need to add the upper (lower) constraint
                pass
            else:
                ## add the ptdf row for this branch, and the flow-tracking expression
                branch['ptdf'] = {bus : PTDFM[i,j] for j, bus in enumerate(buses_idx)}
                expr = libbranch.get_power_flow_expr_ptdf_approx(m, branch, bus_p_loads, gens_by_bus, bus_gs_fixed_shunts, ptdf_tol=None, approximation_type=ApproximationType.PTDF)
                m.pf[bn] = expr
            yield bn, branch

    def _sanity_checks(bn, thermal_limit, constr, ub):
        if thermal_limit is None:
            raise Exception("Branch was found to be violated, but no thermal limits exist, branch={}".format(bn))
        if bn in constr:
            if ub:
                direction, constrn = 'above', 'upper'
            else:
                direction, constrn = 'below', 'lower'
            raise Exception("Branch was found to be violated from {0}, but a {1} bound constraint is already on the model, branch={2}".format(direction, constrn, bn))

    persistent_solver = isinstance(solver, PersistentSolver)

    for i in range(iteration_limit):

        NWV = np.array([pe.value(bus_nw_expr) for bus_nw_expr in bus_nw_exprs])

        PFV  = np.dot(PTDFM, NWV)

        ## get the indices of the violations, but do it in numpy
        gt_viol = np.nonzero(np.greater(PFV, branch_limits+flow_vio_tol))[0]
        lt_viol = np.nonzero(np.less(PFV, -branch_limits-flow_vio_tol))[0]

        viol_num = len(gt_viol)+len(lt_viol)
        print("iteration {0}, found {1} violation(s)".format(i,viol_num))

        if viol_num <= 0:
            ## in this case, there are no violations!
            ## load the duals now too, if we're using a persistent solver
            if persistent_solver:
                solver.load_duals()
            break

        for bn, branch in _iter_over_viol_set(lt_viol):
            constr = m.ineq_pf_branch_thermal_lb
            thermal_limit = branch['rating_long_term']
            _sanity_checks(bn, thermal_limit, constr, False)
            constr[bn] = (-branch['rating_long_term'], m.pf[bn], None)

            if persistent_solver:
                solver.add_constraint(constr[bn])

        for bn, branch in _iter_over_viol_set(gt_viol):
            constr = m.ineq_pf_branch_thermal_ub
            thermal_limit = branch['rating_long_term']
            _sanity_checks(bn, thermal_limit, constr, True)
            constr[bn] = (None, m.pf[bn], branch['rating_long_term'])

            if persistent_solver:
                solver.add_constraint(constr[bn])

        #m.ineq_pf_branch_thermal_lb.pprint()
        #m.ineq_pf_branch_thermal_ub.pprint()

        if persistent_solver:
            solver.solve(m, tee=solver_tee, load_solutions=False, save_results=False)
            solver.load_vars()
        else:
            solver.solve(m, tee=solver_tee, symbolic_solver_labels=symbolic_solver_labels)
    else:
        print('WARNING: Exiting on maximum iterations for lazy PTDF model.')
        print('         Result is not transmission feasible.')


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

    ## we need to create the solver object here for lazy_ptdf
    if isinstance(solver, str):
        solver = po.SolverFactory(solver)
    elif isinstance(solver, po.base.OptSolver):
        pass
    else:
        raise Exception('solver must be string or an instanciated pyomo solver')

    m, results = _solve_model(m,solver,timelimit=timelimit,solver_tee=solver_tee,
                              symbolic_solver_labels=symbolic_solver_labels,options=options)

    if dcopf_model_generator == create_lazy_ptdf_dcopf_model:
        _lazy_ptdf_dcopf_model_solve_loop(m, md, solver, timelimit=timelimit, solver_tee=solver_tee, symbolic_solver_labels=symbolic_solver_labels,iteration_limit=100)

    # save results data to ModelData object
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))

    md.data['system']['total_cost'] = value(m.obj)

    for g,g_dict in gens.items():
        g_dict['pg'] = value(m.pg[g])

    for b,b_dict in buses.items():
        b_dict['pl'] = value(m.pl[b])
        if dcopf_model_generator == create_btheta_dcopf_model:
            b_dict['lmp'] = value(m.dual[m.eq_p_balance[b]])
            b_dict['va'] = value(m.va[b])
        elif dcopf_model_generator == create_ptdf_dcopf_model:
            b_dict['lmp'] = value(m.dual[m.eq_p_balance])
            for k, k_dict in branches.items():
                if k_dict['from_bus'] == b or k_dict['to_bus'] == b:
                    b_dict['lmp'] += k_dict['ptdf'][b]*value(m.dual[m.eq_pf_branch[k]])
        elif dcopf_model_generator == create_lazy_ptdf_dcopf_model:
            b_dict['lmp'] = value(m.dual[m.eq_p_balance])
            for k, k_dict in branches.items():
                if k_dict['from_bus'] == b or k_dict['to_bus'] == b:
                    ## NOTE: if line k's thermal limits are not binding,
                    ##       its dual value is 0.
                    ## NOTE: These need to be - for the dual to be correct.
                    ##       Is this because our ptdf's are "backwards"?
                    if k in m.ineq_pf_branch_thermal_lb:
                        b_dict['lmp'] -= k_dict['ptdf'][b]*value(m.dual[m.ineq_pf_branch_thermal_lb[k]])
                    if k in m.ineq_pf_branch_thermal_ub:
                        b_dict['lmp'] -= k_dict['ptdf'][b]*value(m.dual[m.ineq_pf_branch_thermal_ub[k]])
        else:
            raise Exception("Unrecognized dcopf_mode_generator {}".format(dcopf_model_generator))

    ## calculate the power flows from our PTDF matrix for maximum precision
    if dcopf_model_generator == create_lazy_ptdf_dcopf_model:
        bus_nw_exprs = m._PTDF_bus_nw_exprs
        PTDFM = m._PTDFM
        branches_idx = m._PTDF_branch_idx

        NWV = np.array([pe.value(bus_nw_expr) for bus_nw_expr in bus_nw_exprs])
        PFV  = np.dot(PTDFM, NWV)

        for i,bn in enumerate(branches_idx):
            branches[bn]['pf'] = PFV[i]

    else:
        for k, k_dict in branches.items():
            k_dict['pf'] = value(m.pf[k])

    unscale_ModelData_to_pu(md, inplace=True)

    if return_model and return_results:
        return md, m, results
    elif return_model:
        return md, m
    elif return_results:
        return md, results
    return md


if __name__ == '__main__':
    import os
    from egret.parsers.matpower_parser import create_ModelData
    from pyutilib.misc.timing import TicTocTimer

    timer = TicTocTimer()
    timer.tic('loading instance model_data')
    path = os.path.dirname(__file__)
    #filename = 'pglib_opf_case1354_pegase.m'
    #filename = 'pglib_opf_case240_pserc.m'
    #filename = 'pglib_opf_case300_ieee.m'
    filename = 'pglib_opf_case118_ieee.m'
    #filename = 'pglib_opf_case30_ieee.m'
    filename = 'pglib_opf_case5_pjm.m'
    matpower_file = os.path.join(path, '../../download/pglib-opf/', filename)
    md = create_ModelData(matpower_file)
    timer.toc('data_loaded')

    kwargs = {'include_feasibility_slack':'True'}
    #kwargs = {}
    solver='gurobi_persistent'
    #options = {'optimality_tol':1e-09, 'feasibility_tol':1e-09}
    options = {'method':1}
    #options = {}

    timer.tic('solving btheta DCOPF using '+solver)
    mdo_bt = solve_dcopf(md, solver, options=options, dcopf_model_generator=create_btheta_dcopf_model, **kwargs)
    timer.toc('solved btheta DCOPF using '+solver)

    timer.tic('solving full PTDF DCOPF using '+solver)
    mdo_fp = solve_dcopf(md, solver, options=options, dcopf_model_generator=create_ptdf_dcopf_model, **kwargs)
    timer.toc('solved full PTDF DCOPF using '+solver)

    timer.tic('solving lazy PTDF DCOPF using '+solver)
    mdo_lp = solve_dcopf(md, solver, options=options, dcopf_model_generator=create_lazy_ptdf_dcopf_model, **kwargs)
    timer.toc('solved lazy PTDF DCOPF using '+solver)
