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
import egret.data.data_utils as data_utils

from egret.model_library.defn import CoordinateType, ApproximationType, BasePointType
from egret.data.model_data import map_items, zip_items
from egret.models.copperplate_dispatch import _include_system_feasibility_slack, create_copperplate_dispatch_approx_model
from math import pi, radians


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


def create_ptdf_dcopf_model(model_data, include_feasibility_slack=False,base_point=BasePointType.FLATSTART,
                            ptdf_options=None):
    if ptdf_options is None:
        ptdf_options_dict = dict()
    else:
        ptdf_options_dict = ptdf_options

    lpu.populate_default_ptdf_options(ptdf_options_dict)

    baseMVA = model_data.data['system']['baseMVA']
    lpu.check_and_scale_ptdf_options(ptdf_options_dict, baseMVA)

    rel_ptdf_tol = ptdf_options_dict['rel_ptdf_tol']
    abs_ptdf_tol = ptdf_options_dict['abs_ptdf_tol']

    md = model_data.clone_in_service()

    tx_utils.scale_ModelData_to_pu(md, inplace = True)

    data_utils.create_dicts_of_ptdf(md,base_point=base_point)

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
                                                  buses=buses,
                                                  bus_p_loads=bus_p_loads,
                                                  gens_by_bus=gens_by_bus,
                                                  bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                                  abs_ptdf_tol=abs_ptdf_tol,
                                                  rel_ptdf_tol=rel_ptdf_tol,
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


def create_lazy_ptdf_dcopf_model(model_data, include_feasibility_slack=False, ptdf_options=None):
    
    if ptdf_options is None:
        ptdf_options_dict = dict()
    else:
        ptdf_options_dict = ptdf_options

    lpu.populate_default_ptdf_options(ptdf_options_dict)

    baseMVA = model_data.data['system']['baseMVA']
    lpu.check_and_scale_ptdf_options(ptdf_options_dict, baseMVA)
    
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
    phi_from, phi_to = tx_calc.calculate_phi_constant(branches,branches_idx, buses_idx,)

    phi_adjust_array = np.array([phi_from[i].sum()-phi_to[i].sum() for i,_ in enumerate(buses_idx)])

    branch_list = [ branches[bn] for bn in branches_idx ]

    phase_shift_array = np.array([ -(1/branch['reactance']) * (radians(branch['transformer_phase_shift'])/branch['transformer_tap_ratio']) if (branch['branch_type'] == 'transformer') else 0. for branch in branch_list])

    ## store some information we'll need when iterating on the model object
    model._PTDF_dict = {'PTDFM' : PTDFM,
                        'buses_idx': buses_idx,
                        'branches_idx' : branches_idx,
                        'branch_limits' : np.array([ branches[branch]['rating_long_term'] for branch in branches_idx ]),
                        'branches' : branches,
                        'gens_by_bus' : gens_by_bus,
                        'bus_gs_fixed_shunts' : bus_gs_fixed_shunts,
                        'phi_adjust_array': phi_adjust_array,
                        'phase_shift_array': phase_shift_array,
                        }
    model._PTDF_bus_nw_exprs = [ model.pl[bus] + bus_gs_fixed_shunts[bus] \
                    - sum(model.pg[g] for g in gens_by_bus[bus])
                     for i,bus in enumerate(buses_idx)]

    model._PTDF_bus_p_loads = bus_p_loads

    model._ptdf_options_dict = ptdf_options_dict

    return model, md

def _lazy_ptdf_dcopf_model_solve_loop(m, md, solver, timelimit, solver_tee=True, symbolic_solver_labels=False, iteration_limit=100000):
    from pyomo.solvers.plugins.solvers.persistent_solver import PersistentSolver

    PTDF_dict = m._PTDF_dict
    bus_nw_exprs = m._PTDF_bus_nw_exprs
    bus_p_loads = m._PTDF_bus_p_loads

    ptdf_options_dict = m._ptdf_options_dict


    persistent_solver = isinstance(solver, PersistentSolver)

    lazy_rel_flow_tol = ptdf_options_dict['lazy_rel_flow_tol']

    rel_flow_tol = ptdf_options_dict['rel_flow_tol']
    abs_flow_tol = ptdf_options_dict['abs_flow_tol']

    branch_limits = PTDF_dict['branch_limits']

    ## if lazy_rel_flow_tol < 0, this narrows the branch limits
    PTDF_dict['lazy_branch_limits'] = branch_limits*(1+lazy_rel_flow_tol)

    ## only enforce the relative and absolute, within tollerance
    PTDF_dict['enforced_branch_limits'] = np.maximum(branch_limits*(1+rel_flow_tol), branch_limits+abs_flow_tol)

    for i in range(iteration_limit):

        PFV, viol_num, viols_tup = lpu.check_violations(PTDF_dict, bus_nw_exprs)

        print("iteration {0}, found {1} violation(s)".format(i,viol_num))

        if viol_num <= 0:
            ## in this case, there are no violations!
            ## load the duals now too, if we're using a persistent solver
            if persistent_solver:
                solver.load_duals()
            return lpu.LazyPTDFTerminationCondition.NORMAL

        all_viol_in_model = lpu.add_violations(viols_tup, PFV, m, md, solver, ptdf_options_dict, PTDF_dict, bus_nw_exprs, bus_p_loads,)

        if all_viol_in_model:
            print('WARNING: Terminating with monitored violations!')
            print('         Result is not transmission feasible.')
            if persistent_solver:
                solver.load_duals()
            return lpu.LazyPTDFTerminationCondition.FLOW_VIOLATION

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
                              symbolic_solver_labels=symbolic_solver_labels,options=options, return_solver=True)

    if dcopf_model_generator == create_lazy_ptdf_dcopf_model:
        iter_limit = m._ptdf_options_dict['iteration_limit']
        term_cond = _lazy_ptdf_dcopf_model_solve_loop(m, md, solver, timelimit=timelimit, solver_tee=solver_tee, symbolic_solver_labels=symbolic_solver_labels,iteration_limit=iter_limit)

    # save results data to ModelData object
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))

    md.data['system']['total_cost'] = value(m.obj)

    from pyutilib.misc.timing import TicTocTimer

    timer = TicTocTimer()

    for g,g_dict in gens.items():
        g_dict['pg'] = value(m.pg[g])

    ## calculate the power flows from our PTDF matrix for maximum precision
    ## calculate the LMPC (LMP congestion) using numpy
    if dcopf_model_generator == create_lazy_ptdf_dcopf_model:
        PTDF_dict = m._PTDF_dict
        bus_nw_exprs = m._PTDF_bus_nw_exprs
        PTDFM = PTDF_dict['PTDFM']
        branches_idx = PTDF_dict['branches_idx']

        NWV = np.array([pe.value(bus_nw_expr) for bus_nw_expr in bus_nw_exprs])
        PFV  = np.dot(PTDFM, NWV)
        PFD = np.zeros(len(branches_idx))
        for i,bn in enumerate(branches_idx):
            branches[bn]['pf'] = PFV[i]
            if bn in m.ineq_pf_branch_thermal_lb:
                PFD[i] += value(m.dual[m.ineq_pf_branch_thermal_lb[bn]])
            if bn in m.ineq_pf_branch_thermal_ub:
                PFD[i] += value(m.dual[m.ineq_pf_branch_thermal_ub[bn]])
        ## TODO: PFD is likely to be sparse, implying we just need a few
        ##       rows of the PTDF matrix (or columns in its transpose).
        LMPC = np.dot(-PTDFM.T, PFD)
    else:
        for k, k_dict in branches.items():
            k_dict['pf'] = value(m.pf[k])

    if dcopf_model_generator == create_lazy_ptdf_dcopf_model:
        buses_idx = PTDF_dict['buses_idx']
        LMPE = value(m.dual[m.eq_p_balance])
        for i,b in enumerate(buses_idx):
            b_dict = buses[b]
            b_dict['lmp'] = LMPE + LMPC[i]
            b_dict['pl'] = value(m.pl[b])
    else:
        for b,b_dict in buses.items():
            b_dict['pl'] = value(m.pl[b])
            if dcopf_model_generator == create_btheta_dcopf_model:
                b_dict['lmp'] = value(m.dual[m.eq_p_balance[b]])
                b_dict['va'] = value(m.va[b])
            elif dcopf_model_generator == create_ptdf_dcopf_model:
                b_dict['lmp'] = value(m.dual[m.eq_p_balance])
                for k, k_dict in branches.items():
                    b_dict['lmp'] += k_dict['ptdf'][b]*value(m.dual[m.eq_pf_branch[k]])
            else:
                raise Exception("Unrecognized dcopf_mode_generator {}".format(dcopf_model_generator))

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
#     filename = 'pglib_opf_case300_ieee.m'
#     matpower_file = os.path.join(path, '../../download/pglib-opf/', filename)
#     md = create_ModelData(matpower_file)
#
#     kwargs = {'include_feasibility_slack':False}
#     md_btheta, m_btheta, results_btheta = solve_dcopf(md, "gurobi", dcopf_model_generator=create_btheta_dcopf_model, return_model=True, return_results=True, **kwargs)
#
#     from acopf import solve_acopf
#     model_data, model, results = solve_acopf(md, "ipopt", return_model=True, return_results=True)
#     kwargs = {'include_feasibility_slack':False,'base_point':BasePointType.SOLUTION}
#     md_ptdf, m_ptdf, results_ptdf = solve_dcopf(model_data, "gurobi", dcopf_model_generator=create_ptdf_dcopf_model, return_model=True, return_results=True, **kwargs)
#
