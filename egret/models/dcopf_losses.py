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

Note that since the losses model is quadratic, the create_btheta_losses_dcopf_model and
the create_ptdf_losses_dcopf_model are not equivalent; the former is a QCP and the latter is a LP.

#TODO: document this with examples

"""
import pyomo.environ as pe
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen

import egret.data.ptdf_utils as ptdf_utils
import egret.common.lazy_ptdf_utils as lpu

from egret.model_library.defn import CoordinateType, ApproximationType, RelaxationType, BasePointType
from egret.data.data_utils import map_items, zip_items
from egret.models.copperplate_dispatch import _include_system_feasibility_slack, _validate_and_extract_slack_penalty
from egret.models.dcopf import _include_feasibility_slack
from math import pi, radians


def create_btheta_losses_dcopf_model(model_data, relaxation_type=RelaxationType.SOC, include_angle_diff_limits=False, include_feasibility_slack=False):
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

    ### declare the polar voltages
    va_bounds = {k: (-pi, pi) for k in bus_attrs['va']}
    libbus.declare_var_va(model, bus_attrs['names'], initialize=bus_attrs['va'],
                          bounds=va_bounds
                          )

    dva_initialize = {k: 0.0 for k in branch_attrs['names']}
    libbranch.declare_var_dva(model, branch_attrs['names'],
                              initialize=dva_initialize
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
    pf_bounds = {k: (-p_max[k],p_max[k]) for k in branches.keys()}
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
    pfl_bounds = {k: (0,p_max[k]**2) for k in branches.keys()}
    pfl_init = {k: 0 for k in branches.keys()}

    libbranch.declare_var_pf(model=model,
                             index_set=branch_attrs['names'],
                             initialize=pf_init,
                             bounds=pf_bounds
                             )

    libbranch.declare_var_pfl(model=model,
                              index_set=branch_attrs['names'],
                              initialize=pfl_init,
                              bounds=pfl_bounds
                             )

    ### declare the angle difference constraint
    libbranch.declare_eq_branch_dva(model=model,
                                    index_set=branch_attrs['names'],
                                    branches=branches
                                    )

    ### declare the branch power flow approximation constraints
    libbranch.declare_eq_branch_power_btheta_approx(model=model,
                                                    index_set=branch_attrs['names'],
                                                    branches=branches,
                                                    approximation_type=ApproximationType.BTHETA_LOSSES
                                                    )

    ### declare the branch power loss approximation constraints
    libbranch.declare_eq_branch_loss_btheta_approx(model=model,
                                                    index_set=branch_attrs['names'],
                                                    branches=branches,
                                                    relaxation_type=relaxation_type
                                                    )

    ### declare the p balance
    libbus.declare_eq_p_balance_dc_approx(model=model,
                                          index_set=bus_attrs['names'],
                                          bus_p_loads=bus_p_loads,
                                          gens_by_bus=gens_by_bus,
                                          bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                          inlet_branches_by_bus=inlet_branches_by_bus,
                                          outlet_branches_by_bus=outlet_branches_by_bus,
                                          approximation_type=ApproximationType.BTHETA_LOSSES,
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


def create_ptdf_losses_dcopf_model(model_data, include_feasibility_slack=False, ptdf_options=None):

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
        p_marginal_slack_penalty = _validate_and_extract_slack_penalty(md)
        p_rhs_kwargs, penalty_expr = _include_system_feasibility_slack(model, bus_p_loads, gen_attrs, p_marginal_slack_penalty)

    ### declare net withdraw expression for use in PTDF power flows
    libbus.declare_expr_p_net_withdraw_at_bus(model=model,
                                              index_set=bus_attrs['names'],
                                              bus_p_loads=bus_p_loads,
                                              gens_by_bus=gens_by_bus,
                                              bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                              )

    ### declare the current flows in the branches
    p_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
    pfl_bounds = {k: (-p_max[k]**2,p_max[k]**2) for k in branches.keys()}
    pfl_init = {k: 0 for k in branches.keys()}

    ## Do and store PTDF calculation
    reference_bus = md.data['system']['reference_bus']
    ## We'll assume we have a solution to initialize from
    base_point = BasePointType.SOLUTION

    PTDF = ptdf_utils.PTDFLossesMatrix(branches, buses, reference_bus, base_point, ptdf_options)
    model._PTDF = PTDF
    model._ptdf_options = ptdf_options

    libbranch.declare_expr_pf(model=model,
                             index_set=branch_attrs['names'],
                             )

    libbranch.declare_var_pfl(model=model,
                              index_set=branch_attrs['names'],
                              initialize=pfl_init,
                              bounds=pfl_bounds
                             )

    ### declare the branch power flow approximation constraints
    libbranch.declare_eq_branch_power_ptdf_approx(model=model,
                                                  index_set=branch_attrs['names'],
                                                  PTDF=PTDF,
                                                  abs_ptdf_tol=ptdf_options['abs_ptdf_tol'],
                                                  rel_ptdf_tol=ptdf_options['rel_ptdf_tol'],
                                                  )

    ### declare the branch power loss approximation constraints
    libbranch.declare_eq_branch_loss_ptdf_approx(model=model,
                                                 index_set=branch_attrs['names'],
                                                 PTDF=PTDF,
                                                 abs_ptdf_tol=ptdf_options['abs_ptdf_tol'],
                                                 rel_ptdf_tol=ptdf_options['rel_ptdf_tol'],
                                                 )

    ### declare the p balance
    libbus.declare_eq_p_balance_ed(model=model,
                                   index_set=bus_attrs['names'],
                                   bus_p_loads=bus_p_loads,
                                   gens_by_bus=gens_by_bus,
                                   bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                   include_losses=branch_attrs['names'],
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

def solve_dcopf_losses(model_data,
                solver,
                timelimit = None,
                solver_tee = True,
                symbolic_solver_labels = False,
                options = None,
                dcopf_losses_model_generator = create_btheta_losses_dcopf_model,
                return_model = False,
                return_results = False,
                **kwargs):
    '''
    Create and solve a new dcopf with losses model

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
    from pyomo.environ import value
    from egret.common.solver_interface import _solve_model
    from egret.model_library.transmission.tx_utils import \
        scale_ModelData_to_pu, unscale_ModelData_to_pu

    m, md = dcopf_losses_model_generator(model_data, **kwargs)

    m.dual = pe.Suffix(direction=pe.Suffix.IMPORT)

    m, results = _solve_model(m,solver,timelimit=timelimit,solver_tee=solver_tee,
                              symbolic_solver_labels=symbolic_solver_labels,solver_options=options)

    # save results data to ModelData object
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))

    md.data['system']['total_cost'] = value(m.obj)

    for g,g_dict in gens.items():
        g_dict['pg'] = value(m.pg[g])

    if dcopf_losses_model_generator == create_btheta_losses_dcopf_model:
        for b,b_dict in buses.items():
            b_dict['pl'] = value(m.pl[b])
            b_dict.pop('qlmp',None)
            b_dict['lmp'] = value(m.dual[m.eq_p_balance[b]])
            b_dict['va'] = value(m.va[b])
            if hasattr(m, 'p_load_shed'):
                b_dict['p_balance_violation'] = value(m.p_load_shed[b]) - value(m.p_over_generation[b])
    elif dcopf_losses_model_generator == create_ptdf_losses_dcopf_model:
        PTDF = m._PTDF
        ptdf_r = PTDF.PTDFM
        ldf = PTDF.LDF
        buses_idx = PTDF.buses_keys
        branches_idx = PTDF.branches_keys

        if hasattr(m, 'p_load_shed'):
            md.data['system']['p_balance_violation'] = value(m.p_load_shed) - value(m.p_over_generation)

        for j, b in enumerate(buses_idx):
            b_dict = buses[b]
            b_dict['pl'] = value(m.pl[b])
            b_dict.pop('qlmp',None)

            b_dict['lmp'] = value(m.dual[m.eq_p_balance])
            for i, k in enumerate(branches_idx):
                b_dict['lmp'] += ptdf_r[i,j]*value(m.dual[m.ineq_pf_branch_thermal_lb[k]])
                b_dict['lmp'] += ptdf_r[i,j]*value(m.dual[m.ineq_pf_branch_thermal_ub[k]])
                b_dict['lmp'] += ldf[i,j]*value(m.dual[m.eq_pfl_branch[k]])

    else:
        raise Exception("Unrecognized dcopf_losses_model_generator {}".format(dcopf_losses_model_generator))

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


# if __name__ == '__main__':
#     import os
#     from egret.parsers.matpower_parser import create_ModelData
#
#     path = os.path.dirname(__file__)
#     print(path)
#     filename = 'pglib_opf_case30_ieee.m'
#     test_case = os.path.join(path, '../../download/pglib-opf-master/', filename)
#     md_dict = create_ModelData(test_case)
#     dcopf_losses_model = create_ptdf_losses_dcopf_model
#
#     from egret.models.acopf import solve_acopf
#
#     md_dict, _, _ = solve_acopf(md_dict, "ipopt", solver_tee=False, return_model=True, return_results=True)
#
#     kwargs = {}
#     md, results = solve_dcopf_losses(md_dict, "ipopt", dcopf_losses_model_generator=dcopf_losses_model,
#                                      solver_tee=False, return_results=True, **kwargs)
