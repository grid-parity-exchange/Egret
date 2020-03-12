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
import pao.bilevel as bi
import numpy as np
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen
import egret.model_library.extensions.utils_bilevel_nk as utils
import egret.model_library.extensions.master_bilevel_nk as cons
import egret.model_library.extensions.subproblem_bilevel_nk as subcons
import egret.model_library.decl as decl

from egret.model_library.defn import CoordinateType, ApproximationType, BasePointType
from egret.data.data_utils import map_items, zip_items
from math import pi, radians


def create_bigm(model, md):
    branches = dict(md.elements(element_type='branch'))

    model.BIGM = dict()
    for branch_name, branch in branches.items():
        branch = branches[branch_name]

        tau = 1.0
        shift = 0.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = radians(branch['transformer_phase_shift'])

        x = branch['reactance']
        b = -1/(tau*x)

        model.BIGM[branch_name] = b * (2*pi + shift) + 1.


def create_master(model_data, k=1):
    md = model_data.clone_in_service()
    tx_utils.scale_ModelData_to_pu(md, inplace = True)

    relays = dict(md.elements(element_type='relay'))
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    loads = dict(md.elements(element_type='load'))
    branches = dict(md.elements(element_type='branch'))

    relay_attrs = md.attributes(element_type='relay')
    gen_attrs = md.attributes(element_type='generator')
    load_attrs = md.attributes(element_type='load')
    branch_attrs = md.attributes(element_type='branch')


    model = pe.ConcreteModel()

    ### declare (and fix) the loads at the buses
    bus_p_loads, _ = tx_utils.dict_of_bus_loads(buses, loads)
    buses_with_loads = list(k for k in bus_p_loads.keys() if bus_p_loads[k] != 0.)

    relay_branches = utils.dict_of_relay_branches(relays, branches)
    branch_relays = utils.dict_of_branch_relays(relays, branches)
    relay_branch_tuple = utils.relay_branch_tuple(relay_branches)
    relay_gens = utils.dict_of_relay_gens(relays, gens)
    gen_relays = utils.dict_of_gen_relays(relays, gens)
    relay_gen_tuple = utils.relay_branch_tuple(relay_gens)
    relay_loads = utils.dict_of_relay_loads(relays, loads, buses_with_loads)
    load_relays = utils.dict_of_load_relays(relays, buses_with_loads)
    relay_load_tuple = utils.relay_branch_tuple(relay_loads)

    decl.declare_var('load_shed', model, buses_with_loads, initialize=0.0, domain=pe.NonNegativeReals)
    decl.declare_var('delta', model, relay_attrs['names'], domain=pe.Binary) # relays compromised
    decl.declare_var('u', model, buses_with_loads, domain=pe.Binary) # load available
    decl.declare_var('v', model, gen_attrs['names'], domain=pe.Binary) # generator available
    decl.declare_var('w', model, branch_attrs['names'], domain=pe.Binary) # line available

    cons.declare_budget(model, k, relays)
    cons.declare_load_compromised(model, relay_load_tuple)
    cons.declare_load_uncompromised(model, buses_with_loads, load_relays)
    cons.declare_branch_compromised(model, relay_branch_tuple)
    cons.declare_branch_uncompromised(model, branch_attrs['names'], branch_relays)
    cons.declare_gen_compromised(model, relay_gen_tuple)
    cons.declare_gen_uncompromised(model, gen_attrs['names'], gen_relays)

    model.obj = pe.Objective(expr=sum(model.load_shed[l] for l in buses_with_loads), sense=pe.maximize)

    return model, md


def create_gdp_subproblem(model, model_data, include_angle_diff_limits=False):
    md = model_data
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

    model.subproblem = bi.SubModel(fixed=(model.u, model.v, model.w))

    ### declare (and fix) the loads at the buses
    bus_p_loads, _ = tx_utils.dict_of_bus_loads(buses, loads)
    buses_with_loads = list(k for k in bus_p_loads.keys() if bus_p_loads[k] != 0.)

    libbus.declare_var_pl(model.subproblem, bus_attrs['names'], initialize=bus_p_loads)
    model.subproblem.pl.fix()

    ### declare the fixed shunts at the buses
    _, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    ### declare the polar voltages
    va_bounds = {k: (-pi, pi) for k in bus_attrs['va']}
    libbus.declare_var_va(model.subproblem, bus_attrs['names'], initialize=bus_attrs['va'],
                          bounds=va_bounds
                          )

    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    ref_angle = md.data['system']['reference_bus_angle']
    model.subproblem.va[ref_bus].fix(radians(ref_angle))

    ### declare the generator real power
    pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model.subproblem, gen_attrs['names'], initialize=pg_init,
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

    libbranch.declare_var_pf(model=model.subproblem,
                             index_set=branch_attrs['names'],
                             initialize=pf_init,
                             bounds=pf_bounds
                             )

    # need to include variable references on subproblem to variables, which exist on the master block
    bi.components.varref(model.subproblem)

    ### declare the branch power flow disjuncts (LHS is status quo, RHS is compromised)
    libbranch.declare_eq_branch_power_btheta_approx(model=model.subproblem,
                                                    index_set=branch_attrs['names'],
                                                    branches=branches
                                                    )
    subcons.declare_eq_branch_power_off(model=model.subproblem,
                                        index_set=branch_attrs['names'],
                                        branches=branches
                                        )
    subcons.disjunctify(model=model.subproblem,
                        indicator_name='pf_branch_indicator',
                        disjunct_name='pf_branch_disjunct',
                        LHS_disjunct_set=model.subproblem.eq_pf_branch,
                        RHS_disjunct_set=model.subproblem.eq_pf_branch_off
                        )

    ### declare the load shed disjuncts (LHS is status quo, RHS is compromised)
    subcons.declare_ineq_load_shed_ub(model=model.subproblem,
                                      index_set=buses_with_loads)
    subcons.declare_ineq_load_shed_lb(model=model.subproblem,
                                      index_set=buses_with_loads)
    subcons.declare_ineq_load_shed_lb_off(model=model.subproblem,
                                      index_set=buses_with_loads)
    subcons.disjunctify(model=model.subproblem,
                        indicator_name='load_shed_indicator',
                        disjunct_name='load_shed_disjunct',
                        LHS_disjunct_set=model.subproblem.ineq_load_shed_lb,
                        RHS_disjunct_set=model.subproblem.ineq_load_shed_lb_off
                        )

    ### declare the generator disjuncts (LHS is status quo, RHS is compromised)
    subcons.declare_ineq_gen_on(model=model.subproblem,
                             index_set=gen_attrs['names'],
                             gens=gens)
    subcons.declare_ineq_gen_off(model=model.subproblem,
                                 index_set=gen_attrs['names'],
                                 gens=gens)
    subcons.disjunctify(model=model.subproblem,
                        indicator_name='gen_indicator',
                        disjunct_name='gen_disjunct',
                        LHS_disjunct_set=model.subproblem.ineq_gen,
                        RHS_disjunct_set=model.subproblem.ineq_gen_off
                        )

    ### declare the p balance
    rhs_kwargs = {'include_feasibility_slack_neg':'load_shed'}
    libbus.declare_eq_p_balance_dc_approx(model=model.subproblem,
                                          index_set=bus_attrs['names'],
                                          bus_p_loads=bus_p_loads,
                                          gens_by_bus=gens_by_bus,
                                          bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                          inlet_branches_by_bus=inlet_branches_by_bus,
                                          outlet_branches_by_bus=outlet_branches_by_bus,
                                          approximation_type=ApproximationType.BTHETA,
                                          **rhs_kwargs
                                          )

    ### declare the real power flow limits
    libbranch.declare_ineq_p_branch_thermal_lbub(model=model.subproblem,
                                                 index_set=branch_attrs['names'],
                                                 branches=branches,
                                                 p_thermal_limits=p_max,
                                                 approximation_type=ApproximationType.BTHETA
                                                 )

    ### declare angle difference limits on interconnected buses
    if include_angle_diff_limits:
        libbranch.declare_ineq_angle_diff_branch_lbub(model=model.subproblem,
                                                      index_set=branch_attrs['names'],
                                                      branches=branches,
                                                      coordinate_type=CoordinateType.POLAR
                                                      )

    model.subproblem.obj = pe.Objective(expr=sum(model.load_shed[l] for l in buses_with_loads), sense=pe.minimize)

    return model, md


def create_explicit_subproblem(model, model_data, include_angle_diff_limits=False, include_bigm=False):
    md = model_data
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

    model.subproblem = bi.SubModel(fixed=(model.u, model.v, model.w))

    ### declare (and fix) the loads at the buses
    bus_p_loads, _ = tx_utils.dict_of_bus_loads(buses, loads)
    buses_with_loads = list(k for k in bus_p_loads.keys() if bus_p_loads[k] != 0.)

    libbus.declare_var_pl(model.subproblem, bus_attrs['names'], initialize=bus_p_loads)
    model.subproblem.pl.fix()

    ### declare the fixed shunts at the buses
    _, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    ### declare the polar voltages
    va_bounds = {k: (-pi, pi) for k in bus_attrs['va']}
    libbus.declare_var_va(model.subproblem, bus_attrs['names'], initialize=bus_attrs['va'],
                          bounds=va_bounds
                          )

    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    ref_angle = md.data['system']['reference_bus_angle']
    model.subproblem.va[ref_bus].fix(radians(ref_angle))

    ### declare the generator real power
    pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model.subproblem, gen_attrs['names'], initialize=pg_init,
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

    libbranch.declare_var_pf(model=model.subproblem,
                             index_set=branch_attrs['names'],
                             initialize=pf_init,
                             bounds=pf_bounds
                             )

    # need to include variable references on subproblem to variables, which exist on the master block
    bi.components.varref(model.subproblem)

    if include_bigm:
        # create big-M
        create_bigm(model.subproblem, md)
        ### declare the branch power flow disjuncts
        subcons.declare_eq_branch_power_btheta_approx_bigM(model=model.subproblem,
                                                           index_set=branch_attrs['names'],
                                                           branches=branches
                                                           )
    else:
        subcons.declare_eq_branch_power_btheta_approx_nonlin(model=model.subproblem,
                                                           index_set=branch_attrs['names'],
                                                           branches=branches
                                                           )


    ### declare the load shed disjuncts
    subcons.declare_ineq_load_shed(model=model.subproblem,
                                      index_set=buses_with_loads)

    ### declare the generator disjuncts (LHS is status quo, RHS is compromised)
    subcons.declare_ineq_gen(model=model.subproblem,
                             index_set=gen_attrs['names'],
                             gens=gens)

    ### declare the p balance
    rhs_kwargs = {'include_feasibility_slack_neg':'load_shed'}
    libbus.declare_eq_p_balance_dc_approx(model=model.subproblem,
                                          index_set=bus_attrs['names'],
                                          bus_p_loads=bus_p_loads,
                                          gens_by_bus=gens_by_bus,
                                          bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                          inlet_branches_by_bus=inlet_branches_by_bus,
                                          outlet_branches_by_bus=outlet_branches_by_bus,
                                          approximation_type=ApproximationType.BTHETA,
                                          **rhs_kwargs
                                          )

    ### declare the real power flow limits
    libbranch.declare_ineq_p_branch_thermal_lbub(model=model.subproblem,
                                                 index_set=branch_attrs['names'],
                                                 branches=branches,
                                                 p_thermal_limits=p_max,
                                                 approximation_type=ApproximationType.BTHETA
                                                 )

    ### declare angle difference limits on interconnected buses
    if include_angle_diff_limits:
        libbranch.declare_ineq_angle_diff_branch_lbub(model=model.subproblem,
                                                      index_set=branch_attrs['names'],
                                                      branches=branches,
                                                      coordinate_type=CoordinateType.POLAR
                                                      )

    model.subproblem.obj = pe.Objective(expr=sum(model.load_shed[l] for l in buses_with_loads), sense=pe.minimize)

    return model, md


def solve_bilevel_nk(model_data,
                solver,
                timelimit = None,
                solver_tee = True,
                symbolic_solver_labels = False,
                options = None,
                return_model = False,
                return_results = False,
                **kwargs):
    '''
    Create and solve a new worst-case attacker defender

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

    generate_explicit_subproblem = kwargs.pop('explicit_subproblem',True)

    m, md = create_master(model_data)
    if generate_explicit_subproblem:
        m, md = create_explicit_subproblem(m, model_data)
    else:
        m, md = create_gdp_subproblem(m, model_data)
        pe.TransformationFactory('gdp.bigm').apply_to(m)

    m.dual = pe.Suffix(direction=pe.Suffix.IMPORT)

    opt = pe.SolverFactory('pao.bilevel.ld', solver=solver)#, symbolic_solver_labels=symbolic_solver_labels)
    results = opt.solve(m, tee=solver_tee)

    import pdb
    pdb.set_trace()

    # save results data to ModelData object
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))

    md.data['system']['total_cost'] = value(m.obj)

    for g,g_dict in gens.items():
        g_dict['pg'] = value(m.pg[g])

    for k, k_dict in branches.items():
        k_dict['pf'] = value(m.pf[k])

    for b,b_dict in buses.items():
        b_dict['pl'] = value(m.pl[b])
        b_dict['lmp'] = value(m.dual[m.eq_p_balance[b]])
        b_dict['va'] = value(m.va[b])

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

    path = os.path.dirname(__file__)
    print(path)
    filename = 'pglib_opf_case5_pjm.m'
    test_case = os.path.join(path, '../../download/pglib-opf-master/', filename)
    md_dict = create_ModelData(test_case)

    relays = dict()

    buses = dict(md_dict.elements(element_type='bus'))
    for bus_name, bus in buses.items():
        relay_name = bus_name
        relays[relay_name] = dict()
        relays[relay_name]['branch'] = list()
        relays[relay_name]['gen'] = list()
        relays[relay_name]['load'] = list()

    branches = dict(md_dict.elements(element_type='branch'))
    for branch_name, branch in branches.items():
        relay_name = branch['from_bus']
        relay = relays[relay_name]
        relay['branch'].append(branch_name)
        branch['relay'] = list()
        branch['relay'].append(relay_name)

        relay_name = branch['to_bus']
        relay = relays[relay_name]
        relay['branch'].append(branch_name)
        branch['relay'].append(relay_name)

    gens = dict(md_dict.elements(element_type='generator'))
    for gen_name, gen in gens.items():
        relay_name = gen['bus']
        relay = relays[relay_name]
        relay['gen'].append(gen_name)
        gen['relay'] = list()
        gen['relay'].append(relay_name)

    loads = dict(md_dict.elements(element_type='load'))
    for load_name, load in loads.items():
        relay_name = load['bus']
        relay = relays[relay_name]
        relay['load'].append(load['bus'])
        load['relay'] = list()
        load['relay'].append(relay_name)

    md_dict.data['elements']['relay'] = relays

    kwargs = {'explicit_subproblem': True}
    md_serialization, results = solve_bilevel_nk(md_dict, "gurobi", solver_tee=False,
                                            return_results=True, **kwargs)

