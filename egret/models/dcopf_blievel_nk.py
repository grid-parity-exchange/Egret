#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module provides functions that create the modules for bilevel interdiction
based on Salmeron (2004, 2009) work.
"""

import pyomo.environ as pe
import pao.bilevel as bi
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen
import egret.model_library.extensions.utils_bilevel_nk as utils
import egret.model_library.extensions.master_bilevel_nk as cons
import egret.model_library.extensions.subproblem_bilevel_nk as subcons
import egret.model_library.decl as decl
from egret.model_library.defn import CoordinateType, ApproximationType
from math import pi, radians


def _create_bigm(model, md):
    """
    Create the big-M transformation explicitly
    """
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

        model.BIGM[branch_name] = b * (2*pi + shift) + 1000.


def create_master(model_data, k=1):
    """
    Create the upper-level (master) of the bilevel problem

    Arguments:
        model_data: An Egret dict of dict that stores data for the power system
        k: A positive integer indicating the number of relays that can be attacked

    Returns: Tuple with the following values:
        model: A Pyomo model representing the algebraic form of the bilevel problem
        md: The model_data object associated to the model
    """
    ### power system data
    md = model_data

    ### create dictionaries of object sets
    relays = dict(md.elements(element_type='relay'))
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    loads = dict(md.elements(element_type='load'))
    branches = dict(md.elements(element_type='branch'))

    ### create dictionaries across object attributes for an object of the same set type
    relay_attrs = md.attributes(element_type='relay')
    gen_attrs = md.attributes(element_type='generator')
    branch_attrs = md.attributes(element_type='branch')

    ### declare new Pyomo model
    model = pe.ConcreteModel()

    ### declare (and fix) the loads at the buses
    bus_p_loads, _ = tx_utils.dict_of_bus_loads(buses, loads)
    buses_with_loads = list(k for k in bus_p_loads.keys() if bus_p_loads[k] != 0.)

    ### relay-to-power-device mappings
    relay_branches = utils.dict_of_relay_branches(relays, branches)
    branch_relays = utils.dict_of_branch_relays(relays, branches)
    relay_branch_tuple = utils.relay_branch_tuple(relay_branches)
    relay_gens = utils.dict_of_relay_gens(relays, gens)
    gen_relays = utils.dict_of_gen_relays(relays, gens)
    relay_gen_tuple = utils.relay_branch_tuple(relay_gens)
    relay_loads = utils.dict_of_relay_loads(relays, loads, buses_with_loads)
    load_relays = utils.dict_of_load_relays(relays, buses_with_loads)
    relay_load_tuple = utils.relay_branch_tuple(relay_loads)

    ### upper-level (attacker) variables
    decl.declare_var('load_shed', model, buses_with_loads, initialize=0.0, domain=pe.NonNegativeReals)
    decl.declare_var('delta', model, relay_attrs['names'], domain=pe.Binary) # relays compromised
    decl.declare_var('u', model, buses_with_loads, domain=pe.Binary) # load available
    decl.declare_var('v', model, gen_attrs['names'], domain=pe.Binary) # generator available
    decl.declare_var('w', model, branch_attrs['names'], domain=pe.Binary) # line available

    ### upper-level constraints
    cons.declare_budget(model, k, relays) # note that all k are costed equally in current implementation
    cons.declare_load_compromised(model, relay_load_tuple)
    cons.declare_load_uncompromised(model, buses_with_loads, load_relays)
    cons.declare_branch_compromised(model, relay_branch_tuple)
    cons.declare_branch_uncompromised(model, branch_attrs['names'], branch_relays)
    cons.declare_gen_compromised(model, relay_gen_tuple)
    cons.declare_gen_uncompromised(model, gen_attrs['names'], gen_relays)

    ### upper-level objective for interdiction problem (opposite to lower-level objective)
    model.obj = pe.Objective(expr=sum(model.load_shed[l] for l in buses_with_loads), sense=pe.maximize)

    return model, md


def create_explicit_subproblem(model, model_data, include_angle_diff_limits=False, include_bigm=False):
    ### power system data
    md = model_data

    ### create dictionaries of object sets
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    loads = dict(md.elements(element_type='load'))
    shunts = dict(md.elements(element_type='shunt'))

    ### create dictionaries across object attributes for an object of the same set type
    gen_attrs = md.attributes(element_type='generator')
    bus_attrs = md.attributes(element_type='bus')
    branch_attrs = md.attributes(element_type='branch')

    inlet_branches_by_bus, outlet_branches_by_bus = \
        tx_utils.inlet_outlet_branches_by_bus(branches, buses)
    gens_by_bus = tx_utils.gens_by_bus(buses, gens)

    ### declare lower-level as a PAO (Pyomo-extension) submodel;
    ### be explicit in specifying upper-level variables that appear in this model
    model.subproblem = bi.SubModel(fixed=(model.u, model.v, model.w))

    ### declare (and fix) the loads at the buses
    bus_p_loads, _ = tx_utils.dict_of_bus_loads(buses, loads)
    buses_with_loads = list(k for k in bus_p_loads.keys() if bus_p_loads[k] != 0.)
    #libbus.declare_var_pl(model.subproblem, bus_attrs['names'], initialize=bus_p_loads)
    #model.subproblem.pl.fix()
    model.subproblem.pl = bus_p_loads

    ### declare the fixed shunts at the buses
    _, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    ### declare the polar voltages
    va_bounds = {k: (-pi, pi) for k in bus_attrs['va']}
    va_init = {k: bus_attrs['va'][k]*(pi/180) for k in bus_attrs['va']}
    libbus.declare_var_va(model.subproblem, bus_attrs['names'], initialize=bus_attrs['va'],
                          bounds=va_bounds
                          )

    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    ref_angle = md.data['system']['reference_bus_angle']
    model.subproblem.va[ref_bus].fix(radians(ref_angle))

    ### declare the generator real power
    pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model.subproblem, gen_attrs['names'], initialize=pg_init)

    ### declare the current flows in the branches
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    vj_init = {k: bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]) for k in bus_attrs['vm']}
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
                             initialize=pf_init)

    # need to include variable references on subproblem to variables, which exist on the master block
    bi.components.varref(model.subproblem)

    if include_bigm:
        # create big-M
        _create_bigm(model.subproblem, md)
        ### declare the branch power flow disjuncts
        subcons.declare_eq_branch_power_btheta_approx_bigM(model=model.subproblem,
                                                           index_set=branch_attrs['names'],
                                                           branches=branches
                                                           )

        ### declare the real power flow limits
        p_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
        subcons.declare_ineq_p_branch_thermal_lbub_switch(model=model.subproblem,
                                                          index_set=branch_attrs['names'],
                                                          p_thermal_limits=p_max
                                                          )

    else:
        ### declare the branch power flow with indicator variable in the bilinear term
        subcons.declare_eq_branch_power_btheta_approx_nonlin(model=model.subproblem,
                                                             index_set=branch_attrs['names'],
                                                             branches=branches
                                                             )

        ### declare the real power flow limits
        p_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
        libbranch.declare_ineq_p_branch_thermal_lbub(model=model.subproblem,
                                                     index_set=branch_attrs['names'],
                                                     branches=branches,
                                                     p_thermal_limits=p_max,
                                                     approximation_type=ApproximationType.BTHETA
                                                     )


    ### declare the load shed
    subcons.declare_ineq_load_shed(model=model.subproblem,
                                      index_set=buses_with_loads)

    ### declare the generator compromised
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

    ### declare angle difference limits on interconnected buses
    if include_angle_diff_limits:
        libbranch.declare_ineq_angle_diff_branch_lbub(model=model.subproblem,
                                                      index_set=branch_attrs['names'],
                                                      branches=branches,
                                                      coordinate_type=CoordinateType.POLAR
                                                      )

    ### lower-level objective for interdiction problem (opposite to upper-level objective)
    model.subproblem.obj = pe.Objective(expr=sum(model.load_shed[l] for l in buses_with_loads), sense=pe.minimize)

    return model, md


def solve_bilevel_nk(model_data,
                solver,
                solver_tee = True,
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
    solver_tee : bool (optional)
        Display solver log. Default is True.
    return_model : bool (optional)
        If True, returns the pyomo model object
    return_results : bool (optional)
        If True, returns the pyomo results object
    kwargs : dictionary (optional)
        Additional arguments for building model
    '''

    import pyomo.environ as pe
    from pyomo.environ import value
    from egret.model_library.transmission.tx_utils import \
        scale_ModelData_to_pu, unscale_ModelData_to_pu

    md = model_data.clone_in_service()
    scale_ModelData_to_pu(md, inplace = True)

    ### pop from kwargs the number k for N-k contingency of relay IPs
    attack_budget_k = kwargs.pop('attack_budget_k',1)

    ### create upper-level of the bilevel problem
    m, md = create_master(md,attack_budget_k)
    ### create lower-level of the bilevel problem
    m, md = create_explicit_subproblem(m, md,include_bigm=False)

    ### use PAO (Pyomo-extension) to do the following:
    ### 1. Transform the lower-level primal problem into it's corresponding dual problem
    ### 2. Apply Pyomo.GDP transformations to handle bilinear terms (Big-M)
    ### 3. Solve formulation (upper-level primal with lower-level dual) as a single level MILP
    ### 4. Take optimal solution from MILP, fix upper-level variables that appear in the
    ### lower-level problem, and resolve to determine primal variable solution for the lower-level
    opt = pe.SolverFactory('pao.bilevel.ld', solver=solver)
    ## need to fine-tune bigM and mipgap -- make sure that both the solve and resolve result in the same
    ## best objective
    opt.options.setdefault('bigM', 100)
    opt.options.setdefault('mipgap', 0.001)
    results = opt.solve(m, tee=solver_tee)

    ### save results data to ModelData object
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))

    md.data['system']['total_cost'] = value(m.obj)

    m = m.subproblem
    for g,g_dict in gens.items():
        g_dict['pg'] = value(m.pg[g])

    for k, k_dict in branches.items():
        k_dict['pf'] = value(m.pf[k])

    for b,b_dict in buses.items():
        b_dict['pl'] = value(m.pl[b])
        b_dict['va'] = value(m.va[b])

    unscale_ModelData_to_pu(md, inplace=True)

    ### return model_data (md), model (m), and/or results (results) objects
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

    ### Example code on loading data
    path = os.path.dirname(__file__)
    print(path)
    ### MATPOWER format *.m file for numerous power systems available at github for pglib-opf
    filename = 'case24_ieee_rts_example.m'
    test_case = os.path.join(path, '../../download/', filename)
    md_dict = create_ModelData(test_case)


    ### Generate relay-to-power-device mappings for the synthetic dataset
    relays = dict()

    buses = dict(md_dict.elements(element_type='bus'))
    for bus_name, bus in buses.items():
        relay_name = 'r' + bus_name
        relays[relay_name] = dict()
        relays[relay_name]['branch'] = list()
        relays[relay_name]['gen'] = list()
        relays[relay_name]['load'] = list()

    branches = dict(md_dict.elements(element_type='branch'))
    for branch_name, branch in branches.items():
        relay_name = 'r' + branch['from_bus']
        relay = relays[relay_name]
        relay['branch'].append(branch_name)
        branch['relay'] = list()
        branch['relay'].append(relay_name)

        relay_name = 'r' + branch['to_bus']
        relay = relays[relay_name]
        relay['branch'].append(branch_name)
        branch['relay'].append(relay_name)

    gens = dict(md_dict.elements(element_type='generator'))
    for gen_name, gen in gens.items():
        relay_name = 'r' + gen['bus']
        relay = relays[relay_name]
        relay['gen'].append(gen_name)
        gen['relay'] = list()
        gen['relay'].append(relay_name)

    loads = dict(md_dict.elements(element_type='load'))
    for load_name, load in loads.items():
        relay_name = 'r' + load['bus']
        relay = relays[relay_name]
        relay['load'].append(load['bus'])
        load['relay'] = list()
        load['relay'].append(relay_name)

    md_dict.data['elements']['relay'] = relays

    ### specify the attack budget
    kwargs = {'attack_budget_k': 1}

    ### solve the bilevel interdiction problem
    md_serialization, results = solve_bilevel_nk(md_dict, "gurobi", solver_tee=True,
                                            return_results=True, **kwargs)

