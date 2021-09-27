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
import egret.models.lpac as lpac
from egret.model_library.defn import CoordinateType, ApproximationType
from math import pi, radians
import pdb


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


def create_master(model_data, k=4, attack_blacklist = []):
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
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    loads = dict(md.elements(element_type='load'))
    branches = dict(md.elements(element_type='branch'))

    ### create dictionaries across object attributes for an object of the same set type
    bus_attrs = md.attributes(element_type='bus')
    gen_attrs = md.attributes(element_type='generator')
    branch_attrs = md.attributes(element_type='branch')

    ### create grid relationship dictionaries
    bus_by_gen = {gen: gen_dict['bus'] for gen, gen_dict in gens.items()}
    buses_by_branch = {branch: [branch_dict['from_bus'], branch_dict['to_bus']] for branch, branch_dict in branches.items()}

    ### declare new Pyomo model
    model = pe.ConcreteModel()

    ### declare (and fix) the loads at the buses
    bus_p_loads, _ = tx_utils.dict_of_bus_loads(buses, loads)
    buses_with_loads = list(k for k in bus_p_loads.keys() if bus_p_loads[k] != 0.)
    load_shed_bounds  = {load: (0, bus_p_loads[load]) for load in bus_p_loads}
    
    ### upper-level (attacker) variables
    decl.declare_var('load_shed', model, buses_with_loads, initialize=0.0, domain=pe.NonNegativeReals, bounds = load_shed_bounds)
    decl.declare_var('delta_bus', model, bus_attrs['names'], domain=pe.Binary) # bus attacked
    decl.declare_var('delta_load', model, buses_with_loads, domain=pe.Binary) # load attacked
    decl.declare_var('delta_gen', model, gen_attrs['names'], domain=pe.Binary) # generator attacked
    decl.declare_var('delta_branch', model, branch_attrs['names'], domain=pe.Binary) # line attacked
    decl.declare_var('u', model, bus_attrs['names'], domain=pe.Binary, initialize = 1) # load available
    decl.declare_var('v', model, gen_attrs['names'], domain=pe.Binary, initialize = 1) # generator available
    decl.declare_var('w', model, branch_attrs['names'], domain=pe.Binary, initialize = 1) # line available

    ### upper-level constraints
    cons.declare_physical_budget(model, k)
    cons.declare_physical_load_compromised(model, buses_with_loads)
    cons.declare_physical_load_uncompromised(model, buses_with_loads)
    cons.declare_physical_line_compromised(model, branch_attrs['names'], buses_by_branch)
    cons.declare_physical_line_uncompromised(model, branch_attrs['names'], buses_by_branch)
    cons.declare_physical_gen_compromised(model, gen_attrs['names'], bus_by_gen)
    cons.declare_physical_gen_uncompromised(model, gen_attrs['names'], bus_by_gen)

    ### add cuts
    #model.cut_list = pe.ConstraintList()
    #for attack in attack_blacklist:
    #    expr = sum(1 - model.delta_bus[bus] for bus in attack['bus'])
    #    expr += sum(1 - model.delta_branch[branch] for branch in attack['branch'])
    #    expr += sum(1 - model.delta_gen[gen] for gen in attack['gen'])
    #    expr += sum(1 - model.delta_load[load] for load in attack['load'])
    #    model.cut_list.add(expr = expr >= 1)

    ### upper-level objective for interdiction problem (opposite to lower-level objective)
    model.obj = pe.Objective(expr=sum(model.load_shed[l] for l in buses_with_loads), sense=pe.maximize)

    return model, md


def create_explicit_subproblem(model, model_data, include_angle_diff_limits=False, include_bigm=False):

    ### declare lower-level as a PAO (Pyomo-extension) submodel;
    ### be explicit in specifying upper-level variables that appear in this model
    model.subproblem = bi.SubModel(fixed=(model.u, model.v, model.w))

    lpac.create_cold_start_lpac_model(model_data, cosine_segment_count = 20, lower_bound = -pi/3, upper_bound = pi/3, include_feasibility_slack = True, mode="curvature", parent_model = model)
    return model, model_data


def solve_bilevel_physical_nk(model_data,
                solver,
                solver_tee = True,
                return_model = False,
                return_results = False,
                return_power_flow = False,
                allow_gen = True,
                allow_bus = True,
                allow_load = True,
                allow_branch = True,
                attack_blacklist = [],
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
    m, md = create_master(md,attack_budget_k, attack_blacklist = attack_blacklist)
    ### create lower-level of the bilevel problem
    m, md = create_explicit_subproblem(m, md,include_bigm=False)

    #### disable components that should not be attacked
    if not allow_load:
        m.delta_load.fix(0)
    if not allow_gen:
        m.delta_gen.fix(0)
    if not allow_branch:
        m.delta_branch.fix(0)
    if not allow_bus:
        m.delta_bus.fix(0)

    ### use PAO (Pyomo-extension) to do the following:
    ### 1. Transform the lower-level primal problem into it's corresponding dual problem
    ### 2. Apply Pyomo.GDP transformations to handle bilinear terms (Big-M)
    ### 3. Solve formulation (upper-level primal with lower-level dual) as a single level MILP
    ### 4. Take optimal solution from MILP, fix upper-level variables that appear in the
    ### lower-level problem, and resolve to dectermine primal variable solution for the lower-level
    opt = pe.SolverFactory('pao.bilevel.ld', solver=solver)
    ## need to fine-tune bigM and mipgap -- make sure that both the solve and resolve result in the same
    ## best objective
    opt.options.setdefault('bigM', 10000)
    opt.options.setdefault('mipgap', 0.001)
    opt.options.setdefault('threads', 2)
    results = opt.solve(m, tee=solver_tee)
    
    ### save results data to ModelData object
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))

    md.data['system']['total_cost'] = value(m.obj)

    power_flow = m.subproblem
    for g,g_dict in gens.items():
        g_dict['pg'] = value(power_flow.pg[g])

    for k, k_dict in branches.items():
        k_dict['pf'] = value(power_flow.pf[k])

    for b,b_dict in buses.items():
        b_dict['pl'] = value(power_flow.pl[b])
        b_dict['va'] = value(power_flow.va[b])

    unscale_ModelData_to_pu(md, inplace=True)

    ### return model_data (md), model (m), and/or results (results) objects
    return_item = [md]
    if return_model:
        return_item.append(m)
    if return_power_flow:
        return_item.append(power_flow)
    if return_results:
        return_item.append(results)
    return tuple(return_item)

    '''if return_model and return_results:
        return md, m, results
    elif return_model:
        return md, m
    elif return_results:
        return md, results
    return md'''


if __name__ == '__main__':
    import os
    from egret.parsers.matpower_parser import create_ModelData

    ### Example code on loading data
    path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    ### MATPOWER format *.m file for numerous power systems available at github for pglib-opf
    filename = 'case24_ieee_rts_example.m'
    test_case = os.path.join(path, 'download', filename)
    md_dict = create_ModelData(test_case)
    attack_budget = 10

    ### solve the bilevel interdiction problem
    md_serialization, model, results = solve_bilevel_physical_nk(md_dict, "gurobi", solver_tee=True,
                                            allow_bus = False, return_model = True, return_results=True, attack_budget_k = attack_budget, symbolic_solver_labels = True)

    pdb.set_trace()