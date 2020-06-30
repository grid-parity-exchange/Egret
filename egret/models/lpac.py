"""
This module provides functions for formulations of the LPAC model of Coffrin and Van Hentenryck (2013)

"""
import pyomo.environ as pe
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen

import egret.model_library.decl as decl
from egret.model_library.defn import FlowType, CoordinateType, ApproximationType

from egret.data.data_utils import map_items, zip_items
from math import pi, radians

from egret.models.acopf import _include_feasibility_slack



def create_hot_start_lpac_model(model_data, voltages, cosine_segment_count, include_feasibility_slack = False):
	"""
	The hot start LPAC model assumes that voltages are known, e.g. from an AC base point solution.
	"""

def create_warm_start_lpac_model(model_data, cosine_segment_count, include_feasibility_slack = False):
	"""
	The warm start LPAC model assumes that target voltages can be given, but not all voltages are known. 
	"""

def create_cold_start_lpac_model(model_data, cosine_segment_count, include_feasibility_slack = False):
	"""
	The cold start LPAC model assumes that no target voltages are available and that all voltages are initially approximated as 1 pu. 
	"""
	###Grid data
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

    ### declare the polar voltages
    libbus.declare_var_va(model, bus_attrs['names'], initialize=bus_attrs['va']
                        )

    ### declare the voltage change variables
    decl.declare_var(phi, model, bus_attrs['names'])

    ### declare the cosine approximation variables
    cos_hat_bounds = {k: (0, 1) for k in branch_attrs['names']}
    decl.declare_var(cos_hat, model, branch_attrs['names'], bounds = cos_hat_bounds)

    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    #ref_angle = md.data['system']['reference_bus_angle']
    model.va[ref_bus].fix(radians(0.0))

    ### declare (and fix) the loads at the buses
    bus_p_loads, bus_q_loads = tx_utils.dict_of_bus_loads(buses, loads)

    libbus.declare_var_pl(model, bus_attrs['names'], initialize=bus_p_loads)
    libbus.declare_var_ql(model, bus_attrs['names'], initialize=bus_q_loads)
    model.pl.fix()
    model.ql.fix()

    ### declare the generator real and reactive power
    pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model, gen_attrs['names'], initialize=pg_init,
                          bounds=zip_items(gen_attrs['p_min'], gen_attrs['p_max'])
                          )

    qg_init = {k: (gen_attrs['q_min'][k] + gen_attrs['q_max'][k]) / 2.0 for k in gen_attrs['qg']}
    libgen.declare_var_qg(model, gen_attrs['names'], initialize=qg_init,
                          bounds=zip_items(gen_attrs['q_min'], gen_attrs['q_max'])
                          )

    ### declare the current flows in the branches
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    vj_init = {k: bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    s_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
    s_lbub = dict()
    for k in branches.keys():
        if s_max[k] is None:
            s_lbub[k] = (None, None)
        else:
            s_lbub[k] = (-s_max[k],s_max[k])
    pf_bounds = s_lbub
    #pt_bounds = s_lbub
    qf_bounds = s_lbub
    #qt_bounds = s_lbub
    pf_init = dict()
    #pt_init = dict()
    qf_init = dict()
   # qt_init = dict()
    for branch_name, branch in branches.items():
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']
        y_matrix = tx_calc.calculate_y_matrix_from_branch(branch)
        ifr_init = tx_calc.calculate_ifr(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        ifj_init = tx_calc.calculate_ifj(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        itr_init = tx_calc.calculate_itr(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        itj_init = tx_calc.calculate_itj(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        pf_init[branch_name] = tx_calc.calculate_p(ifr_init, ifj_init, vr_init[from_bus], vj_init[from_bus])
        #pt_init[branch_name] = tx_calc.calculate_p(itr_init, itj_init, vr_init[to_bus], vj_init[to_bus])
        qf_init[branch_name] = tx_calc.calculate_q(ifr_init, ifj_init, vr_init[from_bus], vj_init[from_bus])
        #qt_init[branch_name] = tx_calc.calculate_q(itr_init, itj_init, vr_init[to_bus], vj_init[to_bus])

    libbranch.declare_var_pf(model=model,
                             index_set=branch_attrs['names'],
                             initialize=pf_init,
                             bounds=pf_bounds
                             )
    # libbranch.declare_var_pt(model=model,
    #                          index_set=branch_attrs['names'],
    #                          initialize=pt_init,
    #                          bounds=pt_bounds
    #                         )
    libbranch.declare_var_qf(model=model,
                             index_set=branch_attrs['names'],
                             initialize=qf_init,
                             bounds=qf_bounds
                             )
    # libbranch.declare_var_qt(model=model,
    #                          index_set=branch_attrs['names'],
    #                          initialize=qt_init,
    #                          bounds=qt_bounds
    #                          )


    ### Balance equations at a bus (based on Kirchhoff Current Law)

    #Should be able to just use DC OPF approximation of B-theta type? 

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

    #Need one also for q balance 

    con_set = decl.declare_set('_con_eq_q_balance', model, bus_attrs['names'])

    model.eq_q_balance = pe.Constraint(con_set)

    for bus_name in con_set:
        
        q_expr = -sum([model.qf[branch_name] for branch_name in outlet_branches_by_bus[bus_name]])
        q_expr += sum([model.qf[branch_name] for branch_name in inlet_branches_by_bus[bus_name]])
   

        if bus_q_loads[bus_name] != 0.0: # only applies to fixed loads, otherwise may cause an error
            q_expr -= model.ql[bus_name]

        for gen_name in gens_by_bus[bus_name]:
            q_expr += model.qg[gen_name]

        model.eq_q_balance[bus_name] = \
            q_expr == 0.0


    ### Constraints for power in a branch 

    branch_con_set = decl.declare_set('_con_eq_p_q_lpac_branch_power', model, branch_attrs['names'])

    model.eq_pf_branch_t = pe.Constraint(branch_con_set)
    model.eq_qf_branch_t = pe.Constraint(branch_con_set)

    for branch_name in branch_con_set:
    	branch = branches[branch_name]

        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        g = tx_calc.calculate_conductance(branch)
        b = tx_calc.calculate_susceptance(branch)

        model.eq_pf_branch_t[branch_name] = \
        	model.pf[branch_name] == \
        	g - g * model.cos_hat[branch_name] - b * (model.va[from_bus] - model.va[to_bus])

        model.eq_qf_branch_t[branch_name] = \
        	model.qf[branch_name] == \
        	-b - g*(model.va[from_bus] - model.va[to_bus]) + b*model.cos_hat[branch_name] - b*(model.phi[from_bus] - model.phi[to_bus])

    ### Piecewise linear cosine constraints

    #Cosine segment count is part of function definition. Right now, we give the bounds as (-pi/3, pi/3) for the linearization.
    #Probably worthwhile to make this part of a function definition at some point. 

    ell = -pi/3 #lower bound of domain
    h = pi/3 #upper bound of domain

    increment = (h - ell)/(cosine_segment_count + 1)

    model.cosine_bounds = ConstraintList()

    for branch_name in branch_attrs['names']:
    	branch = branches[branch_name]

    	from_bus = branch['from_bus']
    	to_bus = branch['to_bus']

    	model.cosine_bounds.add(model.cos_hat[branch_name] >= (pe.cos(h) - pe.cos(ell))/(h-ell) * ((model.va[from_bus] - model.va[to_bus]) - ell) + pe.cos(ell))

    	for i in range(1, cosine_segment_count + 1):
    		a = i*(h - ell)/(cosine_segment_count + 1)
    		model.cosine_bounds.add(model.cos_hat[branch_name] <= (-pe.sin(a)*(model.va[from_bus] - model.va[to_bus] - a) + pe.cos(a)))




    ### Objective is to maximize cosine hat variables

    obj_expr = sum(m.cos_hat[branch_name] for branch_name in branch_attrs['names'])

    if include_feasibility_slack:
    	obj_expr += penalty_expr

    model.obj = pe.Objective(expr=obj.expr)

    return model, md


def solve_lpac(model_data,
                solver,
                timelimit = None,
                solver_tee = True,
                symbolic_solver_labels = False,
                options = None,
                lpac_model_generator = create_cold_start_lpac_model,
                return_model = False,
                return_results = False,
                **kwargs):
    '''
    Create and solve a new lpac model

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

    m, md = lpac_model_generator(model_data, **kwargs)

    m.dual = pe.Suffix(direction=pe.Suffix.IMPORT)

    m, results, solver = _solve_model(m,solver,timelimit=timelimit,solver_tee=solver_tee,
                              symbolic_solver_labels=symbolic_solver_labels,solver_options=options, return_solver=True)

    # save results data to ModelData object
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))

    md.data['system']['total_cost'] = value(m.obj)


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
    filename = 'pglib_opf_case14_ieee.m'
    matpower_file = os.path.join(path, '../../download/pglib-opf/', filename)
    model_data = create_ModelData(matpower_file)
    kwargs = {'include_feasibility_slack':False}
    md,m,results = solve_lpac(model_data, "gurobi",lpac_model_generator= ,return_model=True, return_results=True,**kwargs)
    md,m,results = solve_lpac(model_data, "gurobi",lpac_model_generator= ,return_model=True, return_results=True,**kwargs)
    md,m,results = solve_lpac(model_data, "gurobi",lpac_model_generator= ,return_model=True, return_results=True,**kwargs)
