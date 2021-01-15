#Implements the multiple choice model described in Aravena et al. -WL

import pyomo.environ as pe

import json

import numpy as np

from PWL_Approx_Functions import *

import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.decl as decl

import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen


from egret.models.acopf import _include_feasibility_slack
from egret.model_library.defn import FlowType, CoordinateType, ApproximationType
from egret.data.data_utils import map_items, zip_items
import math as math
from collections import OrderedDict


def _create_base_ac_with_pwl_approx_model(model_data, branch_dict, Q, include_feasibility_slack=False):
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

    bus_pairs = zip_items(branch_attrs['from_bus'], branch_attrs['to_bus'])
    unique_bus_pairs = list(OrderedDict((val, None) for idx, val in bus_pairs.items()))

    model = pe.ConcreteModel()

    ### declare (and fix) the loads at the buses
    bus_p_loads, bus_q_loads = tx_utils.dict_of_bus_loads(buses, loads)

    libbus.declare_var_pl(model, bus_attrs['names'], initialize=bus_p_loads)
    libbus.declare_var_ql(model, bus_attrs['names'], initialize=bus_q_loads)
    model.pl.fix()
    model.ql.fix()

    ### declare the fixed shunts at the buses
    bus_bs_fixed_shunts, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    libbus.declare_var_vm(model=model, index_set=bus_attrs['names'], initialize=bus_attrs['vm'], bounds=zip_items(bus_attrs['v_min'], bus_attrs['v_max']))

    libbus.declare_var_vmsq(model=model,
                            index_set=bus_attrs['names'],
                            initialize={k: v**2 for k, v in bus_attrs['vm'].items()},
                            bounds=zip_items({k: v**2 for k, v in bus_attrs['v_min'].items()},
                                             {k: v**2 for k, v in bus_attrs['v_max'].items()}))
    # libbranch.declare_var_c(model=model, index_set=unique_bus_pairs)
    # libbranch.declare_var_s(model=model, index_set=unique_bus_pairs)

    ### declare the polar voltages
    va_bounds = {k: (-math.pi, math.pi) for k in bus_attrs['va']}
    libbus.declare_var_va(model, bus_attrs['names'], initialize=bus_attrs['va'],
                          bounds=va_bounds
                          )

    ###declare the phase angle differences in each branch
    libbranch.declare_var_dva(model, index_set=unique_bus_pairs)

    libbranch.declare_eq_delta_va(model, index_set=unique_bus_pairs)

    ### include the feasibility slack for the bus balances
    p_rhs_kwargs = {}
    q_rhs_kwargs = {}
    if include_feasibility_slack:
        p_rhs_kwargs, q_rhs_kwargs, penalty_expr = _include_feasibility_slack(model, bus_attrs, gen_attrs, bus_p_loads, bus_q_loads)

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
    pt_bounds = s_lbub
    qf_bounds = s_lbub
    qt_bounds = s_lbub
    pf_init = dict()
    pt_init = dict()
    qf_init = dict()
    qt_init = dict()
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
        pt_init[branch_name] = tx_calc.calculate_p(itr_init, itj_init, vr_init[to_bus], vj_init[to_bus])
        qf_init[branch_name] = tx_calc.calculate_q(ifr_init, ifj_init, vr_init[from_bus], vj_init[from_bus])
        qt_init[branch_name] = tx_calc.calculate_q(itr_init, itj_init, vr_init[to_bus], vj_init[to_bus])

    libbranch.declare_var_pf(model=model,
                             index_set=branch_attrs['names'],
                             initialize=pf_init,
                             bounds=pf_bounds
                             )
    libbranch.declare_var_pt(model=model,
                             index_set=branch_attrs['names'],
                             initialize=pt_init,
                             bounds=pt_bounds
                             )
    libbranch.declare_var_qf(model=model,
                             index_set=branch_attrs['names'],
                             initialize=qf_init,
                             bounds=qf_bounds
                             )
    libbranch.declare_var_qt(model=model,
                             index_set=branch_attrs['names'],
                             initialize=qt_init,
                             bounds=qt_bounds
                             )

    ### declare the branch power flow constraints



    ### declare a binary on/off variable for deenergizing a given branch

    decl.declare_var('u', model=model, index_set=branch_attrs['names'], within=pe.Binary)

    model.u.fix(1)

    branch_name_set = decl.declare_set('branch_name', model=model, index_set=branch_attrs['names'])

    model.box_index_set = pe.RangeSet(Q)

    #For active power energization/deenergization
    model.u_branch = pe.Var(branch_name_set, model.box_index_set, within=pe.Binary)

    #For selecting the appropriate interval of the PWL approximation
    model.dva_branch = pe.Var(branch_name_set, model.box_index_set)

    #(5) - Constraints for the on/off variable u

    def u_sum_rule(model, branch_name):
    	return model.u[branch_name] == sum(model.u_branch[branch_name, i] for i in model.box_index_set)

    model.u_sum_Constr = pe.Constraint(branch_name_set, rule=u_sum_rule)


    #(6) - Constraints that sum of dva variables should be equal to total dva

    #Upper bound constraints

    def delta_branch_ub_rule(model, branch_name):
    	branch = branches[branch_name]

    	from_bus = branch['from_bus']
    	to_bus = branch['to_bus']

    	return -model.dva[(from_bus, to_bus)] + sum(model.dva_branch[branch_name, i] for i in model.box_index_set) <= math.pi*(1-model.u[branch_name])

    model.delta_branch_ub_Constr = pe.Constraint(branch_name_set, rule=delta_branch_ub_rule)

    #Lower bound constraints

    def delta_branch_lb_rule(model, branch_name):
    	branch = branches[branch_name]

    	from_bus = branch['from_bus']
    	to_bus = branch['to_bus']

    	return -model.dva[(from_bus, to_bus)] + sum(model.dva_branch[branch_name, i] for i in model.box_index_set) >= -math.pi*(1-model.u[branch_name])

    model.delta_branch_lb_Constr = pe.Constraint(branch_name_set, rule=delta_branch_lb_rule)

    #(7) - Constraints that force dva variable to be in only one interval

    #Upper bound

    def delta_branch_box_ub_rule(model, branch_name, i):
    	delta_ub = branch_dict["Reactive_from_bus"][branch_name]['boxes']['coords'][i-1][7][2]
    	return model.dva_branch[branch_name, i] <= delta_ub*model.u_branch[branch_name, i]

    model.delta_branch_box_ub_Constr = pe.Constraint(branch_name_set, model.box_index_set, rule=delta_branch_box_ub_rule)

    def delta_branch_box_lb_rule(model, branch_name, i):
    	delta_lb = branch_dict["Reactive_from_bus"][branch_name]['boxes']['coords'][i-1][0][2]
    	return model.dva_branch[branch_name, i] >= delta_lb*model.u_branch[branch_name, i]

    model.delta_branch_box_lb_Constr = pe.Constraint(branch_name_set, model.box_index_set, rule=delta_branch_box_lb_rule)

    #(8) - Approximating power flow equation by PWL approximation



    #Active_from_bus
    def pwl_active_from_ub_rule(model, branch_name, i):
    	branch = branches[branch_name]

    	from_bus = branch['from_bus']
    	to_bus = branch['to_bus']

    	g = tx_calc.calculate_conductance(branch)
    	b = tx_calc.calculate_susceptance(branch)
    	coeffs = branch_dict["Active_from_bus"][branch_name]['boxes']['coefficients'][i-1]
    	#M = 10*(g+b) + 4*(coeffs[0]+coeffs[1]+coeffs[2]+coeffs[3])
    	M = 2*s_max[branch_name] + 10*(np.abs(coeffs[0])+np.abs(coeffs[1])+np.abs(coeffs[2])+np.abs(coeffs[3]))
    	return -model.pf[branch_name] + coeffs[0]*model.vm[from_bus] + coeffs[1]*model.vm[to_bus] + coeffs[2]*model.dva_branch[branch_name, i] + coeffs[3] <= M*(1-model.u_branch[branch_name, i])

    model.pwl_active_from_ub_Constr = pe.Constraint(branch_name_set, model.box_index_set, rule=pwl_active_from_ub_rule)

    def pwl_active_from_lb_rule(model, branch_name, i):
    	branch = branches[branch_name]

    	from_bus = branch['from_bus']
    	to_bus = branch['to_bus']

    	g = tx_calc.calculate_conductance(branch)
    	b = tx_calc.calculate_susceptance(branch)

    	coeffs = branch_dict["Active_from_bus"][branch_name]['boxes']['coefficients'][i-1]
    	#M = -(10*(g+b) + 4*(coeffs[0]+coeffs[1]+coeffs[2]+coeffs[3]))
    	M = -2*s_max[branch_name] - 10*(np.abs(coeffs[0])+np.abs(coeffs[1])+np.abs(coeffs[2])+np.abs(coeffs[3]))
    	return -model.pf[branch_name] + coeffs[0]*model.vm[from_bus] + coeffs[1]*model.vm[to_bus] + coeffs[2]*model.dva_branch[branch_name, i] + coeffs[3] >= M*(1-model.u_branch[branch_name, i])

    model.pwl_active_from_lb_Constr = pe.Constraint(branch_name_set, model.box_index_set, rule=pwl_active_from_lb_rule)

    #Active_to_bus
    def pwl_active_to_ub_rule(model, branch_name, i):
    	branch = branches[branch_name]

    	from_bus = branch['from_bus']
    	to_bus = branch['to_bus']

    	g = tx_calc.calculate_conductance(branch)
    	b = tx_calc.calculate_susceptance(branch)

    	coeffs = branch_dict["Active_to_bus"][branch_name]['boxes']['coefficients'][i-1]
    	#M = 10*(g+b) + 4*(coeffs[0]+coeffs[1]+coeffs[2]+coeffs[3])
    	M = 2*s_max[branch_name] + 10*(np.abs(coeffs[0])+np.abs(coeffs[1])+np.abs(coeffs[2])+np.abs(coeffs[3]))
    	return -model.pt[branch_name] + coeffs[0]*model.vm[from_bus] + coeffs[1]*model.vm[to_bus] + coeffs[2]*model.dva_branch[branch_name, i] + coeffs[3] <= M*(1-model.u_branch[branch_name, i])

    model.pwl_active_to_ub_Constr = pe.Constraint(branch_name_set, model.box_index_set, rule=pwl_active_to_ub_rule)

    def pwl_active_to_lb_rule(model, branch_name, i):
    	branch = branches[branch_name]

    	from_bus = branch['from_bus']
    	to_bus = branch['to_bus']

    	g = tx_calc.calculate_conductance(branch)
    	b = tx_calc.calculate_susceptance(branch)

    	coeffs = branch_dict["Active_to_bus"][branch_name]['boxes']['coefficients'][i-1]
    	#M = -(10*(g+b) + 4*(coeffs[0]+coeffs[1]+coeffs[2]+coeffs[3]))
    	M = -2*s_max[branch_name] - 10*(np.abs(coeffs[0])+np.abs(coeffs[1])+np.abs(coeffs[2])+np.abs(coeffs[3]))
    	return -model.pt[branch_name] + coeffs[0]*model.vm[from_bus] + coeffs[1]*model.vm[to_bus] + coeffs[2]*model.dva_branch[branch_name, i] + coeffs[3] >= M*(1-model.u_branch[branch_name, i])

    model.pwl_active_to_lb_Constr = pe.Constraint(branch_name_set, model.box_index_set, rule=pwl_active_to_lb_rule)

    #Reactive_from_bus

    def pwl_reactive_from_ub_rule(model, branch_name, i):
    	branch = branches[branch_name]

    	from_bus = branch['from_bus']
    	to_bus = branch['to_bus']

    	g = tx_calc.calculate_conductance(branch)
    	b = tx_calc.calculate_susceptance(branch)

    	coeffs = branch_dict["Reactive_from_bus"][branch_name]['boxes']['coefficients'][i-1]
    	#M = 10*(g+b) + 4*(coeffs[0]+coeffs[1]+coeffs[2]+coeffs[3])
    	M = 2*s_max[branch_name] + 10*(np.abs(coeffs[0])+np.abs(coeffs[1])+np.abs(coeffs[2])+np.abs(coeffs[3]))
    	return -model.qf[branch_name] + coeffs[0]*model.vm[from_bus] + coeffs[1]*model.vm[to_bus] + coeffs[2]*model.dva_branch[branch_name, i] + coeffs[3] <= M*(1-model.u_branch[branch_name, i])

    model.pwl_reactive_from_ub_Constr = pe.Constraint(branch_name_set, model.box_index_set, rule=pwl_reactive_from_ub_rule)

    def pwl_reactive_from_lb_rule(model, branch_name, i):
    	branch = branches[branch_name]

    	from_bus = branch['from_bus']
    	to_bus = branch['to_bus']

    	g = tx_calc.calculate_conductance(branch)
    	b = tx_calc.calculate_susceptance(branch)

    	coeffs = branch_dict["Reactive_from_bus"][branch_name]['boxes']['coefficients'][i-1]
    	#M = -(10*(g+b) + 4*(coeffs[0]+coeffs[1]+coeffs[2]+coeffs[3]))
    	M = -2*s_max[branch_name] - 10*(np.abs(coeffs[0])+np.abs(coeffs[1])+np.abs(coeffs[2])+np.abs(coeffs[3]))
    	return -model.qf[branch_name] + coeffs[0]*model.vm[from_bus] + coeffs[1]*model.vm[to_bus] + coeffs[2]*model.dva_branch[branch_name, i] + coeffs[3] >= M*(1-model.u_branch[branch_name, i])

    model.pwl_reactive_from_lb_Constr = pe.Constraint(branch_name_set, model.box_index_set, rule=pwl_reactive_from_lb_rule)

    #Reactive_to_bus

    def pwl_reactive_to_ub_rule(model, branch_name, i):
    	branch = branches[branch_name]

    	from_bus = branch['from_bus']
    	to_bus = branch['to_bus']

    	g = tx_calc.calculate_conductance(branch)
    	b = tx_calc.calculate_susceptance(branch)

    	coeffs = branch_dict["Reactive_to_bus"][branch_name]['boxes']['coefficients'][i-1]
    	#M = 10*(g+b) + 4*(coeffs[0]+coeffs[1]+coeffs[2]+coeffs[3])
    	M = 2*s_max[branch_name] + 10*(np.abs(coeffs[0])+np.abs(coeffs[1])+np.abs(coeffs[2])+np.abs(coeffs[3]))
    	return -model.qt[branch_name] + coeffs[0]*model.vm[from_bus] + coeffs[1]*model.vm[to_bus] + coeffs[2]*model.dva_branch[branch_name, i] + coeffs[3] <= M*(1-model.u_branch[branch_name, i])

    model.pwl_reactive_to_ub_Constr = pe.Constraint(branch_name_set, model.box_index_set, rule=pwl_reactive_to_ub_rule)

    def pwl_reactive_to_lb_rule(model, branch_name, i):
    	branch = branches[branch_name]

    	from_bus = branch['from_bus']
    	to_bus = branch['to_bus']

    	g = tx_calc.calculate_conductance(branch)
    	b = tx_calc.calculate_susceptance(branch)

    	coeffs = branch_dict["Reactive_to_bus"][branch_name]['boxes']['coefficients'][i-1]
    	#M = -(10*(g+b) + 4*(coeffs[0]+coeffs[1]+coeffs[2]+coeffs[3]))
    	M = -2*s_max[branch_name] - 10*(np.abs(coeffs[0])+np.abs(coeffs[1])+np.abs(coeffs[2])+np.abs(coeffs[3]))
    	return -model.qt[branch_name] + coeffs[0]*model.vm[from_bus] + coeffs[1]*model.vm[to_bus] + coeffs[2]*model.dva_branch[branch_name, i] + coeffs[3] >= M*(1-model.u_branch[branch_name, i])

    model.pwl_reactive_to_lb_Constr = pe.Constraint(branch_name_set, model.box_index_set, rule=pwl_reactive_to_lb_rule)

    
    ### declare the pq balances
    libbus.declare_eq_p_balance(model=model,
                                index_set=bus_attrs['names'],
                                bus_p_loads=bus_p_loads,
                                gens_by_bus=gens_by_bus,
                                bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                inlet_branches_by_bus=inlet_branches_by_bus,
                                outlet_branches_by_bus=outlet_branches_by_bus,
                                **p_rhs_kwargs
                                )

    libbus.declare_eq_q_balance(model=model,
                                index_set=bus_attrs['names'],
                                bus_q_loads=bus_q_loads,
                                gens_by_bus=gens_by_bus,
                                bus_bs_fixed_shunts=bus_bs_fixed_shunts,
                                inlet_branches_by_bus=inlet_branches_by_bus,
                                outlet_branches_by_bus=outlet_branches_by_bus,
                                **q_rhs_kwargs
                                )

    ### declare the thermal limits
    libbranch.declare_ineq_s_branch_thermal_limit(model=model,
                                                  index_set=branch_attrs['names'],
                                                  branches=branches,
                                                  s_thermal_limits=s_max,
                                                  flow_type=FlowType.POWER
                                                  )

    # declare angle difference limits on interconnected buses
    # libbranch.declare_ineq_angle_diff_branch_lbub_c_s(model=model,
    #                                                   index_set=branch_attrs['names'],
    #                                                   branches=branches
    #                                                   )

    ### declare the generator cost objective
    libgen.declare_expression_pgqg_operating_cost(model=model,
                                                  index_set=gen_attrs['names'],
                                                  p_costs=gen_attrs['p_cost'],
                                                  q_costs=gen_attrs.get('q_cost', None)
                                                  )

    obj_expr = sum(model.pg_operating_cost[gen_name] for gen_name in model.pg_operating_cost)
    if include_feasibility_slack:
        obj_expr += penalty_expr
    if hasattr(model, 'qg_operating_cost'):
        obj_expr += sum(model.qg_operating_cost[gen_name] for gen_name in model.qg_operating_cost)

    model.obj = pe.Objective(expr=obj_expr)

    return model, md


if __name__ == '__main__':
	import os
	from egret.parsers.matpower_parser import create_ModelData

	path = os.path.dirname(__file__)
	case = 'case14_ieee'
	filename = 'pglib_opf_' + case + '.m'
	test_case = os.path.join('c:\\', 'Users', 'wlinz', 'Desktop', 'Restoration', 'Egret', 'egret', 'thirdparty', 'pglib-opf-master', filename) #Better if this isn't so user-dependent
	md_dict = create_ModelData(test_case)
	md = md_dict.clone_in_service()

	
	json_filename = case + '_delta_10_curvature_partition.json'
	with open(json_filename, "r") as read_file:
		branch_dict = json.load(read_file)

	opt = pe.SolverFactory("gurobi")

	MC_model = _create_base_ac_with_pwl_approx_model(md_dict, branch_dict, 10, include_feasibility_slack=False)[0]

	opt.solve(MC_model, tee=True, keepfiles=True, symbolic_solver_labels=True)