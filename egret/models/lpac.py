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
from math import pi, radians, sqrt

from egret.models.acopf import _include_feasibility_slack, solve_acopf, create_psv_acopf_model, create_rsv_acopf_model

def declare_pwl_cosine_bounds(model, index_set, branches, lower_bound, upper_bound, cosine_segment_count):
	"""
	Add piecewise linear constraints for the cosine hat variables. The constraints are indexed over the branches and the number of cosine segments. 
	"""
	cons_index_set = pe.Set(initialize = model.N * index_set)

	increment = (upper_bound - lower_bound)/(cosine_segment_count + 1)
	def pwl_cosine_rule(model, i, branch_name):
		branch = branches[branch_name]

		from_bus = branch['from_bus']
		to_bus = branch['to_bus']
		if i==0:
			return model.cos_hat[branch_name] - ((pe.cos(upper_bound) - pe.cos(lower_bound))/(upper_bound - lower_bound) * ((model.va[from_bus] - model.va[to_bus]) - lower_bound) + pe.cos(lower_bound)) >= 0
		else:
			a = lower_bound + i*increment
			return 0 <= (-pe.sin(a)*(model.va[from_bus] - model.va[to_bus] - a) + pe.cos(a)) - model.cos_hat[branch_name]

	model.pwl_cosine_bounds = pe.Constraint(cons_index_set, rule=pwl_cosine_rule)


def create_hot_start_lpac_model(model_data, voltages, lower_bound = -pi/3, upper_bound = pi/3, cosine_segment_count = 20, include_feasibility_slack = False):
	"""
	The hot start LPAC model assumes that voltages are known, e.g. from an AC base point solution.
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

	###declare (and fix) the voltage magnitudes and squares of voltage magnitudes

	bus_voltage_magnitudes = voltages #Assumes voltages is given as a dictionary 
	libbus.declare_var_vm(model, bus_attrs['names'], initialize=bus_voltage_magnitudes)
	model.vm.fix()

	libbus.declare_var_vmsq(model=model,
                            index_set=bus_attrs['names'],
                            initialize={k: v**2 for k, v in bus_attrs['vm'].items()},
                            bounds=zip_items({k: v**2 for k, v in bus_attrs['v_min'].items()},
                                             {k: v**2 for k, v in bus_attrs['v_max'].items()}))

    ### declare the polar voltages

	libbus.declare_var_va(model, bus_attrs['names'], initialize=bus_attrs['va'])

    ### declare the cosine approximation variables
	cos_hat_bounds = {k: (0, 1) for k in branch_attrs['names']}
	decl.declare_var('cos_hat', model, branch_attrs['names'], bounds = cos_hat_bounds)

    ### fix the reference bus
	ref_bus = md.data['system']['reference_bus']
    #ref_angle = md.data['system']['reference_bus_angle']
	model.va[ref_bus].fix(radians(0.0))

	### declare the fixed shunts at the buses
	bus_bs_fixed_shunts, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    ### declare (and fix) the loads at the buses
	bus_p_loads, bus_q_loads = tx_utils.dict_of_bus_loads(buses, loads)

	libbus.declare_var_pl(model, bus_attrs['names'], initialize=bus_p_loads)
	libbus.declare_var_ql(model, bus_attrs['names'], initialize=bus_q_loads)
	model.pl.fix()
	model.ql.fix()

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


    ####################
    #Constraints
    ####################

    ###Balance equations in a bus 

    #p balance

	libbus.declare_eq_p_balance(model=model,
                                          index_set=bus_attrs['names'],
                                          bus_p_loads=bus_p_loads,
                                          gens_by_bus=gens_by_bus,
                                          bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                          inlet_branches_by_bus=inlet_branches_by_bus,
                                          outlet_branches_by_bus=outlet_branches_by_bus,
                                          **p_rhs_kwargs
                                          )

    #q balance

	libbus.declare_eq_q_balance(model=model, index_set=bus_attrs['names'],
                         		bus_q_loads=bus_q_loads,
                         		gens_by_bus=gens_by_bus,
                         		bus_bs_fixed_shunts=bus_bs_fixed_shunts,
                         		inlet_branches_by_bus=inlet_branches_by_bus, 
                         		outlet_branches_by_bus=outlet_branches_by_bus,
                        		 **q_rhs_kwargs)


    ### Power in a branch

	branch_con_set = decl.declare_set('_con_eq_p_q_lpac_branch_power', model, branch_attrs['names'])

	model.eq_pf_branch_t = pe.Constraint(branch_con_set)
	model.eq_pt_branch_t = pe.Constraint(branch_con_set)
	model.eq_qf_branch_t = pe.Constraint(branch_con_set)
	model.eq_qt_branch_t = pe.Constraint(branch_con_set)

	for branch_name in branch_con_set:
		branch = branches[branch_name]

		from_bus = branch['from_bus']
		to_bus = branch['to_bus']

		g = tx_calc.calculate_conductance(branch)
		b = tx_calc.calculate_susceptance(branch)

		model.eq_pf_branch_t[branch_name] = \
			model.pf[branch_name] == \
			g*model.vmsq[from_bus] - model.vm[from_bus]*model.vm[to_bus]*(g * model.cos_hat[branch_name] + b * (model.va[from_bus] - model.va[to_bus]))

		model.eq_pt_branch_t[branch_name] = \
			model.pt[branch_name] == \
			g*model.vmsq[to_bus] - model.vm[from_bus]*model.vm[to_bus]*(g * model.cos_hat[branch_name] + b * (model.va[to_bus] - model.va[from_bus]))

		model.eq_qf_branch_t[branch_name] = \
			model.qf[branch_name] == \
			-b*model.vmsq[from_bus] - model.vm[from_bus]*model.vm[to_bus]*(g*(model.va[from_bus] - model.va[to_bus]) - b*model.cos_hat[branch_name])

		model.eq_qt_branch_t[branch_name] = \
			model.qt[branch_name] == \
			-b*model.vmsq[to_bus] - model.vm[from_bus]*model.vm[to_bus]*(g*(model.va[to_bus] - model.va[from_bus]) - b*model.cos_hat[branch_name])


    ### Piecewise linear cosine constraints

	model.N = pe.Set(initialize=list(range(cosine_segment_count+1)))

	declare_pwl_cosine_bounds(model = model, 
    							index_set = branch_attrs['names'],
    							 branches=branches, 
    							 lower_bound=lower_bound, 
    							 upper_bound=upper_bound, 
    							 cosine_segment_count=cosine_segment_count)



    ### Objective is to maximize cosine hat variables

	obj_expr = sum(model.cos_hat[branch_name] for branch_name in branch_attrs['names'])

	if include_feasibility_slack:
		obj_expr += penalty_expr

	model.obj = pe.Objective(expr=obj_expr)

	return model, md


def create_warm_start_lpac_model(model_data, voltages, lower_bound = -pi/3, upper_bound = pi/3, cosine_segment_count = 20, include_feasibility_slack = False):
	"""
	The warm start LPAC model assumes that target voltages can be given, but not all voltages are known. 
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

	###declare (and fix) the target voltage magnitudes and squares of voltage magnitudes
	bus_voltage_magnitudes = voltages
	libbus.declare_var_vm(model, bus_attrs['names'], initialize=bus_voltage_magnitudes)
	model.vm.fix()

	libbus.declare_var_vmsq(model=model,
                            index_set=bus_attrs['names'],
                            initialize={k: v**2 for k, v in bus_attrs['vm'].items()},
                            bounds=zip_items({k: v**2 for k, v in bus_attrs['v_min'].items()},
                                             {k: v**2 for k, v in bus_attrs['v_max'].items()}))

    ### declare the polar voltages
	libbus.declare_var_va(model, bus_attrs['names'], initialize=bus_attrs['va'] )

    ### declare the voltage change variables
	decl.declare_var('phi', model, bus_attrs['names'])

    ### declare the cosine approximation variables
	cos_hat_bounds = {k: (0, 1) for k in branch_attrs['names']}
	decl.declare_var('cos_hat', model, branch_attrs['names'], bounds = cos_hat_bounds)

    ### fix the reference bus
	ref_bus = md.data['system']['reference_bus']
    #ref_angle = md.data['system']['reference_bus_angle']
	model.va[ref_bus].fix(radians(0.0))

    ### declare the fixed shunts at the buses
	bus_bs_fixed_shunts, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    ### declare (and fix) the loads at the buses
	bus_p_loads, bus_q_loads = tx_utils.dict_of_bus_loads(buses, loads)

	libbus.declare_var_pl(model, bus_attrs['names'], initialize=bus_p_loads)
	libbus.declare_var_ql(model, bus_attrs['names'], initialize=bus_q_loads)
	model.pl.fix()
	model.ql.fix()

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


    ########################
    #Constraints
    ########################


    ###Balance equations at a bus

    #p balance

	libbus.declare_eq_p_balance(model=model,
                                          index_set=bus_attrs['names'],
                                          bus_p_loads=bus_p_loads,
                                          gens_by_bus=gens_by_bus,
                                          bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                          inlet_branches_by_bus=inlet_branches_by_bus,
                                          outlet_branches_by_bus=outlet_branches_by_bus,
                                          **p_rhs_kwargs
                                          )

    #q balance

	libbus.declare_eq_q_balance(model=model, index_set=bus_attrs['names'],
                         		bus_q_loads=bus_q_loads,
                         		gens_by_bus=gens_by_bus,
                         		bus_bs_fixed_shunts=bus_bs_fixed_shunts,
                         		inlet_branches_by_bus=inlet_branches_by_bus, 
                         		outlet_branches_by_bus=outlet_branches_by_bus,
                        		 **q_rhs_kwargs)


    ###Constraints for power in a branch

	branch_con_set = decl.declare_set('_con_eq_p_q_lpac_branch_power', model, branch_attrs['names'])

	model.eq_pf_branch_t = pe.Constraint(branch_con_set)
	model.eq_pt_branch_t = pe.Constraint(branch_con_set)
	model.eq_qf_branch_t = pe.Constraint(branch_con_set)
	model.eq_qt_branch_t = pe.Constraint(branch_con_set)

	for branch_name in branch_con_set:
		branch = branches[branch_name]

		from_bus = branch['from_bus']
		to_bus = branch['to_bus']

		g = tx_calc.calculate_conductance(branch)
		b = tx_calc.calculate_susceptance(branch)

		model.eq_pf_branch_t[branch_name] = \
			model.pf[branch_name] == \
			g*model.vmsq[from_bus] - model.vm[from_bus]*model.vm[to_bus]*(g * model.cos_hat[branch_name] + b * (model.va[from_bus] - model.va[to_bus]))

		model.eq_pt_branch_t[branch_name] = \
			model.pt[branch_name] == \
			g*model.vmsq[to_bus] - model.vm[from_bus]*model.vm[to_bus]*(g * model.cos_hat[branch_name] + b * (model.va[to_bus] - model.va[from_bus]))

		model.eq_qf_branch_t[branch_name] = \
			model.qf[branch_name] == \
			-b*model.vmsq[from_bus] - model.vm[from_bus]*model.vm[to_bus]*(g*(model.va[from_bus] - model.va[to_bus]) - b*model.cos_hat[branch_name]) - model.vm[from_bus]*b*(model.phi[from_bus] - model.phi[to_bus]) - (model.vm[from_bus] - model.vm[to_bus])*model.phi[from_bus]

		model.eq_qt_branch_t[branch_name] = \
			model.qt[branch_name] == \
			-b*model.vmsq[to_bus] - model.vm[from_bus]*model.vm[to_bus]*(g*(model.va[to_bus] - model.va[from_bus]) - b*model.cos_hat[branch_name]) - model.vm[to_bus]*b*(model.phi[to_bus] - model.phi[from_bus]) - (model.vm[to_bus] - model.vm[from_bus])*model.phi[to_bus]

    ### Piecewise linear cosine constraints

	model.N = pe.Set(initialize=list(range(cosine_segment_count+1)))

	declare_pwl_cosine_bounds(model = model, 
    							index_set = branch_attrs['names'],
    							 branches=branches, 
    							 lower_bound=lower_bound, 
    							 upper_bound=upper_bound, 
    							 cosine_segment_count=cosine_segment_count)



    ### Objective is to maximize cosine hat variables

	obj_expr = sum(model.cos_hat[branch_name] for branch_name in branch_attrs['names'])

	if include_feasibility_slack:
		obj_expr += penalty_expr

	model.obj = pe.Objective(expr=obj_expr)

	return model, md


def create_cold_start_lpac_model(model_data, cosine_segment_count = 20, lower_bound = -pi/3, upper_bound = pi/3, include_feasibility_slack = False):
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

	libbus.declare_var_vmsq(model=model,
                            index_set=bus_attrs['names'],
                            initialize={k: v**2 for k, v in bus_attrs['vm'].items()},
                            bounds=zip_items({k: v**2 for k, v in bus_attrs['v_min'].items()},
                                             {k: v**2 for k, v in bus_attrs['v_max'].items()}))

    ### declare the voltage change variables
	decl.declare_var('phi', model, bus_attrs['names'])

    ### declare the cosine approximation variables
	cos_hat_bounds = {k: (0, 1) for k in branch_attrs['names']}
	decl.declare_var('cos_hat', model, branch_attrs['names'], bounds = cos_hat_bounds)

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

	### declare the fixed shunts at the buses
	bus_bs_fixed_shunts, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

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



	################################
	#Constraints
	################################

    ### Balance equations at a bus (based on Kirchhoff Current Law)

    #Should be able to just use DC OPF approximation of B-theta type? 

    ### declare the p balance
	libbus.declare_eq_p_balance(model=model,
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

	libbus.declare_eq_q_balance(model=model, index_set=bus_attrs['names'],
                         		bus_q_loads=bus_q_loads,
                         		gens_by_bus=gens_by_bus,
                         		bus_bs_fixed_shunts=bus_bs_fixed_shunts,
                         		inlet_branches_by_bus=inlet_branches_by_bus, 
                         		outlet_branches_by_bus=outlet_branches_by_bus,
                        		 **q_rhs_kwargs)

	

    ### Constraints for power in a branch 

	branch_con_set = decl.declare_set('_con_eq_p_q_lpac_branch_power', model, branch_attrs['names'])

	model.eq_pf_branch_t = pe.Constraint(branch_con_set)
	model.eq_pt_branch_t = pe.Constraint(branch_con_set)
	model.eq_qf_branch_t = pe.Constraint(branch_con_set)
	model.eq_qt_branch_t = pe.Constraint(branch_con_set)

	for branch_name in branch_con_set:
		branch = branches[branch_name]

		from_bus = branch['from_bus']
		to_bus = branch['to_bus']

		g = tx_calc.calculate_conductance(branch)
		b = tx_calc.calculate_susceptance(branch)

		model.eq_pf_branch_t[branch_name] = \
        	model.pf[branch_name] == \
        	g - g * model.cos_hat[branch_name] - b * (model.va[from_bus] - model.va[to_bus])

		model.eq_pt_branch_t[branch_name] = \
        	model.pt[branch_name] == \
        	g - g * model.cos_hat[branch_name] - b * (model.va[to_bus] - model.va[from_bus])

		model.eq_qf_branch_t[branch_name] = \
        	model.qf[branch_name] == \
        	-b - g*(model.va[from_bus] - model.va[to_bus]) + b*model.cos_hat[branch_name] - b*(model.phi[from_bus] - model.phi[to_bus])

		model.eq_qt_branch_t[branch_name] = \
        	model.qt[branch_name] == \
        	-b - g*(model.va[to_bus] - model.va[from_bus]) +b*model.cos_hat[branch_name] - b*(model.phi[to_bus] - model.phi[from_bus])

    ### Piecewise linear cosine constraints

	model.N = pe.Set(initialize=list(range(cosine_segment_count+1)))

	declare_pwl_cosine_bounds(model = model, 
    							index_set = branch_attrs['names'],
    							 branches=branches, 
    							 lower_bound=lower_bound, 
    							 upper_bound=upper_bound, 
    							 cosine_segment_count=cosine_segment_count)
	
    ### Objective is to maximize cosine hat variables

	# obj_expr = sum(model.cos_hat[branch_name] for branch_name in branch_attrs['names'])

	# if include_feasibility_slack:
	# 	obj_expr += penalty_expr

	# model.obj = pe.Objective(expr=obj_expr)

	###Objective to match with acopf.py

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


def solve_lpac(model_data,
                solver,
                ac_solver = None,
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
    ac_solver : str or pyomo.opt.base.solvers.OptSolver (optional)
    	Either a string specifying a pyomo solver name, or an instantiated pyomo solver.
    	Default is None for the cold start lpac model. 
    timelimit : float (optional)
        Time limit for dcopf run. Default of None results in no time
        limit being set.
    solver_tee : bool (optional)
        Display solver log. Default is True.
    symbolic_solver_labels : bool (optional)
        Use symbolic solver labels. Useful for debugging; default is False.
    options : dict (optional)
        Other options to pass into the solver. Default is dict().
    lpac_model_generator : function (optional)
        Function for generating the lpac model. Default is
        the cold start lpac model
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


    if lpac_model_generator == create_hot_start_lpac_model or lpac_model_generator == create_warm_start_lpac_model:
    	if ac_solver != None:
    		ac_md, ac_m, ac_results = solve_acopf(model_data, ac_solver, options=options,acopf_model_generator=create_psv_acopf_model,return_model=True, return_results=True,**kwargs)
    	else:
    		ac_md, ac_m, ac_results = solve_acopf(model_data, solver, options=options,acopf_model_generator=create_psv_acopf_model,return_model=True, return_results=True,**kwargs)
    	voltages = dict({})
    	for bus in ac_md.elements(element_type="bus"):
    		voltages[bus[0]] = bus[1]['vm']
    	#print(voltages)
    	m, md = lpac_model_generator(model_data, voltages, **kwargs)
    else:
    	m, md = lpac_model_generator(model_data, **kwargs)


    m, results, solver = _solve_model(m, solver, timelimit=timelimit, solver_tee=solver_tee, \
    									symbolic_solver_labels = symbolic_solver_labels, solver_options=options, return_solver=True)



	# save results data to ModelData object

    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))

    md.data['system']['total_cost'] = value(m.obj)

    for g,g_dict in gens.items():
	    g_dict['pg'] = value(m.pg[g])
	    g_dict['qg'] = value(m.qg[g])

    for b,b_dict in buses.items():
	#b_dict['lmp'] = value(m.dual[m.eq_p_balance[b]])
	#b_dict['qlmp'] = value(m.dual[m.eq_q_balance[b]])
        b_dict['pl'] = value(m.pl[b])
        #if hasattr(m, 'vj'):
            #b_dict['vm'] = tx_calc.calculate_vm_from_vj_vr(value(m.vj[b]), value(m.vr[b]))
            #b_dict['va'] = tx_calc.calculate_va_from_vj_vr(value(m.vj[b]), value(m.vr[b]))
        #else:
            #b_dict['vm'] = value(m.vm[b])
            #b_dict['va'] = value(m.va[b])

    for k, k_dict in branches.items():
        if hasattr(m,'pf'):
            k_dict['pf'] = value(m.pf[k])
            k_dict['pt'] = value(m.pt[k])
            k_dict['qf'] = value(m.qf[k])
            k_dict['qt'] = value(m.qt[k])
        if hasattr(m,'irf'):
            b = k_dict['from_bus']
            k_dict['pf'] = value(tx_calc.calculate_p(value(m.ifr[k]), value(m.ifj[k]), value(m.vr[b]), value(m.vj[b])))
            k_dict['qf'] = value(tx_calc.calculate_q(value(m.ifr[k]), value(m.ifj[k]), value(m.vr[b]), value(m.vj[b])))
            b = k_dict['to_bus']
            k_dict['pt'] = value(tx_calc.calculate_p(value(m.itr[k]), value(m.itj[k]), value(m.vr[b]), value(m.vj[b])))
            k_dict['qt'] = value(tx_calc.calculate_q(value(m.itr[k]), value(m.itj[k]), value(m.vr[b]), value(m.vj[b])))


    unscale_ModelData_to_pu(md, inplace=True)

    #print(buses)
    #print(gens)
    # print(branches)

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

    filename = 'pglib_opf_case14_ieee.m'
    test_case = os.path.join('c:\\', 'Users', 'wlinz', 'Desktop', 'Restoration', 'Egret', 'egret', 'thirdparty', 'pglib-opf-master', filename) #Better if this isn't so user-dependent
    model_data = create_ModelData(test_case)
    baron_options = {'LPSol': 3, 'CplexLibName': "C:/Program Files/IBM/ILOG/CPLEX_Studio129/cplex/bin/x64_win64/cplex1290.dll", 'summary': 1, 'SumName': "summary"}
    knitro_options = {'maxit': 20000}
    kwargs = {'include_feasibility_slack':False}
    #md,m,results = solve_lpac(model_data, "baron", lpac_model_generator=create_cold_start_lpac_model,return_model=True, return_results=True,**kwargs)
    md,m,results = solve_lpac(model_data, "cplex", ac_solver = "knitroampl",lpac_model_generator=create_hot_start_lpac_model,return_model=True, return_results=True,**kwargs)
    md,m,results = solve_lpac(model_data, "knitroampl",lpac_model_generator=create_warm_start_lpac_model,return_model=True, return_results=True,**kwargs)
    md,m,results = solve_lpac(model_data, "knitroampl",options=knitro_options,lpac_model_generator=create_cold_start_lpac_model,return_model=True, return_results=True,**kwargs)
    #ac_md, ac_m, ac_results = solve_acopf(model_data, "knitroampl", options=knitro_options,acopf_model_generator=create_rsv_acopf_model,return_model=True, return_results=True,**kwargs)
    
    
    


    
