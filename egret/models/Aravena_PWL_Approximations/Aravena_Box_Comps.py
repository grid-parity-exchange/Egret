#Implements the inner and outer boxes described by Aravena


import numpy as np

import pyomo.environ as pe

import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc


def power_flow_through_branch(Vi, Vj, delta, branch, bus_type = "from_bus", power_type = "Reactive"):
    if not (power_type =="Active" or power_type =="Reactive"):
        raise ValueError('Power type must be "Active" (for p) or "Reactive" (for q)')

    if not (bus_type == "from_bus" or bus_type == "to_bus"):
        raise ValueError('Bus type must be "from_bus" (for f) or "to_bus" (for t)')

    g = tx_calc.calculate_conductance(branch)
    b = tx_calc.calculate_susceptance(branch)

    if power_type == "Active":
        if bus_type == "from_bus":
        	return Vi**2*g - Vi*Vj*g*pe.cos(delta) - Vi*Vj*b*pe.sin(delta)
        else:
            return Vj**2*g - Vi*Vj*g*pe.cos(delta) - Vi*Vj*b*pe.sin(delta)

    else:
        if bus_type == "from_bus":
    	    return -Vi**2*b + Vi*Vj*b*pe.cos(delta) - Vi*Vj*g*pe.sin(delta)
        else:
            return -Vj**2*b + Vi*Vj*b*pe.cos(delta) - Vi*Vj*g*pe.sin(delta)

def inner_box(branch):
	from_bus = branch['from_bus']
	to_bus = branch['to_bus']

	v_min = bus_attrs['v_min'][from_bus]
	v_max = bus_attrs['v_max'][from_bus]

	s_max = branch['rating_long_term']

	model=pe.ConcreteModel()

	model.Vi_lower = pe.Var(within =pe.Reals, bounds=(v_min, v_max))
	model.Vi_upper = pe.Var(within =pe.Reals, bounds=(v_min, v_max))
	model.Vj_lower = pe.Var(within=pe.Reals, bounds=(v_min, v_max))
	model.Vj_upper = pe.Var(within=pe.Reals,bounds=(v_min, v_max))
	model.delta_lower = pe.Var(within=pe.Reals)
	model.delta_upper = pe.Var(within=pe.Reals)

	def obj_rule(model):
		return pe.log(model.Vi_upper-model.Vi_lower) + pe.log(model.Vj_upper-model.Vi_lower) + pe.log(model.delta_upper-model.delta_lower)

	model.obj = pe.Objective(rule=obj_rule, sense=pe.maximize)

	# #Constraints that each of the vertices of the box must be in the feasible region

	#Vi_lower, Vj_lower, delta_lower
	def f_VilVjldul_rule(model):
		return (power_flow_through_branch(model.Vi_lower, model.Vj_lower, model.delta_lower, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_lower, model.Vj_lower, model.delta_lower, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_VilVjldul_constr = pe.Constraint(rule=f_VilVjldul_rule)

	def t_VilVjldul_rule(model):
		return (power_flow_through_branch(model.Vi_lower, model.Vj_lower, model.delta_lower, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_lower, model.Vj_lower, model.delta_lower, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_VilVjldul_constr = pe.Constraint(rule=t_VilVjldul_rule)

	#Vi_lower, Vj_lower, delta_upper
	def f_VilVjlduu_rule(model):
		return (power_flow_through_branch(model.Vi_lower, model.Vj_lower, model.delta_upper, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_lower, model.Vj_lower, model.delta_upper, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_VilVjlduu_constr = pe.Constraint(rule=f_VilVjlduu_rule)

	def t_VilVjlduu_rule(model):
		return (power_flow_through_branch(model.Vi_lower, model.Vj_lower, model.delta_upper, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_lower, model.Vj_lower, model.delta_upper, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_VilVjlduu_constr = pe.Constraint(rule=t_VilVjlduu_rule)

	#Vi_lower, Vj_upper, delta_lower
	def f_VilVjudul_rule(model):
		return (power_flow_through_branch(model.Vi_lower, model.Vj_upper, model.delta_lower, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_lower, model.Vj_upper, model.delta_lower, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_VilVjudul_constr = pe.Constraint(rule=f_VilVjudul_rule)

	def t_VilVjudul_rule(model):
		return (power_flow_through_branch(model.Vi_lower, model.Vj_upper, model.delta_lower, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_lower, model.Vj_upper, model.delta_lower, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_VilVjudul_constr = pe.Constraint(rule=t_VilVjudul_rule)

	#Vi_lower, Vj_upper, delta_upper
	def f_VilVjuduu_rule(model):
		return (power_flow_through_branch(model.Vi_lower, model.Vj_upper, model.delta_upper, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_lower, model.Vj_upper, model.delta_upper, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_VilVjuduu_constr = pe.Constraint(rule=f_VilVjuduu_rule)

	def t_VilVjuduu_rule(model):
		return (power_flow_through_branch(model.Vi_lower, model.Vj_upper, model.delta_upper, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_lower, model.Vj_upper, model.delta_upper, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_VilVjuduu_constr = pe.Constraint(rule=t_VilVjuduu_rule)

	#Vi_upper, Vj_lower, delta_lower
	def f_ViuVjldul_rule(model):
		return (power_flow_through_branch(model.Vi_upper, model.Vj_lower, model.delta_lower, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_upper, model.Vj_lower, model.delta_lower, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_ViuVjldul_constr = pe.Constraint(rule=f_ViuVjldul_rule)

	def t_ViuVjldul_rule(model):
		return (power_flow_through_branch(model.Vi_upper, model.Vj_lower, model.delta_lower, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_upper, model.Vj_lower, model.delta_lower, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_ViuVjldul_constr = pe.Constraint(rule=t_ViuVjldul_rule)

	#Vi_upper, Vj_lower, delta_upper
	def f_ViuVjlduu_rule(model):
		return (power_flow_through_branch(model.Vi_upper, model.Vj_lower, model.delta_upper, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_upper, model.Vj_lower, model.delta_upper, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_ViuVjlduu_constr = pe.Constraint(rule=f_ViuVjlduu_rule)

	def t_ViuVjlduu_rule(model):
		return (power_flow_through_branch(model.Vi_upper, model.Vj_lower, model.delta_upper, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_upper, model.Vj_lower, model.delta_upper, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_ViuVjlduu_constr = pe.Constraint(rule=t_ViuVjlduu_rule)

	#Vi_upper, Vj_upper, delta_lower
	def f_ViuVjudul_rule(model):
		return (power_flow_through_branch(model.Vi_upper, model.Vj_upper, model.delta_lower, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_upper, model.Vj_upper, model.delta_lower, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_ViuVjudul_constr = pe.Constraint(rule=f_ViuVjudul_rule)

	def t_ViuVjudul_rule(model):
		return (power_flow_through_branch(model.Vi_upper, model.Vj_upper, model.delta_lower, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_upper, model.Vj_upper, model.delta_lower, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_ViuVjudul_constr = pe.Constraint(rule=t_ViuVjudul_rule)

	#Vi_upper, Vj_upper, delta_upper
	def f_ViuVjuduu_rule(model):
		return (power_flow_through_branch(model.Vi_upper, model.Vj_upper, model.delta_upper, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_upper, model.Vj_upper, model.delta_upper, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_ViuVjuduu_constr = pe.Constraint(rule=f_ViuVjuduu_rule)

	def t_ViuVjuduu_rule(model):
		return (power_flow_through_branch(model.Vi_upper, model.Vj_upper, model.delta_upper, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi_upper, model.Vj_upper, model.delta_upper, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_ViuVjuduu_constr = pe.Constraint(rule=t_ViuVjuduu_rule)

	model.pprint()

	opt = pe.SolverFactory("knitroampl")

	opt.solve(model, tee=True)

	return [pe.value(model.Vi_lower), pe.value(model.Vi_upper), pe.value(model.Vj_lower), pe.value(model.Vj_upper), pe.value(model.delta_lower), pe.value(model.delta_upper)]


def ob_Vi_min(v_min, v_max, s_max):
	model=pe.ConcreteModel()

	model.Vi = pe.Var(bounds=(v_min, v_max))
	model.Vj = pe.Var(bounds=(v_min, v_max))
	model.delta = pe.Var()

	#Objective is to minimize arg Vi
	def obj_rule(model):
		return model.Vi

	model.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

	#Constraint is (vi, vj, delta_ij) should be in Omega_ij

	def f_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_omega_constr = pe.Constraint(rule=f_omega_rule)

	def t_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_omega_constr = pe.Constraint(rule=t_omega_rule)

	opt=pe.SolverFactory("knitroampl")

	opt.solve(model, tee=True)

	return pe.value(model.obj)

def ob_Vi_max(v_min, v_max, s_max):
	model=pe.ConcreteModel()

	model.Vi = pe.Var(bounds=(v_min, v_max))
	model.Vj = pe.Var(bounds=(v_min, v_max))
	model.delta = pe.Var()

	#Objective is to minimize arg Vi
	def obj_rule(model):
		return model.Vi

	model.obj = pe.Objective(rule=obj_rule, sense=pe.maximize)

	#Constraint is (vi, vj, delta_ij) should be in Omega_ij

	def f_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_omega_constr = pe.Constraint(rule=f_omega_rule)

	def t_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_omega_constr = pe.Constraint(rule=t_omega_rule)

	opt=pe.SolverFactory("knitroampl")

	opt.solve(model, tee=True)

	return pe.value(model.obj)

def ob_Vj_min(v_min, v_max, s_max):
	model=pe.ConcreteModel()

	model.Vi = pe.Var(bounds=(v_min, v_max))
	model.Vj = pe.Var(bounds=(v_min, v_max))
	model.delta = pe.Var()

	#Objective is to minimize arg Vi
	def obj_rule(model):
		return model.Vj

	model.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

	#Constraint is (vi, vj, delta_ij) should be in Omega_ij

	def f_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_omega_constr = pe.Constraint(rule=f_omega_rule)

	def t_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_omega_constr = pe.Constraint(rule=t_omega_rule)

	opt=pe.SolverFactory("knitroampl")

	opt.solve(model, tee=True)

	return pe.value(model.obj)

def ob_Vj_max(v_min, v_max, s_max):
	model=pe.ConcreteModel()

	model.Vi = pe.Var(bounds=(v_min, v_max))
	model.Vj = pe.Var(bounds=(v_min, v_max))
	model.delta = pe.Var()

	#Objective is to minimize arg Vi
	def obj_rule(model):
		return model.Vj

	model.obj = pe.Objective(rule=obj_rule, sense=pe.maximize)

	#Constraint is (vi, vj, delta_ij) should be in Omega_ij

	def f_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_omega_constr = pe.Constraint(rule=f_omega_rule)

	def t_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_omega_constr = pe.Constraint(rule=t_omega_rule)

	opt=pe.SolverFactory("knitroampl")

	opt.solve(model, tee=True)

	return pe.value(model.obj)

def ob_delta_min(v_min, v_max, s_max):
	model=pe.ConcreteModel()

	model.Vi = pe.Var(bounds=(v_min, v_max))
	model.Vj = pe.Var(bounds=(v_min, v_max))
	model.delta = pe.Var()

	#Objective is to minimize arg Vi
	def obj_rule(model):
		return model.delta

	model.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

	#Constraint is (vi, vj, delta_ij) should be in Omega_ij

	def f_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_omega_constr = pe.Constraint(rule=f_omega_rule)

	def t_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_omega_constr = pe.Constraint(rule=t_omega_rule)

	opt=pe.SolverFactory("knitroampl")

	opt.solve(model, tee=True)

	return pe.value(model.obj)

def ob_delta_max(v_min, v_max, s_max):
	model=pe.ConcreteModel()

	model.Vi = pe.Var(bounds=(v_min, v_max))
	model.Vj = pe.Var(bounds=(v_min, v_max))
	model.delta = pe.Var()

	#Objective is to minimize arg Vi
	def obj_rule(model):
		return model.delta

	model.obj = pe.Objective(rule=obj_rule, sense=pe.maximize)

	#Constraint is (vi, vj, delta_ij) should be in Omega_ij

	def f_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="from_bus", power_type="Reactive"))**2 <= s_max**2
	model.f_omega_constr = pe.Constraint(rule=f_omega_rule)

	def t_omega_rule(model):
		return (power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Active"))**2 + \
				(power_flow_through_branch(model.Vi, model.Vj, model.delta, branch, bus_type="to_bus", power_type="Reactive"))**2 <= s_max**2
	model.t_omega_constr = pe.Constraint(rule=t_omega_rule)

	opt=pe.SolverFactory("knitroampl")

	opt.solve(model, tee=True)

	return pe.value(model.obj)


def outer_box(branch):
	from_bus = branch['from_bus']
	to_bus = branch['to_bus']

	v_min = bus_attrs['v_min'][from_bus]
	v_max = bus_attrs['v_max'][from_bus]

	s_max = branch['rating_long_term']

	Vi_lower = ob_Vi_min(v_min, v_max, s_max)
	Vi_upper = ob_Vi_max(v_min, v_max, s_max)
	Vj_lower = ob_Vj_min(v_min, v_max, s_max)
	Vj_upper = ob_Vj_max(v_min, v_max, s_max)
	delta_lower = max(ob_delta_min(v_min, v_max, s_max), branch['angle_diff_min']*np.pi/180)
	delta_upper = min(ob_delta_max(v_min, v_max, s_max), branch['angle_diff_max']*np.pi/180)

	return [Vi_lower, Vi_upper, Vj_lower, Vj_upper, delta_lower, delta_upper]


if __name__ == '__main__':
    import os
    from egret.parsers.matpower_parser import create_ModelData

    path = os.path.dirname(__file__)
    power_types = ["Reactive", "Active"]
    bus_types = ["from_bus", "to_bus"]
    case = 'case14_ieee'
    filename = 'pglib_opf_' + case + '.m'
    test_case = os.path.join('c:\\', 'Users', 'wlinz', 'Desktop', 'Restoration', 'Egret', 'egret', 'thirdparty', 'pglib-opf-master', filename) #Better if this isn't so user-dependent
    md_dict = create_ModelData(test_case)
    md = md_dict.clone_in_service()

    branches = dict(md.elements(element_type='branch'))
    branch_attrs = md.attributes(element_type='branch')

    buses = dict(md.elements(element_type='bus'))
    bus_attrs = md.attributes(element_type='bus')

    branch = branches['1']

    print(branch)

    inner_box = inner_box(branch)

    outer_box = outer_box(branch)

    print(inner_box)

    print(outer_box)

