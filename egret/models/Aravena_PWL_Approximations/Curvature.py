#Code to produce values for cumulated curvature of power equations

import numpy as np

import scipy.integrate as integrate

import egret.model_library.transmission.tx_calc as tx_calc

import matplotlib.pyplot as plt


#############
#Functions for Curvature for Power Equations
#############

#Assumes delta can be changed, but Vi, Vj are fixed
def power_equation(delta, branch, bus_type = "from_bus", power_type = "Reactive"):
	if not (power_type =="Active" or power_type =="Reactive"):
		raise ValueError('Power type must be "Active" (for p) or "Reactive" (for q)')

	if not (bus_type == "from_bus" or bus_type == "to_bus"):
		raise ValueError('Bus type must be "from_bus" (for f) or "to_bus" (for t)')

	g = tx_calc.calculate_conductance(branch)
	b = tx_calc.calculate_susceptance(branch)

	
	if power_type == "Active":
		return g - g*np.cos(delta) - b*np.sin(delta)
	else:
		return -b + b*np.cos(delta) - g*np.sin(delta)

def power_deriv(delta, branch, bus_type = "from_bus", power_type = "Reactive"):
	if not (power_type =="Active" or power_type =="Reactive"):
	    raise ValueError('Power type must be "Active" (for p) or "Reactive" (for q)')

	if not (bus_type == "from_bus" or bus_type == "to_bus"):
	    raise ValueError('Bus type must be "from_bus" (for f) or "to_bus" (for t)')

	g = tx_calc.calculate_conductance(branch)
	b = tx_calc.calculate_susceptance(branch)

	
	if power_type == "Active":
		return g*np.sin(delta) - b*np.cos(delta)
	else:
		return -b*np.sin(delta) - g*np.cos(delta)

def power_second_deriv(delta, branch, bus_type = "from_bus", power_type = "Reactive"):
	g = tx_calc.calculate_conductance(branch)
	b = tx_calc.calculate_susceptance(branch)

	
	if power_type == "Active":
		return g*np.cos(delta) + b*np.sin(delta)
	else:
		return -b*np.cos(delta) + g*np.sin(delta)

def curvature(delta, branch, bus_type = "from_bus", power_type = "Reactive"):
	return np.absolute(power_second_deriv(delta, branch, bus_type, power_type))/(1+(power_deriv(delta, branch, bus_type, power_type))**2)

#########################
#Functions for Partitioning the domain into Q pieces of equal curvature
#########################

def curvature_target_x_value(integration_lb, interval_lb, interval_ub, target_value, branch, bus_type, power_type="Reactive", eps = 0.0001):
	#Finds an approximation for a value x in the domain where accumulated curvature is equal to a target value
	#The strategy is to recursively narrow down the given interval into smaller and smaller intervals until an acceptable point is found. 
	x = np.linspace(interval_lb, interval_ub)
	i = 0
	for value in x:
		cumulated_curvature = integrate.quad(curvature, integration_lb, value, args=(branch, bus_type, power_type))[0]
		if np.abs(cumulated_curvature - target_value) <= eps:
			return value
		if target_value - cumulated_curvature > eps:
			i = i + 1
			continue
		if cumulated_curvature - target_value > eps:
			return curvature_target_x_value(integration_lb, x[i-1], x[i], target_value, branch, bus_type, power_type, eps)

def eq_curvature_partition(lb, ub, Q, branch, bus_type, power_type="Reactive", eps = 0.0001):
	#Divides the domain [lb, ub] into Q pieces of (approximately) equal curvature
	if Q == 1:
		return [lb, ub]
	breakpoints = [lb]
	target_value = (integrate.quad(curvature, lb, ub, args=(branch, bus_type, power_type))[0])/Q
	for i in range(Q-1):
		breakpoints.append(curvature_target_x_value(breakpoints[i], breakpoints[i], ub, target_value, branch, bus_type, power_type, eps))
	breakpoints.append(ub)
	return breakpoints

def dev_from_linear(lb, ub, branch, bus_type, power_type="Reactive", lin_tol=0.1):
	#Finds an x-value where deviation from normal is more than an epsilon tolerance
	x = np.linspace(lb, ub)
	slope = (power_equation(x[1], branch, bus_type, power_type) - power_equation(x[0], branch, bus_type, power_type))/(x[1] - x[0])
	i = 0
	for value in x: 
		if np.abs(power_equation(value, branch, bus_type, power_type) - (slope*(value - x[0])+power_equation(x[0], branch, bus_type, power_type))) <= lin_tol:
			i = i + 1
			continue
		else:
			return x[i-1]
	return ub

def close_to_linear_cuts(lb, ub, branch, bus_type, power_type="Reactive", lin_tol=0.1):
	#Refines a given interval based on how much the power_equation deviates from the linear approximation. 
	partition = [lb]
	x = dev_from_linear(lb, ub, branch, bus_type, power_type, lin_tol)
	if x == ub:
		partition.append(x)
		return partition
	else:
		return partition + close_to_linear_cuts(x, ub, branch, bus_type, power_type, lin_tol)
	

def refined_eq_curvature_partition(lb, ub, Q, branch, bus_type, power_type = "Reactive", eps=0.0001, lin_tol=0.1):
	partition = eq_curvature_partition(lb, ub, Q, branch, bus_type, power_type, eps)
	refined_partition = [lb]
	for i in range(Q):
		refined_partition = refined_partition + [partition[i]] + close_to_linear_cuts(partition[i], partition[i+1], branch, bus_type, power_type, lin_tol) + [partition[i+1]]
	no_dup_refined_partition = []
	for j in range(len(refined_partition) - 1):
		if refined_partition[j] == refined_partition[j+1]:
			continue
		else:
			no_dup_refined_partition.append(refined_partition[j])
	no_dup_refined_partition.append(ub)
	return no_dup_refined_partition

if __name__ == '__main__':
    import os
    import pyomo.environ as pe
    import egret.model_library.transmission.tx_utils as tx_utils
    import egret.model_library.transmission.tx_calc as tx_calc
    import egret.model_library.transmission.bus as libbus
    import egret.model_library.transmission.branch as libbranch
    import egret.model_library.transmission.gen as libgen
    from egret.parsers.matpower_parser import create_ModelData

    from egret.data.data_utils import map_items, zip_items

    path = os.path.dirname(__file__)
    filename = 'pglib_opf_case30_ieee.m'
    test_case = os.path.join('c:\\', 'Users', 'wlinz', 'Desktop', 'Restoration', 'Egret', 'egret', 'thirdparty', 'pglib-opf-master', filename) #Better if this isn't so user-dependent
    md_dict = create_ModelData(test_case)
    md = md_dict.clone_in_service()

    branches = dict(md.elements(element_type='branch'))
    branch_attrs = md.attributes(element_type='branch')

    buses = dict(md.elements(element_type='bus'))
    bus_attrs = md.attributes(element_type='bus')

    model=pe.ConcreteModel()

    ### declare the polar voltages
    libbus.declare_var_va(model, bus_attrs['names'], initialize=bus_attrs['va'])

    libbus.declare_var_vm(model, bus_attrs['names'], initialize=bus_attrs['vm'])

    libbus.declare_var_vmsq(model=model,
                            index_set=bus_attrs['names'],
                            initialize={k: v**2 for k, v in bus_attrs['vm'].items()},
                            bounds=zip_items({k: v**2 for k, v in bus_attrs['v_min'].items()},
                                             {k: v**2 for k, v in bus_attrs['v_max'].items()}))


   	#print(branches)
    #print(branches.keys())

    branch = branches['10']
    print(branch)

    from_partition = eq_curvature_partition(-np.pi/6, np.pi/6, 20, branch, bus_type="from_bus")

    to_partition = eq_curvature_partition(-np.pi/6, np.pi/6, 20, branch, bus_type="to_bus")

    from_refined_partition = refined_eq_curvature_partition(-np.pi/6, np.pi/6, 20, branch, bus_type="from_bus", power_type="Reactive", eps=0.0001, lin_tol = 0.5)

    to_refined_partition = refined_eq_curvature_partition(-np.pi/6, np.pi/6, 20, branch, bus_type="to_bus", power_type="Reactive", eps=0.0001, lin_tol = 0.5)

    print("Unrefined from partition:", from_partition)

    print("Unrefined to partition:", to_partition)

    print("Refined from partition:", from_refined_partition)

    print("Refined to partition: ", to_refined_partition)

    x = np.linspace(-np.pi/6, np.pi/6)

    y1 = power_equation(x, branch)

    #plt.plot(x, y1)
    #plt.plot(x, y2)
    #plt.show()

    for i in range(len(from_partition) - 1):
    	x = np.linspace(from_partition[i], from_partition[i+1])
    	y = power_equation(x, branch)
    	plt.plot(x, y)

    plt.show()

    for i in range(len(from_refined_partition) - 1):
    	x = np.linspace(from_refined_partition[i], from_refined_partition[i+1])
    	y = power_equation(x, branch)
    	plt.plot(x, y)

    plt.show()

	




###########
#Plotting
##########

#Some quick code to generate a few plots for certain test cases

# x2 = np.linspace(-np.pi/3, np.pi/3)

# y3 = f(x2)

# y4 = curvature(x2)

# y5 = []

# for value in x2: 
# 	y5.append(integrate.quad(lambda z: curvature(z), -np.pi/3, value)[0])

# #print(x2)

# #print(y4)

# plt.plot(x2, y3, color = 'r')

# plt.show()

# plt.plot(x2, y4, color = 'g')

# plt.show()

# plt.plot(x2, y5, color = 'y')
# plt.show()

# for i in range(len(partition) - 1):
# 	x = np.linspace(partition[i], partition[i+1])
# 	y = []
# 	for value in x:
# 		y.append(integrate.quad(lambda z: curvature(z), -np.pi/3, value)[0])
# 	plt.plot(x, y)

# plt.show()

###############
#Square Example
###############

# def square_curvature_target_x_value(integration_lb, interval_lb, interval_ub, target_value, eps = 0.0001):
# 	x = np.linspace(interval_lb, interval_ub)
# 	i = 0
# 	for value in x:
# 		cumulated_curvature = integrate.quad(lambda z: 2.0/(1+4*z**2), integration_lb, value)[0]
# 		if np.abs(cumulated_curvature - target_value) <= eps:
# 			return value
# 		if target_value - cumulated_curvature > eps:
# 			i = i + 1
# 			continue
# 		if cumulated_curvature - target_value > eps:
# 			return square_curvature_target_x_value(integration_lb, x[i-1], x[i], target_value, eps)

# def square_eq_curvature_partition(lb, ub, Q = 20, eps = 0.0001):
# 	if Q == 1:
# 		return [lb, ub]
# 	breakpoints = [lb]
# 	target_value = (integrate.quad(lambda z: 2.0/(1+4*z**2), lb, ub)[0])/Q
# 	for i in range(Q-1):
# 		breakpoints.append(square_curvature_target_x_value(breakpoints[i], breakpoints[i], ub, target_value, eps))
# 	breakpoints.append(ub)
# 	return breakpoints

# partition = square_eq_curvature_partition(-5.0, 5.0, 100)

# x = np.linspace(-5.0, 5.0)

# y1 = x**2

# plt.plot(x, y1, color="g")
# plt.show()

# # y2 = 2.0/np.sqrt(1+4*x**2)
# # plt.plot(x, y2)
# # plt.show()

# # y3 = []

# # for value in x:
# # 	y3.append(integrate.quad(lambda z: 2.0/(1+4*z**2), -5.0, value)[0])

# # plt.plot(x, y3)
# # plt.show()

# for i in range(len(partition) - 1):
# 	x = np.linspace(partition[i], partition[i+1])
# 	y = x**2
# 	# for value in x:
# 	# 	y.append(integrate.quad(lambda z: 2.0/(1+4*z**2), -5.0, value)[0])
# 	plt.plot(x, y)

# plt.show()


			

