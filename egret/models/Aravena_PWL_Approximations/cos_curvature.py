#Code to produce values for cumulated curvature of the cosine function

import numpy as np

import scipy.integrate as integrate

import matplotlib.pyplot as plt


#############
#Functions for Curvature for Cosine
#############

def cos_curvature(delta):
	return np.absolute(-np.cos(delta))/(1+(np.sin(delta)**2))


#########################
#Functions for Partitioning the domain into Q pieces of equal curvature
#########################

def curvature_target_x_value(integration_lb, interval_lb, interval_ub, target_value, eps = 0.000001):
	#Finds an approximation for a value x in the domain where accumulated curvature is equal to a target value
	#The strategy is to recursively narrow down the given interval into smaller and smaller intervals until an acceptable point is found. 
	x = np.linspace(interval_lb, interval_ub)
	i = 0
	for value in x:
		cumulated_curvature = integrate.quad(cos_curvature, integration_lb, value)[0]
		if np.abs(cumulated_curvature - target_value) <= eps:
			return value
		if target_value - cumulated_curvature > eps:
			i = i + 1
			continue
		if cumulated_curvature - target_value > eps:
			return curvature_target_x_value(integration_lb, x[i-1], x[i], target_value, eps)

def eq_curvature_partition(lb, ub, Q, eps = 0.000001):
	#Divides the domain [lb, ub] into Q pieces of (approximately) equal curvature
	if Q == 1:
		return [lb, ub]
	breakpoints = [lb]
	target_value = (integrate.quad(cos_curvature, lb, ub)[0])/Q
	for i in range(Q-1):
		breakpoints.append(curvature_target_x_value(breakpoints[i], breakpoints[i], ub, target_value, eps))
	breakpoints.append(ub)
	return breakpoints

def eq_cosine_partition(lb, ub, Q):
	#Divides the domain [lb, ub] into Q pieces of equal curvature for the cosine function. Domain must be contained in (-pi/2, pi/2)
	if Q==1:
		return [lb, ub]
	total_curvature = np.arctan(np.sin(ub)) - np.arctan(np.sin(lb))
	breakpoints = [lb]
	for i in range(Q-1):
		breakpoints.append(np.arcsin(np.tan(total_curvature/Q + np.arctan(np.sin(breakpoints[i])))))
	breakpoints.append(ub)
	return breakpoints



if __name__ == '__main__':

    partition = eq_curvature_partition(-np.pi/6, np.pi/6, 16)

    cos_partition = eq_cosine_partition(-np.pi/6, np.pi/6, 16)

    print("Unrefined partition:", partition)

    print("Cosine partition:", cos_partition)

    x = np.linspace(-np.pi/6, np.pi/6)

    y1 = np.cos(x)

    plt.plot(x, y1)

    plt.show()

    for i in range(len(partition) - 1):
        x = np.linspace(partition[i], partition[i+1])
        y = np.cos(x)
        plt.plot(x, y)

    plt.show()

    for i in range(len(cos_partition) - 1):
    	x = np.linspace(cos_partition[i], cos_partition[i+1])
    	y = np.cos(x)
    	plt.plot(x, y)

    plt.show()