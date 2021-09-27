#Code to produce results from Aravena (2018) paper. 
import pdb
import numpy as np

import pyomo.environ as pe

import matplotlib.pyplot as plt

import matplotlib.patches as mpatches

from scipy.special.orthogonal import p_roots

from mpl_toolkits import mplot3d

import itertools as it

from Aravena_PWL_Approximations.Curvature import *

import json

#import pyomo.environ as pe
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
# import egret.model_library.transmission.bus as libbus
# import egret.model_library.transmission.branch as libbranch
# import egret.model_library.transmission.gen as libgen

# import egret.model_library.decl as decl
# from egret.model_library.defn import FlowType, CoordinateType, ApproximationType

# from egret.data.data_utils import map_items, zip_items
# from math import pi, radians



def power_flow_through_branch(Vi, Vj, delta, branch, bus_type = "from_bus", power_type = "Reactive"):
    if not (power_type =="Active" or power_type =="Reactive"):
        raise ValueError('Power type must be "Active" (for p) or "Reactive" (for q)')

    if not (bus_type == "from_bus" or bus_type == "to_bus"):
        raise ValueError('Bus type must be "from_bus" (for f) or "to_bus" (for t)')

    g = tx_calc.calculate_conductance(branch)
    b = tx_calc.calculate_susceptance(branch)

    if power_type == "Active":
        if bus_type == "from_bus":
        	return Vi**2*g - Vi*Vj*g*np.cos(delta) - Vi*Vj*b*np.sin(delta)
        else:
            return Vj**2*g - Vi*Vj*g*np.cos(delta) - Vi*Vj*b*np.sin(delta)

    else:
        if bus_type == "from_bus":
    	    return -Vi**2*b + Vi*Vj*b*np.cos(delta) - Vi*Vj*g*np.sin(delta)
        else:
            return -Vj**2*b + Vi*Vj*b*np.cos(delta) - Vi*Vj*g*np.sin(delta)



def triple_integral_coefs(Vi_lower, Vi_upper, Vj_lower, Vj_upper, delta_lower, delta_upper, branch, bus_type = "from_bus", power_type = "Active"):
    #Returns list of coefficients from the triple integral for a piecewise linear approximation of the power flow equations
    if not (power_type =="Active" or power_type =="Reactive"):
        raise ValueError('Power type must be "Active" (for p) or "Reactive" (for q)')

    if not (bus_type == "from_bus" or bus_type == "to_bus"):
        raise ValueError('Bus type must be "from_bus" (for f) or "to_bus" (for t)')
    
    g = tx_calc.calculate_conductance(branch)
    b = tx_calc.calculate_susceptance(branch)

    #For testing
    # g = 5
    # b = -15

    ######
    #Squared Terms
    ######
    
    #a_{i, 1}^2
    sq_term_1 = ((Vi_upper**3 - Vi_lower**3)/3)*(Vj_upper - Vj_lower)*(delta_upper - delta_lower)

    #a_{i, 2}^2
    sq_term_2 = (Vi_upper - Vi_lower)*((Vj_upper**3 - Vj_lower**3)/3)*(delta_upper - delta_lower)

    #a_{i, 3}^2
    sq_term_3 = (Vi_upper - Vi_lower)*(Vj_upper - Vj_lower)*((delta_upper**3 - delta_lower**3)/3)

    #b_i^2
    sq_term_4 = (Vi_upper - Vi_lower)*(Vj_upper - Vj_lower)*(delta_upper - delta_lower)

    #constant terms
    if power_type == "Active":
        if bus_type == "from_bus":
            sq_term_5 = g**2*((Vi_upper**5-Vi_lower**5)/5)*(Vj_upper - Vj_lower)*(delta_upper - delta_lower) + \
                        g**2*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**3-Vj_lower**3)/3)*((np.sin(2*delta_upper)/4 + delta_upper/2) - (np.sin(2*delta_lower)/4 + delta_lower/2)) + \
                        b**2*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**3-Vj_lower**3)/3)*((delta_upper/2 - np.sin(2*delta_upper)/4) - (delta_lower/2 - np.sin(2*delta_lower)/4))
        else:
            sq_term_5 = g**2*((Vj_upper**5-Vj_lower**5)/5)*(Vi_upper - Vi_lower)*(delta_upper - delta_lower) + \
                        g**2*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**3-Vj_lower**3)/3)*((np.sin(2*delta_upper)/4 + delta_upper/2) - (np.sin(2*delta_lower)/4 + delta_lower/2)) + \
                        b**2*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**3-Vj_lower**3)/3)*((delta_upper/2 - np.sin(2*delta_upper)/4) - (delta_lower/2 - np.sin(2*delta_lower)/4))


    else:
        if bus_type == "from_bus":
            sq_term_5 = b**2*((Vi_upper**5-Vi_lower**5)/5)*(Vj_upper - Vj_lower)*(delta_upper - delta_lower) + \
                        b**2*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**3-Vj_lower**3)/3)*((np.sin(2*delta_upper)/4 + delta_upper/2) - (np.sin(2*delta_lower)/4 + delta_lower/2)) + \
                        g**2*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**3-Vj_lower**3)/3)*((delta_upper/2 - np.sin(2*delta_upper)/4) - (delta_lower/2 - np.sin(2*delta_lower)/4))
        else:
            sq_term_5 = b**2*((Vj_upper**5-Vj_lower**5)/5)*(Vi_upper - Vi_lower)*(delta_upper - delta_lower) + \
                        b**2*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**3-Vj_lower**3)/3)*((np.sin(2*delta_upper)/4 + delta_upper/2) - (np.sin(2*delta_lower)/4 + delta_lower/2)) + \
                        g**2*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**3-Vj_lower**3)/3)*((delta_upper/2 - np.sin(2*delta_upper)/4) - (delta_lower/2 - np.sin(2*delta_lower)/4))
    #####
    #Cross Terms
    #####

    #a_{i,1}*a_{i,2}
    cross_term_1 = 2*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**2-Vj_lower**2)/2)*(delta_upper-delta_lower)

    #a_{i,1}*a_{i,3}
    cross_term_2 = 2*((Vi_upper**2-Vi_lower**2)/2)*(Vj_upper-Vj_lower)*((delta_upper**2-delta_lower**2)/2)

    #a_{i,1}*b_{i}
    cross_term_3 = 2*((Vi_upper**2-Vi_lower**2)/2)*(Vj_upper-Vj_lower)*(delta_upper-delta_lower)

    #a_{i,1}
    if power_type == "Active":
        if bus_type == "from_bus":
            cross_term_4 = -2*g*((Vi_upper**4-Vi_lower**4)/4)*(Vj_upper-Vj_lower)*(delta_upper-delta_lower)
        else:
            cross_term_4 = -2*g*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**3-Vj_lower**3)/3)*(delta_upper-delta_lower)
    else:
        if bus_type == "from_bus":
            cross_term_4 = 2*b*((Vi_upper**4-Vi_lower**4)/4)*(Vj_upper-Vj_lower)*(delta_upper-delta_lower)
        else:
            cross_term_4 = 2*b*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**3-Vj_lower**3)/3)*(delta_upper-delta_lower)

    #a_{i,1}
    if power_type == "Active":
        cross_term_5 = 2*g*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**2-Vj_lower**2)/2)*(np.sin(delta_upper) - np.sin(delta_lower))
    else:
        cross_term_5 = -2*b*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**2-Vj_lower**2)/2)*(np.sin(delta_upper) - np.sin(delta_lower))

    #a_{i,1}
    if power_type == "Active":
        cross_term_6 = 2*b*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**2-Vj_lower**2)/2)*(-np.cos(delta_upper) + np.cos(delta_lower))
    else:
        cross_term_6 = 2*g*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**2-Vj_lower**2)/2)*(-np.cos(delta_upper) + np.cos(delta_lower))
        

    #a_{i,2}*a_{i,3}
    cross_term_7 = 2*(Vi_upper-Vi_lower)*((Vj_upper**2-Vj_lower**2)/2)*((delta_upper**2-delta_lower**2)/2)

    #a_{i, 2}*b_i
    cross_term_8 = 2*(Vi_upper-Vi_lower)*((Vj_upper**2-Vj_lower**2)/2)*(delta_upper - delta_lower)

    #a_{i, 2}
    if power_type == "Active":
        if bus_type == "from_bus":
            cross_term_9 = -2*g*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**2-Vj_lower**2)/2)*(delta_upper - delta_lower)
        else:
            cross_term_9 = -2*g*(Vi_upper-Vi_lower)*((Vj_upper**4-Vj_lower**4)/4)*(delta_upper - delta_lower)
    else:
        if bus_type == "from_bus":
            cross_term_9 = 2*b*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**2-Vj_lower**2)/2)*(delta_upper - delta_lower)
        else:
            cross_term_9 = 2*b*(Vi_upper-Vi_lower)*((Vj_upper**4-Vj_lower**4)/4)*(delta_upper - delta_lower)

    #a_{i, 2}
    if power_type == "Active":
        cross_term_10 = 2*g*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**3-Vj_lower**3)/3)*(np.sin(delta_upper) - np.sin(delta_lower))
    else:
        cross_term_10 = -2*b*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**3-Vj_lower**3)/3)*(np.sin(delta_upper) - np.sin(delta_lower))

    #a_{i, 2}
    if power_type == "Active":
        cross_term_11 = 2*b*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**3-Vj_lower**3)/3)*(-np.cos(delta_upper) - (-np.cos(delta_lower)))
        
    else:
        cross_term_11 = 2*g*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**3-Vj_lower**3)/3)*(-np.cos(delta_upper) - (-np.cos(delta_lower)))


    #a_{i, 3}*b_i
    cross_term_12 = 2*(Vi_upper-Vi_lower)*(Vj_upper-Vj_lower)*((delta_upper**2-delta_lower**2)/2)

    #a_{i, 3}
    if power_type == "Active":
        if bus_type == "from_bus":
            cross_term_13 = -2*g*((Vi_upper**3-Vi_lower**3)/3)*(Vj_upper-Vj_lower)*((delta_upper**2-delta_lower**2)/2)
        else:
            cross_term_13 = -2*g*((Vj_upper**3-Vj_lower**3)/3)*(Vi_upper-Vi_lower)*((delta_upper**2-delta_lower**2)/2)
    else:
        if bus_type == "from_bus":
            cross_term_13 = 2*b*((Vi_upper**3-Vi_lower**3)/3)*(Vj_upper-Vj_lower)*((delta_upper**2-delta_lower**2)/2)
        else:
            cross_term_13 = 2*b*((Vj_upper**3-Vj_lower**3)/3)*(Vi_upper-Vi_lower)*((delta_upper**2-delta_lower**2)/2)

    #a_{i, 3}
    if power_type == "Active":
        cross_term_14 = 2*g*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**2-Vj_lower**2)/2)*((delta_upper*np.sin(delta_upper) + np.cos(delta_upper)) - (delta_lower*np.sin(delta_lower) + np.cos(delta_lower)))
    else:
        cross_term_14 = -2*b*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**2-Vj_lower**2)/2)*((delta_upper*np.sin(delta_upper) + np.cos(delta_upper)) - (delta_lower*np.sin(delta_lower) + np.cos(delta_lower)))

    #a_{i, 3}
    if power_type == "Active":
        cross_term_15 = 2*b*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**2-Vj_lower**2)/2)*((-delta_upper*np.cos(delta_upper) + np.sin(delta_upper)) - (-delta_lower*np.cos(delta_lower) + np.sin(delta_lower)))
    else:
        cross_term_15 = 2*g*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**2-Vj_lower**2)/2)*((-delta_upper*np.cos(delta_upper) + np.sin(delta_upper)) - (-delta_lower*np.cos(delta_lower) + np.sin(delta_lower)))


    #b_i
    if power_type == "Active":
        if bus_type == "from_bus":
            cross_term_16 = -2*g*((Vi_upper**3-Vi_lower**3)/3)*(Vj_upper - Vj_lower)*(delta_upper - delta_lower)
        else:
            cross_term_16 = -2*g*((Vj_upper**3-Vj_lower**3)/3)*(Vi_upper - Vi_lower)*(delta_upper - delta_lower)
    else:
        if bus_type == "from_bus":
            cross_term_16 = 2*b*((Vi_upper**3-Vi_lower**3)/3)*(Vj_upper - Vj_lower)*(delta_upper - delta_lower)
        else:
            cross_term_16 = 2*b*((Vj_upper**3-Vj_lower**3)/3)*(Vi_upper - Vi_lower)*(delta_upper - delta_lower) 

    #b_i
    if power_type == "Active":
        cross_term_17 = 2*g*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**2-Vj_lower**2)/2)*(np.sin(delta_upper)-np.sin(delta_lower))
    else:
        cross_term_17 = -2*b*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**2-Vj_lower**2)/2)*(np.sin(delta_upper)-np.sin(delta_lower))

    #b_i
    if power_type == "Active":
        cross_term_18 = 2*b*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**2-Vj_lower**2)/2)*(-np.cos(delta_upper) + np.cos(delta_lower))
    else:
        cross_term_18 = 2*g*((Vi_upper**2-Vi_lower**2)/2)*((Vj_upper**2-Vj_lower**2)/2)*(-np.cos(delta_upper) + np.cos(delta_lower))


    #constant
    if power_type == "Active":
        if bus_type == "from_bus":
            cross_term_19 = -2*g**2*((Vi_upper**4-Vi_lower**4)/4)*((Vj_upper**2-Vj_lower**2)/2)*(np.sin(delta_upper) - np.sin(delta_lower))
        else:
            cross_term_19 = -2*g**2*((Vj_upper**4-Vj_lower**4)/4)*((Vi_upper**2-Vi_lower**2)/2)*(np.sin(delta_upper) - np.sin(delta_lower))
    else:
        if bus_type == "from_bus":
            cross_term_19 = -2*b**2*((Vi_upper**4-Vi_lower**4)/4)*((Vj_upper**2-Vj_lower**2)/2)*(np.sin(delta_upper) - np.sin(delta_lower))
        else:
            cross_term_19 = -2*b**2*((Vj_upper**4-Vj_lower**4)/4)*((Vi_upper**2-Vi_lower**2)/2)*(np.sin(delta_upper) - np.sin(delta_lower))

    #constant
    if power_type == "Active":
        if bus_type == "from_bus":
            cross_term_20 = -2*g*b*((Vi_upper**4-Vi_lower**4)/4)*((Vj_upper**2-Vj_lower**2)/2)*(-np.cos(delta_upper) - (-np.cos(delta_lower)))
        else:
            cross_term_20 = -2*g*b*((Vj_upper**4-Vj_lower**4)/4)*((Vi_upper**2-Vi_lower**2)/2)*(-np.cos(delta_upper) - (-np.cos(delta_lower)))
    else:
        if bus_type == "from_bus":
            cross_term_20 = 2*g*b*((Vi_upper**4-Vi_lower**4)/4)*((Vj_upper**2-Vj_lower**2)/2)*(-np.cos(delta_upper) - (-np.cos(delta_lower)))
        else:
            cross_term_20 = 2*g*b*((Vj_upper**4-Vj_lower**4)/4)*((Vi_upper**2-Vi_lower**2)/2)*(-np.cos(delta_upper) - (-np.cos(delta_lower)))

    #constant
    if power_type == "Active":
        cross_term_21 = g*b*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**3-Vj_lower**3)/3)*(-np.cos(2*delta_upper)/2 - (-np.cos(2*delta_lower)/2))
    else:
        cross_term_21 = -g*b*((Vi_upper**3-Vi_lower**3)/3)*((Vj_upper**3-Vj_lower**3)/3)*(-np.cos(2*delta_upper)/2 - (-np.cos(2*delta_lower)/2))


    #Return list of coefficients [a_{i,1}^2, a_{i, 2}^2, a_{i, 3}^2, b_i^2, a_{i,1}a_{i,2}, a_{i,1}a_{i,3}, a_{i,1}b_i, a_{i,2}a_{i,3}, a_{i,2}b_{i}, a_{i,3}b_i, a_{i,1}, a_{i,2}, a_{i,3}, b_i, constant]
    #Could also return as dict instead? -WL

    ai1_sq_coef = sq_term_1
    ai2_sq_coef = sq_term_2
    ai3_sq_coef = sq_term_3
    bi_sq_coef = sq_term_4
    ai1ai2_coef = cross_term_1
    ai1ai3_coef = cross_term_2
    ai1bi_coef = cross_term_3
    ai2ai3_coef = cross_term_7
    ai2bi_coef = cross_term_8
    ai3bi_coef = cross_term_12
    ai1_coef = cross_term_4 + cross_term_5 + cross_term_6
    ai2_coef = cross_term_9 + cross_term_10 + cross_term_11
    ai3_coef = cross_term_13 + cross_term_14 + cross_term_15
    bi_coef = cross_term_16 + cross_term_17 + cross_term_18
    constant_coef = sq_term_5 + cross_term_19 + cross_term_20 + cross_term_21

    return [ai1_sq_coef, ai2_sq_coef, ai3_sq_coef, bi_sq_coef, ai1ai2_coef, ai1ai3_coef, ai1bi_coef, ai2ai3_coef, ai2bi_coef, ai3bi_coef, ai1_coef, ai2_coef, ai3_coef, bi_coef, constant_coef]


def generate_pwl_model(Vi_lower, Vi_upper, Vj_lower, Vj_upper, delta_lower, delta_upper, branch, num_delta_boxes, bus_type = "from_bus", power_type = "Active"):
    #Given a box, solve for the piecewise linear approximation over this box.
    #I'd like to give the user functions to find partitions of the domains of Vi, Vj and Delta,
    #but for now I will specify them manually.

    Vi_partition = [Vi_lower, Vi_upper]
    Vj_partition = [Vj_lower, Vj_upper]
    delta_partition = eq_curvature_partition(delta_lower, delta_upper, num_delta_boxes, branch, bus_type = bus_type, power_type=power_type)
    model = pe.ConcreteModel()

    #Create a list of boxes from the elements of the given partition. They will be specified by a triple of coordinates
    Vi_intervals = [(Vi_partition[i], Vi_partition[i+1]) for i in range(0, len(Vi_partition) - 1)]
    Vj_intervals = [(Vj_partition[i], Vj_partition[i+1]) for i in range(0, len(Vj_partition) - 1)]
    delta_intervals = [(delta_partition[i], delta_partition[i+1]) for i in range(0, len(delta_partition) - 1)]
    boxes = list(it.product(Vi_intervals, Vj_intervals, delta_intervals))
    #print(boxes)
    model.boxes = pe.Set(initialize=range(len(boxes)))

    facets = []

    for i in it.product(Vi_intervals, Vj_intervals):
        for j in range(len(delta_intervals) - 1):
            facets.append(('delta', (i[0], i[1], delta_intervals[j]), (i[0], i[1], delta_intervals[j+1])))

    for i in it.product(Vj_intervals, delta_intervals):
        for j in range(len(Vi_intervals) - 1):
            facets.append(('Vi', (i[0], i[1], Vi_intervals[j]), (i[0], i[1], Vi_intervals[j+1])))

    for i in it.product(Vi_intervals, delta_intervals):
        for j in range(len(Vj_intervals) - 1):
            facets.append(('Vj', (i[0], i[1], Vj_intervals[j]), (i[0], i[1], Vj_intervals[j+1])))
    # print(facets)
    # print(len(facets))

    model.facets = pe.Set(initialize=range(len(facets)))
    model.N = pe.Set(initialize=[0,1,2,3])

    #Create variables ai1, ai2, ai3, bi for each i in boxes (these are the coefficients for Vi, Vj, delta_ij and the constant term respectively)
    model.ai1 = pe.Var(model.boxes)
    model.ai2 = pe.Var(model.boxes)
    model.ai3 = pe.Var(model.boxes)
    model.bi = pe.Var(model.boxes)

    #Objective is to minimize sum over all boxes with coefs derived
    
    def box_obj_rule(model, i):
        box = boxes[i]
        box_coefs = triple_integral_coefs(box[0][0], box[0][1], box[1][0], box[1][1], box[2][0], box[2][1], branch, bus_type=bus_type, power_type=power_type)
        return box_coefs[0]*model.ai1[i]*model.ai1[i] + box_coefs[1]*model.ai2[i]*model.ai2[i] + box_coefs[2]*model.ai3[i]*model.ai3[i] + box_coefs[3]*model.bi[i]*model.bi[i] + \
                        box_coefs[4]*model.ai1[i]*model.ai2[i] +box_coefs[5]*model.ai1[i]*model.ai3[i] + box_coefs[6]*model.ai1[i]*model.bi[i] + box_coefs[7]*model.ai2[i]*model.ai3[i] + \
                        box_coefs[8]*model.ai2[i]*model.bi[i] + box_coefs[9]*model.ai3[i]*model.bi[i] + box_coefs[10]*model.ai1[i] + box_coefs[11]*model.ai2[i] + box_coefs[12]*model.ai3[i] + \
                        box_coefs[13]*model.bi[i] + box_coefs[14]

    def obj_rule(model):
        expr = 0
        for i in range(len(boxes)):
            expr += box_obj_rule(model, i)
        return expr

    model.obj = pe.Objective(rule=obj_rule)

    #Constraints to enforce continuity on each facet 

    def facet_cons_rule(model, i, vertex):
        facet = facets[i]
        if facet[0] == 'delta':
            l_box = boxes.index(facet[1])
            r_box = boxes.index(facet[2])
            if vertex == 0:
                return model.ai1[l_box]*facet[1][0][0] + model.ai2[l_box]*facet[1][1][0] + model.ai3[l_box]*facet[1][2][1] + model.bi[l_box] == \
                        model.ai1[r_box]*facet[2][0][0] + model.ai2[r_box]*facet[2][1][0] + model.ai3[r_box]*facet[2][2][0] + model.bi[r_box]
            if vertex == 1:
                return model.ai1[l_box]*facet[1][0][0] + model.ai2[l_box]*facet[1][1][1] + model.ai3[l_box]*facet[1][2][1] + model.bi[l_box] == \
                        model.ai1[r_box]*facet[2][0][0] + model.ai2[r_box]*facet[2][1][1] + model.ai3[r_box]*facet[2][2][0] + model.bi[r_box]
            if vertex == 2:
                return model.ai1[l_box]*facet[1][0][1] + model.ai2[l_box]*facet[1][1][0] + model.ai3[l_box]*facet[1][2][1] + model.bi[l_box] == \
                        model.ai1[r_box]*facet[2][0][1] + model.ai2[r_box]*facet[2][1][0] + model.ai3[r_box]*facet[2][2][0] + model.bi[r_box]
            if vertex == 3:
                return model.ai1[l_box]*facet[1][0][1] + model.ai2[l_box]*facet[1][1][1] + model.ai3[l_box]*facet[1][2][1] + model.bi[l_box] == \
                        model.ai1[r_box]*facet[2][0][1] + model.ai2[r_box]*facet[2][1][1] + model.ai3[r_box]*facet[2][2][0] + model.bi[r_box]
            return pe.Constraint.Skip

        if facet[0] == 'Vi':
            pdb.set_trace()
            l_box = boxes.index((facet[1][2], facet[1][0], facet[1][1]))
            r_box = boxes.index((facet[2][2], facet[2][0], facet[2][1]))
            if vertex == 0:
                return model.ai2[l_box]*facet[1][0][0] + model.ai3[l_box]*facet[1][1][0] + model.ai1[l_box]*facet[1][2][1] + model.bi[l_box] == \
                        model.ai2[r_box]*facet[1][0][0] + model.ai3[r_box]*facet[1][1][0] + model.ai1[r_box]*facet[1][2][1] + model.bi[r_box]
            if vertex == 1:
                return model.ai2[l_box]*facet[1][0][0] + model.ai3[l_box]*facet[1][1][1] + model.ai1[l_box]*facet[1][2][1] + model.bi[l_box] == \
                        model.ai2[r_box]*facet[1][0][0] + model.ai3[r_box]*facet[1][1][1] + model.ai1[r_box]*facet[1][2][1] + model.bi[r_box]
            if vertex == 2:
                return model.ai2[l_box]*facet[1][0][1] + model.ai3[l_box]*facet[1][1][0] + model.ai1[l_box]*facet[1][2][1] + model.bi[l_box] == \
                        model.ai2[r_box]*facet[1][0][1] + model.ai3[r_box]*facet[1][1][0] + model.ai1[r_box]*facet[1][2][1] + model.bi[r_box]
            #if vertex == 3:
            #    return model.ai2[l_box]*facet[1][0][1] + model.ai3[l_box]*facet[1][1][1] + model.ai1[l_box]*facet[1][2][1] + model.bi[l_box] == \
            #            model.ai2[r_box]*facet[1][0][1] + model.ai3[r_box]*facet[1][1][1] + model.ai1[r_box]*facet[1][2][1] + model.bi[r_box]

        if facet[0] == 'Vj':
            pdb.set_trace()
            l_box = boxes.index((facet[1][0], facet[1][2], facet[1][1]))
            r_box = boxes.index((facet[2][0], facet[2][2], facet[2][1]))
            if vertex == 0:
                return model.ai1[l_box]*facet[1][0][0] + model.ai3[l_box]*facet[1][1][0] + model.ai2[l_box]*facet[1][2][1] + model.bi[l_box] == \
                        model.ai1[r_box]*facet[1][0][0] + model.ai3[r_box]*facet[1][1][0] + model.ai2[r_box]*facet[1][2][1] + model.bi[r_box]
            if vertex == 1:
                return model.ai1[l_box]*facet[1][0][0] + model.ai3[l_box]*facet[1][1][1] + model.ai2[l_box]*facet[1][2][1] + model.bi[l_box] == \
                        model.ai1[r_box]*facet[1][0][0] + model.ai3[r_box]*facet[1][1][1] + model.ai2[r_box]*facet[1][2][1] + model.bi[r_box]
            if vertex == 2:
                return model.ai1[l_box]*facet[1][0][1] + model.ai3[l_box]*facet[1][1][0] + model.ai2[l_box]*facet[1][2][1] + model.bi[l_box] == \
                        model.ai1[r_box]*facet[1][0][1] + model.ai3[r_box]*facet[1][1][0] + model.ai2[r_box]*facet[1][2][1] + model.bi[r_box]
            #if vertex == 3:
            #    return model.ai1[l_box]*facet[1][0][1] + model.ai3[l_box]*facet[1][1][1] + model.ai2[l_box]*facet[1][2][1] + model.bi[l_box] == \
            #            model.ai1[r_box]*facet[1][0][1] + model.ai3[r_box]*facet[1][1][1] + model.ai2[r_box]*facet[1][2][1] + model.bi[r_box]


    #model.facet_Constr = pe.Constraint(model.facets, model.N, rule=facet_cons_rule)


    return [model, boxes]

def pwl_model_branch_preprocessing(Vi_lower, Vi_upper, Vj_lower, Vj_upper, delta_lower, delta_upper, branch, num_delta_boxes, bus_type = "from_bus", power_type = "Active"):
    #Obtains information for a given branch based on the branch parameters
    m = generate_pwl_model(Vi_lower, Vi_upper, Vj_lower, Vj_upper, delta_lower, delta_upper, branch, num_delta_boxes = num_delta_boxes, bus_type=bus_type, power_type = power_type)
    pwl_model = m[0]
    boxes = m[1]

    branch_dict={}

    box_coords_list = []
    for box in boxes:
        box_coords = [(box[0][0], box[1][0], box[2][0]), (box[0][0], box[1][0], box[2][1]), (box[0][0], box[1][1], box[2][0]), (box[0][0], box[1][1], box[2][1]), \
                        (box[0][1], box[1][0], box[2][0]), (box[0][1], box[1][0], box[2][1]), (box[0][1], box[1][1], box[2][0]), (box[0][1], box[1][1], box[2][1])]
        box_coords_list.append(box_coords)

    branch_dict['from_bus'] = branch['from_bus']
    branch_dict['to_bus'] = branch['to_bus']

    opt = pe.SolverFactory("baron")

    #Knitro_options 

    opt.solve(pwl_model, tee=True)

    coeffs = []

    for i in range(len(boxes)):
        box_coeffs = [round(pe.value(pwl_model.ai1[i]), 5), round(pe.value(pwl_model.ai2[i]), 5), round(pe.value(pwl_model.ai3[i]), 5), round(pe.value(pwl_model.bi[i]), 5)]
        coeffs.append(box_coeffs)

    branch_dict['boxes'] = {'coords': box_coords_list, 'coefficients': coeffs}

    return branch_dict

if __name__ == '__main__':
    import os
    from egret.parsers.matpower_parser import create_ModelData

    path = os.path.dirname(os.path.dirname(os.getcwd()))
    power_types = ["Reactive", "Active"]
    bus_types = ["from_bus", "to_bus"]
    case = 'case14_ieee'
    filename = 'pglib_opf_' + case + '.m'
    test_case = os.path.join(path, 'thirdparty', 'pglib-opf-master', filename)    
    md_dict = create_ModelData(test_case)
    md = md_dict.clone_in_service()
    tx_utils.scale_ModelData_to_pu(md, inplace = True)

    branches = dict(md.elements(element_type='branch'))
    branch_attrs = md.attributes(element_type='branch')

    buses = dict(md.elements(element_type='bus'))
    bus_attrs = md.attributes(element_type='bus')

    num_delta_boxes = 4
    #s_max = branch['rating_long_term']

    #print(triple_integral_coefs(0.95, 1.05, 0.95, 1.05, 0, np.pi/6, branch, power_type = "Active"))

    #m = generate_pwl_model(0.95, 1.05, 0.95, 1.05, 0, np.pi/6, branch, power_type = "Reactive")

    all_branches_coefs_dict = dict([])
    for bus_type in bus_types:
            for power_type in power_types:
                branch_p_b_dict = dict([])
                for branch_name in branches.keys():
                    branch = branches[branch_name]
                    from_bus = branch['from_bus']
                    to_bus = branch['to_bus']
                    if bus_attrs['v_min'][from_bus] == None:
                        from_bus_lb = 0.95
                    else:
                        from_bus_lb = bus_attrs['v_min'][from_bus]
                    if bus_attrs['v_max'][from_bus] == None:
                        from_bus_ub = 1.05
                    else:
                        from_bus_ub = bus_attrs['v_max'][from_bus]
                    if bus_attrs['v_min'][to_bus] == None:
                        to_bus_lb = 0.95
                    else:
                        to_bus_lb = bus_attrs['v_min'][to_bus]    
                    if bus_attrs['v_max'][to_bus] == None:
                        to_bus_ub = 1.05
                    else:
                        to_bus_ub = bus_attrs['v_max'][to_bus]
                    if branch['angle_diff_min'] == None:
                        delta_lb = (-np.pi)/6
                    else:
                        delta_lb = branch['angle_diff_min']*(np.pi)/180.0
                    if branch['angle_diff_max'] == None:
                        delta_ub = (np.pi)/6
                    else:
                        delta_ub = branch['angle_diff_max']*(np.pi)/180.0
                    branch_p_b_dict[branch_name] = pwl_model_branch_preprocessing(from_bus_lb, from_bus_ub, to_bus_lb, to_bus_ub, delta_lb, delta_ub, branch, num_delta_boxes, bus_type, power_type)
                all_branches_coefs_dict[power_type + '_' + bus_type] = branch_p_b_dict
    
    json_filename = case + '_delta_' + str(num_delta_boxes) + '_curvature_partition.json'
    with open(json_filename, 'w') as file:
        json.dump(all_branches_coefs_dict, file, indent=2)

    

    #opt = pe.SolverFactory("knitroampl")

    #Knitro_options 

    #opt.solve(m, tee=True)

    # print(branch)
    # print(from_bus_lb)
    # print(from_bus_ub)
    # print(to_bus_lb)
    # print(to_bus_ub)
    # print(s_max)

    # for v in m.component_objects(pe.Var):
    #     for index in v:
    #         print('{0} = {1}'.format(v[index], pe.value(v[index])))



    # x = np.linspace(from_bus_lb, from_bus_ub)
    # y = np.linspace(to_bus_lb, to_bus_ub)

    # delta = np.linspace(-np.pi/6, np.pi/6)

    # X, DELTA = np.meshgrid(x, delta)

    # z = power_flow_through_branch(X, 1, DELTA, branch)

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')

    # ax.plot_surface(X, DELTA, z, cmap='viridis', edgecolor='none')
    # ax.set_title('Power function with constant voltage')
    # plt.show()

    

    # print((power_flow_through_branch(X, 1, DELTA, branch, power_type="Active")**2 + power_flow_through_branch(X, 1, DELTA, branch, power_type="Reactive")**2 <= s_max**2).astype(int))



    # plt.imshow( (power_flow_through_branch(X, 1, DELTA, branch, power_type="Active")**2 + power_flow_through_branch(X, 1, DELTA, branch, power_type="Reactive")**2 <= s_max**2).astype(int),  # \
				# #(power_flow_through_branch(Y, X, -delta, branch, power_type="Active")**2 + power_flow_through_branch(Y, X, -delta, branch, power_type="Reactive")**2 <= s_max**2)).astype(int), \
				# extent=(x.min(),x.max(),delta.min(),delta.max()),origin="lower", cmap="Blues")

    # plt.xlim(from_bus_lb, from_bus_ub)
    # plt.ylim(to_bus_lb, to_bus_ub)

    # plt.show()
















