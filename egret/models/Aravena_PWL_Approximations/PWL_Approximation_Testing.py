#Testing functions for the PWL approximation -WL

import numpy as np 

import itertools as it

from PWL_Approx_Functions import *

import json





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


	
	json_filename = case + '_delta_10_curvature_partition.json'
	with open(json_filename, "r") as read_file:
		full_branches_dict = json.load(read_file)

	branch_av_sq_diff_dict = dict([])

	for power_type in power_types:
		for bus_type in bus_types:
			p_b_branches_dict = full_branches_dict[power_type + '_' + bus_type]
			p_b_branch_av_sq_diff_dict = dict([])
			for branch_name in p_b_branches_dict.keys():
				branch_sq_diffs = []
				for i in range(len(p_b_branches_dict[branch_name]['boxes']['coords'])):
					box = p_b_branches_dict[branch_name]['boxes']['coords'][i]
					Vi = np.linspace(box[0][0], box[7][0], num=5)
					Vj = np.linspace(box[0][1], box[7][1], num=5)
					delta = np.linspace(box[0][2], box[7][2], num=5)
					box_sq_diffs = []
					for point in it.product(Vi, Vj, delta):
						branch = branches[branch_name]
						power_flow = power_flow_through_branch(point[0], point[1], point[2], branch, bus_type = bus_type, power_type=power_type) #Get actual power flow value
						coeffs_list = p_b_branches_dict[branch_name]['boxes']['coefficients'][i]
						pwl_approx = coeffs_list[0]*point[0]+coeffs_list[1]*point[1]+coeffs_list[2]*point[2]+coeffs_list[3]
						box_sq_diffs.append((power_flow - pwl_approx)**2)
					branch_sq_diffs.append((sum(box_sq_diffs))/(len(box_sq_diffs)))
				p_b_branch_av_sq_diff_dict[branch_name] = (sum(branch_sq_diffs))/(len(branch_sq_diffs))
			branch_av_sq_diff_dict[power_type + '_' + bus_type] = p_b_branch_av_sq_diff_dict

	print(branch_av_sq_diff_dict) 


