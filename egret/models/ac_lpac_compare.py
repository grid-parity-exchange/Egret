"""
This module provides means of comparison for LPAC and ACOPF formulations

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

from egret.models.lpac import create_cold_start_lpac_model, solve_lpac


if __name__ == '__main__':
    import os
    from egret.parsers.matpower_parser import create_ModelData

    filename = 'pglib_opf_case30_ieee.m'
    test_case = os.path.join('c:\\', 'Users', 'wlinz', 'Desktop', 'Restoration', 'Egret', 'egret', 'thirdparty', 'pglib-opf-master', filename) #Better if this isn't so user-dependent
    model_data = create_ModelData(test_case)
    baron_options = {'LPSol': 8}# 'CplexLibName': "C:/Program Files/IBM/ILOG/CPLEX_Studio129/cplex/bin/x64_win64/cplex1290.dll"}
    knitro_options = {'maxit': 20000}
    kwargs = {'include_feasibility_slack':False}
    kwargs_for_lpac={}
    #md,m,results = solve_lpac(model_data, "baron",options=baron_options,lpac_model_generator=create_cold_start_lpac_model,return_model=True, return_results=True,**kwargs)
    #md,m,results = solve_lpac(model_data, "knitroampl",options=knitro_options,lpac_model_generator=create_hot_start_lpac_model,return_model=True, return_results=True,**kwargs)
    #md,m,results = solve_lpac(model_data, "knitroampl",options=knitro_options,lpac_model_generator=create_warm_start_lpac_model,return_model=True, return_results=True,**kwargs)
    md,m,results = solve_lpac(model_data, "cplex",lpac_model_generator=create_cold_start_lpac_model,return_model=True, return_results=True, kwargs=kwargs, kwargs_for_lpac=kwargs_for_lpac)
    updated_knitro_options = {'maxit': 20000, 'strat_warm_start': 1}
    ac_md, ac_m, ac_results = solve_acopf(model_data, "knitroampl",options=knitro_options,acopf_model_generator=create_rsv_acopf_model,return_model=True, return_results=True,**kwargs)
    ac_md, ac_m, ac_results = solve_acopf(model_data, "baron",options=baron_options,acopf_model_generator=create_rsv_acopf_model,return_model=True, return_results=True,**kwargs)
    
    #print(ac_results
    lpac_pg_gens = []
    lpac_qg_gens = []
    for gen in md.elements(element_type="generator"):
    	lpac_pg_gens.append(gen[1]['pg'])
    	lpac_qg_gens.append(gen[1]['qg'])
    
    acopf_pg_gens = []
    acopf_qg_gens = []
    for gen in ac_md.elements(element_type="generator"):
    	acopf_pg_gens.append(gen[1]['pg'])
    	acopf_qg_gens.append(gen[1]['qg'])
    
    # for bus in ac_md.elements(element_type="bus"):
    # 	print(bus[0])
    # 	print(bus[1])


    print("LPAC objective: ", pe.value(m.obj()))
    print("AC objective: ", pe.value(ac_m.obj()))
    print("\n")
    print("Real Power\n")
    pg_differences = [lpac_pg_gens[i] - acopf_pg_gens[i] for i in range(len(lpac_pg_gens))]
    #print("Actual differences between lpac and ac: ", pg_differences)
    abs_pg_differences = [abs(pg_differences[i]) for i in range(len(lpac_pg_gens))]
    #print("Absolute value of differences: ", abs_pg_differences)
    print("L_1 norm: ", sum(abs_pg_differences))
    average_diff = sum(abs_pg_differences)/(len(abs_pg_differences))
    print("Average difference of absolute values: ", average_diff)
    max_abs_diff = max(abs_pg_differences)
    print("L_infty norm: ", max_abs_diff)
    square_diff = sqrt(sum(abs_pg_differences[i]**2 for i in range(len(lpac_pg_gens))))
    print("L_2 norm: ", square_diff)

    print("\n")

    print("Reactive Power\n")
    qg_differences = [lpac_qg_gens[i] - acopf_qg_gens[i] for i in range(len(lpac_qg_gens))]
    #print("Actual differences between lpac and ac: ", qg_differences)
    abs_qg_differences = [abs(qg_differences[i]) for i in range(len(lpac_qg_gens))]
    #print("Absolute value of differences: ", abs_qg_differences)
    print("L_1 norm: ", sum(abs_qg_differences))
    average_diff = sum(abs_qg_differences)/(len(abs_qg_differences))
    print("Average difference of absolute values: ", average_diff)
    max_abs_diff = max(abs_qg_differences)
    print("L_infty norm: ", max_abs_diff)
    square_diff = sqrt(sum(abs_qg_differences[i]**2 for i in range(len(lpac_qg_gens))))
    print("L_2 norm: ", square_diff)