from egret.data.model_data import ModelData
 
import egret.model_library.transmission.tx_calc as tx_calc
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg

import sys
import time

from pyutilib.misc.timing import TicTocTimer

tt = TicTocTimer()

error = False

md = ModelData.read(sys.argv[1])

branches = dict(md.elements(element_type='branch')); buses = dict(md.elements(element_type='bus')); index_set_bus = tuple(buses.keys()); index_set_branch = tuple(branches.keys()); reference_bus = md.data['system']['reference_bus']; mapping_bus_to_idx = { b : idx for idx,b in enumerate(index_set_bus)}; mapping_branch_to_idx = {b : idx for idx,b in enumerate(index_set_branch) };

graph = tx_calc.construct_connection_graph(branches, mapping_bus_to_idx)
tt.tic("Using networkx for N-1 connection check")
branches_not_disconnecting = tx_calc.get_N_minus_1_branches(graph, branches, mapping_bus_to_idx)
num_connected = len(branches_not_disconnecting)
print("Number of N-1 contingencies monitorable", num_connected)
tt.toc("Finished N-1 connection check")

contingencies = dict()
for bn in branches_not_disconnecting:
    contingencies[f'cont_{bn}'] = {'branch_contingency': bn}

MLU_MP, B_dA, ref_bus_mask, contingency_compensators = tx_calc.calculate_ptdf_factorization(branches, buses, index_set_branch, index_set_bus, reference_bus, mapping_bus_to_idx=mapping_bus_to_idx, mapping_branch_to_idx=mapping_branch_to_idx, contingencies=contingencies)
tt.toc("Finished pre-computations")
#time.sleep(60)
#print(contingency_compensators.keys())
