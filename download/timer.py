import sys
import gc

import numpy as np

from pyutilib.misc.timing import TicTocTimer

fn = sys.argv[1]
tt_timer = TicTocTimer()

#if sys.argv[2] == 'factor':
#    use_factorization = True
#elif sys.argv[2] == 'ptdf':
#    use_factorization = False
#else:
#    sys.exit(1)
use_factorization = False

from egret.data.model_data import ModelData
from egret.model_library.transmission.tx_calc import calculate_ptdf, _calculate_J11, calculate_adjacency_matrix_transpose, check_network_connection
from egret.data.ptdf_utils import PTDFMatrix
from egret.model_library.defn import BasePointType, ApproximationType
from egret.common.lazy_ptdf_utils import populate_default_ptdf_options

tt_timer.tic(f'reading file {fn}')
md = ModelData.read(fn)

tt_timer.tic(f'setting up dictionaries')
buses = dict(md.elements(element_type = 'bus'))
branches = dict(md.elements(element_type = 'branch'))

reference_bus = md.data['system']['reference_bus']

index_set_bus = tuple(buses.keys())
index_set_branch = tuple(branches.keys())

mapping_bus_to_idx = { bus_n : i for i, bus_n in enumerate(index_set_bus) }

tt_timer.tic(f'skeleton pyomo model')
import pyomo.environ as pe
mb = pe.ConcreteModel()
_p_nw = { bus : 0. for bus in index_set_bus }

for g, g_dict in md.elements(element_type='generator'):
    _p_nw[g_dict['bus']] -= g_dict['pg']

for l, l_dict in md.elements(element_type='load'):
    _p_nw[l_dict['bus']] += l_dict['p_load']

mb.p_nw = pe.Param(index_set_bus, mutable=True, initialize=_p_nw)

tt = tt_timer
#print(md.data['elements']['bus'])
#print(md.data['elements']['branch'])
tt.tic(f"computing full PTDF")
PTDF = PTDFMatrix(branches, buses, reference_bus, BasePointType.FLATSTART, populate_default_ptdf_options(None))
tt.toc(f"COMPUTED full PTDF")

tt.tic(f"computing flows with full PTDF")
PFV, _ = PTDF.calculate_all_flows(mb)
tt.toc(f"COMPUTED flows with full PTDF")
#print(f"PTDF: {PTDF.PTDFM_masked}")
print(f"masked: {PTDF.masked}")
#print(f"PFV: {PFV}")
del PTDF
gc.collect()

## try masked
md.data['elements']['bus']['1']['base_kv'] = 1000.

ptdf_options = populate_default_ptdf_options(None)
ptdf_options['branch_kv_threshold'] = 999.
ptdf_options['kv_threshold_type'] = 'one'
tt.tic(f"computing partial PTDF")
PTDF = PTDFMatrix(branches, buses, reference_bus, BasePointType.FLATSTART, ptdf_options)
tt.toc(f"COMPUTED partial PTDF")

tt.tic(f'computing monitored flows')
PTDF.calculate_monitored_flows(mb)
tt.toc(f'COMPUTED monitored flows')

tt.tic(f"computing flows with lu_solve")
PFV_BS, _ = PTDF.calculate_all_flows(mb)
tt.toc(f"COMPUTED flows with lu_solve")
#print(f"PTDF: {PTDF.PTDFM_masked}")
print(f"masked: {PTDF.masked}")
#print(f"PFV: {PFV_BS}")

if np.allclose(PFV, PFV_BS):
    print(f"Both flows the same")
else:
    print(f"BAAAHHHHHH some issue with flows")

'''
def get_ptdf(branches,buses,index_set_branch,index_set_bus,reference_bus,mapping_bus_to_idx):
    #PTDF = calculate_ptdf(branches,buses,index_set_branch,index_set_bus,reference_bus,mapping_bus_to_idx=mapping_bus_to_idx, sparse_index_set_branch=[index_set_branch[0], index_set_branch[1]])
    PTDF = calculate_ptdf(branches,buses,index_set_branch,index_set_bus,reference_bus,mapping_bus_to_idx=mapping_bus_to_idx)
    return PTDF

if False:
    import numpy as np
    base_point = BasePointType.FLATSTART
    _len_bus = len(index_set_bus)

    if mapping_bus_to_idx is None:
        mapping_bus_to_idx = {bus_n: i for i, bus_n in enumerate(index_set_bus)}

    _len_branch = len(index_set_branch)

    _ref_bus_idx = mapping_bus_to_idx[reference_bus]

    ## check if the network is connected
    tt.tic("checking connection")
    connected = check_network_connection(branches, index_set_branch, index_set_bus, mapping_bus_to_idx)

    tt.tic("getting J matrix")
    J = _calculate_J11(branches,buses,index_set_branch,index_set_bus,mapping_bus_to_idx,base_point,approximation_type=ApproximationType.PTDF)
    tt.toc("GOT J matrix")
    tt.tic("getting A matrix")
    A = calculate_adjacency_matrix_transpose(branches,index_set_branch,index_set_bus,mapping_bus_to_idx)
    tt.toc("GOT A matrix")

    tt.tic("computing M matrix")
    M = A@J
    tt.toc("COMPUTED M matrix")
    if connected:

        tt.tic("getting ref_bus_mask")
        ref_bus_mask = np.ones(_len_bus, dtype=bool)
        ref_bus_mask[_ref_bus_idx] = False
        tt.toc("GOT ref_bus_mask")

        # M is now (A^T B_d A) with
        # row and column of reference
        # bus removed
        tt.tic("getting J0")
        J0 = M[ref_bus_mask,:][:,ref_bus_mask]
        tt.toc("GOT J0")

        # (B_d A) with reference bus column removed
        tt.tic("getting B_dA")
        B_dA = J[:,ref_bus_mask].A
        tt.toc("GOT B_dA")

        tt.tic("transposing/densifying J0")
        J0_T_A = J0.T.A
        tt.toc("TRANSPOSED/DENSIFYED J0") 

        tt.tic("transposing B_dA")
        B_dA_T = B_dA.T
        tt.toc("TRANSPOSED B_dA") 

        #tt.tic("computing PTDF")
        #PTDF = np.linalg.solve(J0_T_A, B_dA_T)
        #tt.toc("COMPUTED PTDF")

        #tt.tic("transposing PTDF")
        #PTDF = PTDF.T
        #tt.toc("TRANSPOSED PTDF") 

        #tt.tic("inserting 0 column in PTDF")
        #PTDF = np.insert(PTDF, _ref_bus_idx, np.zeros(_len_branch), axis=1)
        #tt.toc("INSERTED 0 column in PTDF")

else:
    tt.tic(f'calculating ptdf matrix')
    PTDF = get_ptdf(branches,buses,index_set_branch,index_set_bus,reference_bus,mapping_bus_to_idx)
    tt.toc(f'calculated ptdf matrix')

tt_timer.tic("DONE")
'''
