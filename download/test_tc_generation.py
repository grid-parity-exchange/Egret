from egret.data.model_data import ModelData
 
import egret.model_library.transmission.tx_calc as tx_calc
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg

import sys

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

'''
tt.tic("Starting N-1 connection check")
num_connected = 0
branches_unmonitored_sp = set()
for branch_name in branches:
    
    check_contingency_connection = tx_calc.check_contingency_connection(graph, branches, [branch_name], mapping_bus_to_idx)
    if check_contingency_connection:
        num_connected += 1
    #    print(branch_name)
    #    branchn = (branches[branch_name]['from_bus'], branches[branch_name]['to_bus'])
    #    branches_unmonitored_sp.add(branch_name)
    #else:
    #    print("Contingency network not connected using check_contingency_connection")
    
print("Number of N-1 contingencies monitorable", num_connected)
tt.toc("Finished N-1 connection check")

#print(branches_unmonitored_nx-branches_unmonitored_sp)
#print(branches_unmonitored_sp-branches_unmonitored_nx)
'''
    
mdn2 = md.clone()

branch_name = sys.argv[2]
del mdn2.data['elements']['branch'][branch_name]

branchesn2 = mdn2.data['elements']['branch']; busesn2 = mdn2.data['elements']['bus']
index_set_branchn2 = tuple(branchesn2); index_set_busn2 = tuple(busesn2)
Bdn2 = tx_calc._calculate_Bd(branchesn2, index_set_branchn2) 

mapping_bus_to_idxn2 = {b : idx for idx,b in enumerate(index_set_busn2) }
mapping_branch_to_idxn2 = {b : idx for idx,b in enumerate(index_set_branchn2) }

graph2 = tx_calc.construct_connection_graph(branchesn2, mapping_bus_to_idxn2)
connected = tx_calc.check_network_connection(graph2, index_set_busn2)

if connected:
    print("Contingency network connected")
    assert branch_name in branches_not_disconnecting
else:
    print("Contingency network not connected")

if not connected:
    print(f"NETWORK is not connected with branch {branch_name} removed")

assert connected

contingencies = { 'c0' : {'branch_contingency': branch_name} }

MLU_MP, B_dA, ref_bus_mask, contingency_compensators = tx_calc.calculate_ptdf_factorization(branches, buses, index_set_branch, index_set_bus, reference_bus, mapping_bus_to_idx=mapping_bus_to_idx, mapping_branch_to_idx=mapping_branch_to_idx, contingencies=contingencies)

phi_adjust = tx_calc.calculate_phi_adjust(branches,index_set_branch,index_set_bus,mapping_bus_to_idx=mapping_bus_to_idx)
phi_adjust = phi_adjust[ref_bus_mask]

print(contingency_compensators)

Bd = tx_calc._calculate_Bd(branches, index_set_branch)

At = tx_calc.calculate_adjacency_matrix_transpose(branches,index_set_branch, index_set_bus, mapping_bus_to_idx) 

A = At.T

Y = At@Bd@A


Bdn2 = tx_calc._calculate_Bd(branchesn2, index_set_branchn2)
Atn2 = tx_calc.calculate_adjacency_matrix_transpose(branchesn2, index_set_branchn2, index_set_busn2, mapping_bus_to_idxn2)

An2 = Atn2.T

Yn2 = Atn2@Bdn2@An2

# modification line/removal
branch_idx = mapping_branch_to_idx[branch_name]
M1 = A[branch_idx].T

## the *change* is suspectance is -Bd
delY = M1*(-Bd[branch_idx,branch_idx])*M1.T

if not np.allclose((Y+delY).A, Yn2.A):
    print('BASSDFS, mod matrix is not working')
    error = True

M = Y[ref_bus_mask,:][:,ref_bus_mask]

Mn2 = Yn2[ref_bus_mask,:][:,ref_bus_mask]

delM = delY[ref_bus_mask,:][:,ref_bus_mask]

if not np.allclose((M+delM).A, Mn2.A):
    print('BASSDFS, mod matrix without reference bus is not working')
    error = True

MLU_MPn2, B_dAn2, ref_bus_maskn2, _ = tx_calc.calculate_ptdf_factorization(branchesn2, busesn2, index_set_branchn2, index_set_busn2, reference_bus, mapping_bus_to_idx=mapping_bus_to_idxn2, mapping_branch_to_idx=mapping_branch_to_idxn2)

phi_adjust_2 = tx_calc.calculate_phi_adjust(branchesn2,index_set_branchn2,index_set_busn2,mapping_bus_to_idx=mapping_bus_to_idxn2)

phi_adjust_2 = phi_adjust_2[ref_bus_maskn2]

tt.tic()
buff = np.zeros((len(index_set_bus)-1,1))
'''
## need L/U factors like in table I scheme 3
_pr_pc_len = len(index_set_bus) - 1
assert _pr_pc_len == B_dA.shape[1]
Pr = sp.csc_matrix((np.ones(_pr_pc_len), (MLU_MP.perm_r, np.arange(_pr_pc_len))))
Pc = sp.csc_matrix((np.ones(_pr_pc_len), (np.arange(_pr_pc_len), MLU_MP.perm_c)))

#L = Pr.T * MLU_MP.L
#U = MLU_MP.U * Pc.T

#tt.toc('got basic L/U factors')

#assert (np.allclose((L*U).A, M.A))

tt.toc('setup Pr/Pc')
## The functions above are generic for every branch out

## import for this to be negative!!
dely = -Bd[branch_idx, branch_idx]
## Table I Scheme 3
M1_masked = M1[ref_bus_mask]
## shouldn't need to re-order
splu_options = {
                 "Equil":False,
                 "ColPerm":"NATURAL",
                 #"DiagPivotThresh":0.0,
               }
L_factor = sp.linalg.splu(MLU_MP.L,options=splu_options)
U_factor = sp.linalg.splu(MLU_MP.U,options=splu_options)
tt.toc('computed silly L/U factors')
print(f'MLU_MP.L.nnz: {MLU_MP.L.nnz}')
print(f'L_factor.L.nnz: {L_factor.L.nnz}')
print(f'L_factor.U.nnz: {L_factor.U.nnz}')

print(f'MLU_MP.U.nnz: {MLU_MP.U.nnz}')
print(f'U_factor.L.nnz: {U_factor.L.nnz}')
print(f'U_factor.U.nnz: {U_factor.U.nnz}')
#W = sp.linalg.spsolve(MLU_MP.L,Pr@M1_masked,permc_spec="NATURAL")
#tt.toc('allocated array')
#print(f'buff.shape: {buff.shape}')
#in_ = (Pr@M1_masked).toarray(out=buff)
#tt.toc('computed in_ for L')
#
#in_csc = L_factor.solve(in_)
#tt.toc('computed in_csc for csc')
#W = sp.csc_matrix( in_csc )
W = sp.csc_matrix( L_factor.solve((Pr@M1_masked).toarray(out=buff)) )
tt.toc('computed sparse W')
print(f'   W.nnz: {W.nnz}')
#W = L_factor.solve((Pr@M1_masked).A)
## NOTE: In principle, the triangular solve should be more efficient, but
##       scipy's implementation is probably too python-based. Should look
##       into cythonized version of a triangular solve, which should be quicker
##       than spsolve. Further, W/Wbar are often sparse because of the way
##       SuperLU computes L/U (i.e., on a 24464-bus network W/Wbar only have
##       254 nnz
#W = sp.linalg.spsolve_triangular(MLU_MP.L.tocsr(), (Pr@M1_masked).A, unit_diagonal=True, overwrite_b=True)
Wbar = sp.csc_matrix( U_factor.solve((Pc.T@M1_masked).toarray(out=buff), 'T') )
print(f'Wbar.nnz: {Wbar.nnz}')
#Wbar = U_factor.solve((Pc.T@M1_masked).A, 'T')
#Wbar = sp.linalg.spsolve(MLU_MP.U.T, Pc.T@M1_masked,permc_spec="NATURAL")
#Wbar = sp.linalg.spsolve_triangular(MLU_MP.U.T, (Pc.T@M1_masked).A, overwrite_b=True)
tt.toc('computed W/Wbar')

print(f'   W min: {np.min(np.abs(W.data))}')
print(f'Wbar min: {np.min(np.abs(Wbar.data))}')

#z = Wbar.T@W
z = (Wbar.T@W)[0,0]

tt.toc('computed z')

## (7a)
c = 1/(1/dely + z)

tt.toc('computed c')
'''
comp = contingency_compensators['c0']
'''
assert all(W == comp.W)
assert all(Wbar == comp.Wbar)
assert c == comp.c
assert all(M1_masked == comp.M)
'''

W = comp.W
Wbar = comp.Wbar
c = comp.c
M1_masked = comp.M

L_factor = contingency_compensators.L
U_factor = contingency_compensators.U
Pc = contingency_compensators.Pc
Pr = contingency_compensators.Pr

#tt.toc('verified compensator')

# NOTE: IEEE 300 branch 390 has a non-trivial compensator
print(comp.phi_compensator)
assert all(phi_adjust_2 == (phi_adjust + comp.phi_compensator))
tt.toc('verified phi compensator')

p_nw = { bus : 0. for bus in index_set_bus }
for g, g_dict in md.elements(element_type='generator'):
    p_nw[g_dict['bus']] -= g_dict['pg']

for l, l_dict in md.elements(element_type='load'):
    p_nw[l_dict['bus']] += l_dict['p_load']

NWV = np.fromiter((p_nw[b] for b in index_set_bus), float, count=len(index_set_bus))

NWV = NWV[ref_bus_mask]


VA0 = MLU_MP.solve(NWV)

VA0n2 = MLU_MPn2.solve(NWV)

tt.tic('setting up adjust array')

adjust_array = M1_masked*((-c)*(M1_masked.T@VA0))
tt.toc('computed adjust_array')

VA_delta = MLU_MP.solve(adjust_array)
tt.toc('computed VA_delta')

if not np.allclose(VA0+VA_delta, VA0n2):
    print("ADASD problem with angle adjustment")
    error = True
tt.tic('verified VA')

PF0 = B_dA@VA0

tt.toc('computed base-case flows')

PF_delta = B_dA@VA_delta
## zero-out flow on the line
PF_delta[branch_idx] = -PF0[branch_idx]
tt.toc('computed PF_delta')

PF0n2 = B_dAn2@VA0n2

## insert 0 for the line taken out
PF0n2 = np.insert(PF0n2, branch_idx, 0.)

if not np.allclose(PF0+PF_delta, PF0n2, atol=1e-06):
    print("ASDASD problem with flow adjustement")
    error = True

if error:
    print("SOME ISSUE ABOVE")
else:
    print(f"For case {sys.argv[1]}, everything matches with flows!")

error = False

# try using mid-compensation for VA0n2
hatF = L_factor.solve(Pr@NWV)
delF = W*((-c)*(Wbar.T@hatF))
VA0n2_new = Pc@(U_factor.solve(hatF+delF))
#VA0n2_new = Pc@(U_factor.solve((hatF+delF.T)[0]))

if not np.allclose(VA0n2_new, VA0n2):
    print("Problem with basic mid-compensentaation")
    error = True

# branch PTDFs
row_idx = 2
assert branch_idx != row_idx
if branch_idx > row_idx:
    offset = 0
else:
    offset = -1

assert np.allclose(B_dA[row_idx].A[0], B_dAn2[row_idx+offset].A[0])

tt.tic('looking at PTDFs')

buff.shape = (1,len(index_set_bus)-1)
PTDF_row = MLU_MP.solve(B_dA[row_idx].toarray(out=buff)[0],'T')
tt.toc('got row from base case')

PTDF_rown2 = MLU_MPn2.solve(B_dA[row_idx].toarray(out=buff)[0],'T')
tt.toc('got row from augmented case')

hatF = U_factor.solve((B_dA[row_idx]@Pc).toarray(out=buff)[0], 'T')
#tt.toc('U_factor')
delF = Wbar*((-c)*(W.T@hatF))
#tt.toc('delF')

PTDF_rown2_new = Pr.T@L_factor.solve(hatF+delF, 'T')
#PTDF_rown2_new = Pr.T@L_factor.solve((hatF+delF.T)[0], 'T')
#tt.toc('L_factor')
tt.toc('got augmented case row from base case')

if not np.allclose(PTDF_rown2_new, PTDF_rown2):
    print("Problem with PTDF construction")
    error = True

if error:
    print("Some issue with PTDF construction")
else:
    print("PTDFs match!")
