import sys
import gc

from pyutilib.misc.timing import TicTocTimer

fn = sys.argv[1]
tt_timer = TicTocTimer()

if sys.argv[2] == 'factor':
    use_factorization = True
elif sys.argv[2] == 'ptdf':
    use_factorization = False
else:
    sys.exit(1)

from egret.data.model_data import ModelData
from egret.model_library.transmission.tx_calc import calculate_ptdf_factorization, calculate_ptdf
from egret.data.ptdf_utils import PTDFMatrix, VirtualPTDFMatrix
from egret.model_library.defn import BasePointType
from egret.common.lazy_ptdf_utils import populate_default_ptdf_options

tt_timer.tic(f'reading file {fn}')
md = ModelData.read(fn)

tt_timer.tic(f'setting up dictionaries')
buses = dict(md.elements(element_type = 'bus'))
branches = dict(md.elements(element_type = 'branch'))

reference_bus = md.data['system']['reference_bus']

index_set_bus = tuple(buses.keys())
index_set_branch = tuple(branches.keys())

bus_index_mapping = { bus_n : i for i, bus_n in enumerate(index_set_bus) }

tt_timer.tic(f'skeleton pyomo model')
import pyomo.environ as pe
mb = pe.ConcreteModel()
_p_nw = { bus : 0. for bus in index_set_bus }

for g, g_dict in md.elements(element_type='generator'):
    _p_nw[g_dict['bus']] -= g_dict['pg']

for l, l_dict in md.elements(element_type='load'):
    _p_nw[l_dict['bus']] += l_dict['p_load']

mb.p_nw = pe.Param(index_set_bus, mutable=True, initialize=_p_nw)

if use_factorization:
    #tt_timer.tic(f'calculating ptdf factorization')
    #calculate_ptdf_factorization(branches,buses,index_set_branch,index_set_bus,reference_bus,mapping_bus_to_idx=bus_index_mapping)
    #tt_timer.toc(f'calculated ptdf factorization')
    import scipy.linalg as la

    tt_timer.tic(f'Calculating Virtual PTDF Matrix; then PTDF matrix')
    PTDF = VirtualPTDFMatrix(branches, buses, reference_bus, BasePointType.FLATSTART, populate_default_ptdf_options(None))
    self = PTDF
    b_da = self.B_dA.A.T
    #tt_timer.tic(f'Calling lu_solve')
    PTDFM = la.lu_solve(self.MLU, b_da, trans=1, overwrite_b=True, check_finite=False)
    tt_timer.toc(f'Calculated PTDF Matrix!')

    '''

    tt_timer.tic(f'getting a PTDF row')
    branch_name = index_set_branch[0]
    tt_timer.tic(f'looking up index')
    branch_idx = self._branchname_to_index_map[branch_name]
    tt_timer.tic(f'munging B_dA')
    b_da = self.B_dA[branch_idx].A[0]
    tt_timer.toc(f'<- TIME FOR making B_dA dense')
    tt_timer.tic(f'Calling lu_solve')
    PTDF_row = la.lu_solve(self.MLU, b_da, trans=1, overwrite_b=True, check_finite=False)
    tt_timer.toc(f'<- TIME FOR lu_solve')
    tt_timer.tic(f'Transposing row')
    PTDF_row = PTDF_row.T
    tt_timer.toc(f'Transposed')
    tt_timer.tic(f'adding to dictionary cache')
    self._ptdf_rows[branch_name] = PTDF_row
    tt_timer.tic(f'making dict for pyomo')
    { bus_n : val for bus_n, val in zip(self.buses_keys, PTDF_row) }
    tt_timer.toc(f'GOT PTDF row!')

    tt_timer.tic(f'calculating system flows')
    PFV, _ = PTDF.calculate_PFV(mb)
    tt_timer.toc(f'calculated system flows')

    tt_timer.tic(f'calculating whole PTDF matrix using back solve')
    tt_timer.tic(f'transposing/densifying B_dA')
    b_da = self.B_dA.A.T
    tt_timer.tic(f'Calling lu_solve')
    PTDF = la.lu_solve(self.MLU, b_da, trans=1, overwrite_b=True, check_finite=False)
    tt_timer.toc(f'<- TIME FOR lu_solve')
    '''

else:
    #tt_timer.tic(f'calculating ptdf matrix')
    #calculate_ptdf(branches,buses,index_set_branch,index_set_bus,reference_bus,mapping_bus_to_idx=bus_index_mapping)
    #tt_timer.toc(f'calculated ptdf matrix')

    tt_timer.tic(f'Calculating PTDF Matrix; initializing data')
    PTDF = PTDFMatrix(branches, buses, reference_bus, BasePointType.FLATSTART, populate_default_ptdf_options(None))
    tt_timer.toc(f'Calculated PTDF Matrix!')

    print(PTDF.PTDFM)

    '''
    tt_timer.tic(f'getting a PTDF row')
    { bus_n : val for bus_n, val in PTDF.get_branch_ptdf_iterator(index_set_branch[0]) }
    tt_timer.toc(f'GOT PTDF row!')

    tt_timer.tic(f'calculating system flows')
    PFV, _ = PTDF.calculate_PFV(mb)
    tt_timer.toc(f'calculated system flows')
    '''

tt_timer.tic("DONE")
