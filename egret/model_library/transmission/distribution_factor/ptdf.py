from pyomo.core.base.block import _BlockData, Block
import pyomo.environ as pe
from collections import Iterable
from pyomo.core.kernel.component_map import ComponentMap
from pyomo.core.kernel.component_set import ComponentSet
import egret.model_library.transmission.tx_utils as tx_utils
from egret.model_library.defn import DistributionFactorType, BasePointType, ApproximationType
import sys
from .distribution_factor_base import BaseDistributionFactorData
from egret.model_library.transmission.tx_calc import calculate_J11, calculate_adjacency_matrix
import numpy as np
import scipy as sp
"""
Base classes for distribution factor approximations
"""

def _calculate_ptdf(branches,buses,reference_bus,base_point=BasePointType.FLATSTART,sparse_index_set_branch=None):
    """
    Calculates the sensitivity of the voltage angle to real power injections

    Parameters
    ----------
    branches: dict{}
        The dictionary of branches for the test case
    buses: dict{}
        The dictionary of buses for the test case
    reference_bus: key value
        The reference bus key value
    base_point: egret.model_library_defn.BasePointType
        The base-point type for calculating the PTDF matrix
    sparse_index_set_branch: list
        The list of keys for branches needed to compute a sparse PTDF matrix
    """

    index_set_bus = tuple(buses.keys())
    _len_bus = len(index_set_bus)
    _mapping_bus = {i: index_set_bus[i] for i in list(range(0,_len_bus))}

    index_set_branch = tuple(branches.keys())
    _len_branch = len(index_set_branch)
    _mapping_branch = {i: index_set_branch[i] for i in list(range(0,_len_branch))}

    _ref_bus_idx = [key for key, value in _mapping_bus.items() if value == reference_bus][0]

    J = calculate_J11(branches,buses,base_point,approximation_type=ApproximationType.PTDF)
    A = calculate_adjacency_matrix(branches,index_set_branch,index_set_bus)
    M = np.matmul(A.transpose(),J)

    J0 = np.zeros((_len_bus + 1, _len_bus + 1))
    J0[:-1, :-1] = M
    J0[-1][_ref_bus_idx] = 1
    J0[_ref_bus_idx][-1] = 1

    if sparse_index_set_branch is None or len(sparse_index_set_branch) == _len_branch:
        try:
            SENSI = sp.linalg.inv(J0)
        except sp.linalg.LinAlgError:
            print("Matrix not invertible. Calculating pseudo-inverse instead.")
            SENSI = sp.linalg.pinv(J0,rcond=1e-7)
            pass
        SENSI = SENSI[:-1,:-1]

        PTDF = sp.sparse.lil_matrix(sp.matmul(J,SENSI))
    elif len(sparse_index_set_branch) < _len_branch:
        n, m = M.shape
        B = np.array([], dtype=np.int64).reshape(_len_bus + 1,0)
        _sparse_mapping_branch = {i: index_set_branch[i] for i in list(range(0, _len_branch)) if index_set_branch[i] in sparse_index_set_branch}

        for idx, branch_name in _sparse_mapping_branch.items():
            b = np.zeros((_len_branch,1))
            b[idx] = 1
            _tmp = np.matmul(J.transpose(),b)
            _tmp = np.vstack([_tmp,0])
            B = np.concatenate((B,_tmp), axis=1)
        PTDF = sp.sparse.lil_matrix((_len_branch,_len_bus))
        row_idx = list(_sparse_mapping_branch.keys())
        _ptdf = sp.sparse.linalg.spsolve(sp.sparse.csr_matrix(J0.transpose()), B).T
        PTDF[row_idx] = _ptdf[:,:-1]

    return PTDF


@declare_custom_block(name='PTDF')
class PTDFData(BaseDistributionFactorData):
    def __init__(self, component, model_data, base_point_type):
        BaseDistributionFactorData.__init__(self, component)

        self.__distribution_factor_type = DistributionFactorType.PTDF

        if base_point_type not in BasePointType:
            raise ValueError('{0} is not a valid member of BasePointType'.format(base_point_type))
        self.__base_point_type = base_point_type

        md = model_data.clone_in_service()
        tx_utils.scale_ModelData_to_pu(md, inplace = True)
        self.__reference_bus = md.data['system']['reference_bus']
        self.__buses = dict(md.elements(element_type='bus'))
        self.__branches = dict(md.elements(element_type='branch'))
        self._sparse_branches = None

        self._ptdf_m = None

    @property
    def distribution_factor_type(self):
        return self.__distribution_factor_type

    @property
    def base_point_type(self):
        return self.__base_point_type

    @property
    def reference_bus(self):
        return self.__reference_bus

    @property
    def buses(self):
        return self.__buses

    @property
    def branches(self):
        return self.__branches

    @property
    def sparse_branches(self):
        return self._sparse_branches

    @sparse_branches.setter
    def sparse_branches(self, val):
        if tuple(self._sparse_branches.keys()) not in tuple(self._branches.keys()):
            raise ValueError('Sparse branches not a subset of branches.')
        self._sparse_branches = val

    @property
    def ptdf_mat(self):
        return self._ptdf_m

    def build(self):
        self._ptdf_m = _calculate_ptdf(self._branches, self.__buses, self.__base_point_type, self._sparse_branches)
