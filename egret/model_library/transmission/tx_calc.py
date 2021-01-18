#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module collects some helper functions useful for performing
different computations for transmission models
"""
import math
import weakref
import collections.abc as abc
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
import networkx as nx

from math import cos, sin
from egret.model_library.defn import BasePointType, ApproximationType
from egret.common.log import logger

def calculate_conductance(branch):
    rs = branch['resistance']
    xs = branch['reactance']
    return rs / (rs**2 + xs**2)


def calculate_susceptance(branch):
    rs = branch['resistance']
    xs = branch['reactance']
    return -xs / (rs**2 + xs**2)


def calculate_y_matrix_from_branch(branch):
    rs = branch['resistance']
    xs = branch['reactance']
    bs = branch['charging_susceptance']
    tau = 1.0
    shift = 0.0
    if branch['branch_type'] == 'transformer':
        tau = branch['transformer_tap_ratio']
        shift = branch['transformer_phase_shift']
    return calculate_y_matrix(rs, xs, bs, tau, shift)


def calculate_y_matrix(rs, xs, bc, tau, shift):
    """
    Compute the y matrix from various branch properties

    Parameters
    ----------
    rs : float
        Branch resistance
    xs : float
        Branch reactance
    bc : float
        Branch charging susceptance
    tau : float
        Branch transformer tap ratio
    shift : float
        Branch transformer phase shift

    Returns
    -------
        list : list of floats representing the y matrix
               [Y(ifr,vfr), Y(ifr,vfj), Y(ifr,vtr), Y(ifr,vtj),
               Y(ifj,vfr), Y(ifj,vfj), Y(ifj,vtr), Y(ifj,vtj),
               Y(itr,vfr), Y(itr,vfj), Y(itr,vtr), Y(itr,vtj),
               Y(itj,vfr), Y(itj,vfj), Y(itj,vtr), Y(itj,vtj)]
    """
    bc = bc/2
    tr = tau * math.cos(math.radians(shift))
    tj = tau * math.sin(math.radians(shift))
    mag = rs**2 + xs**2

    a = rs/(tau**2*mag)                    # c1
    b = (1/tau**2) * (xs/mag - bc)         # c2
    c = (-rs*tr - xs*tj)/(tau**2 * mag)    # c3
    d = (rs*tj - xs*tr)/(tau**2 * mag)     # c4
    e = -b                                 # -c2
    f = a                                  # c1
    g = -d                                 # -c4
    h = c                                  # c3
    i = (xs*tj - rs*tr)/(tau**2 * mag)     # c7
    j = (-rs*tj - xs*tr)/(tau**2 * mag)    # c8
    k = rs/mag                             # c5
    l = xs/mag - bc                        # c6
    m = -j                                 # -c8
    n = i                                  # c7
    o = -l                                 # -c6
    p = k                                  # c5

    # y = [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p]
    y_dict = {}
    y_dict[('ifr', 'vfr')] = a
    y_dict[('ifr', 'vfj')] = b
    y_dict[('ifr', 'vtr')] = c
    y_dict[('ifr', 'vtj')] = d
    
    y_dict[('ifj', 'vfr')] = e
    y_dict[('ifj', 'vfj')] = f
    y_dict[('ifj', 'vtr')] = g
    y_dict[('ifj', 'vtj')] = h
    
    y_dict[('itr', 'vfr')] = i
    y_dict[('itr', 'vfj')] = j
    y_dict[('itr', 'vtr')] = k
    y_dict[('itr', 'vtj')] = l

    y_dict[('itj', 'vfr')] = m
    y_dict[('itj', 'vfj')] = n
    y_dict[('itj', 'vtr')] = o
    y_dict[('itj', 'vtj')] = p

    return y_dict

def calculate_ifr(vfr, vfj, vtr, vtj, y_matrix):
    """
    Compute ifr from voltages and the y_matrix (computed
    from the branch properties using :py:meth:`calculate_branch_y_matrix`)
    """
    ifr = y_matrix['ifr', 'vfr'] * vfr + y_matrix['ifr', 'vfj'] * vfj + \
        y_matrix['ifr', 'vtr'] * vtr + y_matrix['ifr', 'vtj'] * vtj
    return ifr


def calculate_ifj(vfr, vfj, vtr, vtj, y_matrix):
    """
    Compute ify from voltages and the y_matrix (computed
    from the branch properties using :py:meth:`calculate_branch_y_matrix`)
    """
    ifj = y_matrix['ifj', 'vfr'] * vfr + y_matrix['ifj', 'vfj'] * vfj + \
        y_matrix['ifj', 'vtr'] * vtr + y_matrix['ifj', 'vtj'] * vtj
    return ifj


def calculate_itr(vfr, vfj, vtr, vtj, y_matrix):
    """
    Compute itr from voltages and the y_matrix (computed
    from the branch properties using :py:meth:`calculate_branch_y_matrix`)
    """
    itr = y_matrix['itr', 'vfr'] * vfr + y_matrix['itr', 'vfj'] * vfj + \
        y_matrix['itr', 'vtr'] * vtr + y_matrix['itr', 'vtj'] * vtj
    return itr


def calculate_itj(vfr, vfj, vtr, vtj, y_matrix):
    """
    Compute itj from voltages and the y_matrix (computed
    from the branch properties using :py:meth:`calculate_branch_y_matrix`)
    """
    itj = y_matrix['itj', 'vfr'] * vfr + y_matrix['itj', 'vfj'] * vfj + \
        y_matrix['itj', 'vtr'] * vtr + y_matrix['itj', 'vtj'] * vtj
    return itj


def calculate_ir(p, q, vr, vj):
    """
    Compute ir from power flows and voltages
    """
    ir = (q*vj+p*vr)/(vj**2 + vr**2)
    return ir


def calculate_ij(p, q, vr, vj):
    """
    Compute ij from power flows and voltages
    """
    ij = (p*vj-q*vr)/(vj**2 + vr**2)
    return ij


def calculate_p(ir, ij, vr, vj):
    """
    Compute real power flow from currents and voltages
    """
    p = vr * ir + vj * ij
    return p


def calculate_q(ir, ij, vr, vj):
    """
    Compute reactive power flow from currents and voltages
    """
    q = vj * ir - vr * ij
    return q


def calculate_vr_from_vm_va(vm, va):
    """
    Compute the value of vr from vm and va

    Parameters
    ----------
    vm : float
        The value of voltage magnitude (per)
    va : float
        The value of voltage angle (degrees)

    Returns
    -------
        float : the value of vr or None if
           either vm or va (or both) is None
    """
    if vm is not None and va is not None:
        vr = vm * math.cos(va*math.pi/180)
        return vr
    return None


def calculate_vj_from_vm_va(vm, va):
    """
    Compute the value of vj from vm and va

    Parameters
    ----------
    vm : float
        The value of voltage magnitude (per)
    va : float
        The value of voltage angle (degrees)

    Returns
    -------
        float : the value of vj or None if
           either vm or va (or both) is None
    """
    if vm is not None and va is not None:
        vj = vm * math.sin(va*math.pi/180)
        return vj
    return None


def calculate_vm_from_vj_vr(vj,vr):
    """
    Compute the value of vm from vj and vr

    Parameters
    ----------
    vj : float
        The value of the imaginary part of the voltage phasor (per)
    vr : float
        The value of the real part of the voltage phasor (per)

    Returns
    -------
        float : the value of the voltage magnitude vm or None if
           either vj or vr (or both) is None
    """
    if vj is not None and vr is not None:
        vm = math.sqrt(vj**2 + vr**2)
        return vm
    return None


def calculate_va_from_vj_vr(vj, vr):
    """
    Compute the value of va from vj and vr

    Parameters
    ----------
    vj : float
        The value of the imaginary part of the voltage phasor (per)
    vr : float
        The value of the real part of the voltage phasor (per)

    Returns
    -------
        float : the value of the voltage angle va in degrees or None if
           either vj or vr (or both) is None
    """
    if vj is not None and vr is not None:
        va = math.degrees(math.atan(vj/vr))
        return va
    return None

def _get_susceptance(branch, approximation_type):
    if branch['branch_type'] == 'transformer':
        tau = branch['transformer_tap_ratio']
    else:
        tau = 1.

    if approximation_type == ApproximationType.PTDF:
        x = branch['reactance']
        b = -1./(tau*x)
    elif approximation_type == ApproximationType.PTDF_LOSSES:
        b = calculate_susceptance(branch)/tau
    else:
        raise RuntimeError("Could not find appropriate susceptance value")
    return b

def _calculate_J11(branches,buses,index_set_branch,index_set_bus,mapping_bus_to_idx,base_point=BasePointType.FLATSTART,approximation_type=ApproximationType.PTDF):
    """
    Compute the power flow Jacobian for partial derivative of real power flow to voltage angle
    """
    _len_bus = len(index_set_bus)
    _len_branch = len(index_set_branch)

    data = []
    row = []
    col = []

    for idx_row, branch_name in enumerate(index_set_branch):
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        b = _get_susceptance(branch, approximation_type)

        if base_point == BasePointType.FLATSTART:
            val = -b

        elif base_point == BasePointType.SOLUTION: # TODO: check that we are loading the correct values (or results)
            vn = buses[from_bus]['vm']
            vm = buses[to_bus]['vm']
            tn = buses[from_bus]['va']
            tm = buses[to_bus]['va']

            val = -b * vn * vm * cos(tn - tm)

        idx_col = mapping_bus_to_idx[from_bus]
        row.append(idx_row)
        col.append(idx_col)
        data.append(val)

        idx_col = mapping_bus_to_idx[to_bus]
        row.append(idx_row)
        col.append(idx_col)
        data.append(-val)

    J11 = sp.coo_matrix( (data, (row,col)), shape=(_len_branch, _len_bus))
    return J11.tocsc()

def _calculate_Bd(branches,index_set_branch,base_point=BasePointType.FLATSTART,approximation_type=ApproximationType.PTDF):
    """
    Compute the power flow Jacobian for partial derivative of real power flow to voltage angle
    """
    _len_branch = len(index_set_branch)

    data = []
    row = []
    col = []

    for idx_row, branch_name in enumerate(index_set_branch):
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        b = _get_susceptance(branch, approximation_type)

        if base_point == BasePointType.FLATSTART:
            val = b

        elif base_point == BasePointType.SOLUTION: # TODO: check that we are loading the correct values (or results)
            vn = buses[from_bus]['vm']
            vm = buses[to_bus]['vm']
            tn = buses[from_bus]['va']
            tm = buses[to_bus]['va']

            val = b * vn * vm * cos(tn - tm)

        data.append(val)
        row.append(idx_row)
        col.append(idx_row)

    Bd = sp.coo_matrix( (data, (row,col)), shape=(_len_branch, _len_branch))
    return Bd.tocsc()


def _calculate_L11(branches,buses,index_set_branch,index_set_bus,mapping_bus_to_idx,base_point=BasePointType.FLATSTART):
    """
    Compute the power flow Jacobian for partial derivative of real power losses to voltage angle
    """
    _len_bus = len(index_set_bus)
    _len_branch = len(index_set_branch)

    row = []
    col = []
    data = []

    for idx_row, branch_name in enumerate(index_set_branch):
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        tau = 1.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
        g = calculate_conductance(branch)/tau

        if base_point == BasePointType.FLATSTART:
            vn = 1.
            vm = 1.
            tn = 0.
            tm = 0.
        elif base_point == BasePointType.SOLUTION: # TODO: check that we are loading the correct values (or results)
            vn = buses[from_bus]['vm']
            vm = buses[to_bus]['vm']
            tn = buses[from_bus]['va']
            tm = buses[to_bus]['va']

        val = 2 * g * vn * vm * sin(tn - tm)

        idx_col = mapping_bus_to_idx[from_bus]
        row.append(idx_row)
        col.append(idx_col)
        data.append(val)

        idx_col = mapping_bus_to_idx[to_bus]
        row.append(idx_row)
        col.append(idx_col)
        data.append(-val)

    L11 = sp.coo_matrix((data,(row,col)),shape=(_len_branch,_len_bus))
    return L11.tocsr()

def calculate_phase_shift_flow_adjuster(branches, index_set_branch):
    _len_branch = len(index_set_branch)
    row = []
    col = []
    data = []

    idx_col = 0

    for idx_row, branch_name in enumerate(index_set_branch):
        branch = branches[branch_name]
        if branch['branch_type'] == 'transformer' and branch['transformer_phase_shift'] != 0.:
            val = -(1./branch['reactance']) * (math.radians(branch['transformer_phase_shift'])/branch['transformer_tap_ratio'])

            row.append(idx_row)
            col.append(idx_col)
            data.append(val)

    phase_shift_flow_adjuster = sp.coo_matrix((data, (row,col)), shape=(_len_branch,1))
    return phase_shift_flow_adjuster.tocsc()

def calculate_phi_constant(branches,index_set_branch,index_set_bus,approximation_type=ApproximationType.PTDF, mapping_bus_to_idx=None):
    """
    Compute the phase shifter constant for fixed phase shift transformers
    """
    _len_bus = len(index_set_bus)

    if mapping_bus_to_idx is None:
        mapping_bus_to_idx = {bus_n: i for i, bus_n in enumerate(index_set_bus)}

    _len_branch = len(index_set_branch)

    row_from = []
    row_to = []
    col = []
    data = []

    for idx_col, branch_name in enumerate(index_set_branch):
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        if branch['branch_type'] == 'transformer':
            shift = math.radians(branch['transformer_phase_shift'])
        else: # shift == 0
            continue

        b = _get_susceptance(branch, approximation_type)
        b *= shift

        row_from.append(mapping_bus_to_idx[from_bus])
        row_to.append(mapping_bus_to_idx[to_bus])
        col.append(idx_col)
        data.append(b)

    phi_from = sp.coo_matrix((data,(row_from,col)), shape=(_len_bus,_len_branch))
    phi_to = sp.coo_matrix((data,(row_to,col)), shape=(_len_bus,_len_branch))

    return phi_from.tocsr(), phi_to.tocsr()

def calculate_phi_adjust(branches,index_set_branch,index_set_bus,approximation_type=ApproximationType.PTDF, mapping_bus_to_idx=None):
    """
    Compute the phase shifter constant for fixed phase shift transformers
    """
    _len_bus = len(index_set_bus)

    if mapping_bus_to_idx is None:
        mapping_bus_to_idx = {bus_n: i for i, bus_n in enumerate(index_set_bus)}

    _len_branch = len(index_set_branch)

    row = []
    col = []
    data = []

    for idx_col, branch_name in enumerate(index_set_branch):
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        if branch['branch_type'] == 'transformer' and branch['transformer_phase_shift'] != 0.:
            shift = math.radians(branch['transformer_phase_shift'])
        else: # shift == 0
            continue

        b = _get_susceptance(branch, approximation_type)
        b *= shift

        row.append(mapping_bus_to_idx[from_bus])
        col.append(0)
        data.append(b)

        row.append(mapping_bus_to_idx[to_bus])
        col.append(0)
        data.append(-b)

    phi_adjust = sp.coo_matrix((data,(row,col)), shape=(_len_bus,1))

    return phi_adjust.tocsc()


def calculate_phi_loss_constant(branches,index_set_branch,index_set_bus,approximation_type=ApproximationType.PTDF_LOSSES, mapping_bus_to_idx=None):
    """
    Compute the phase shifter constant for fixed phase shift transformers
    """
    _len_bus = len(index_set_bus)

    if mapping_bus_to_idx is None:
        mapping_bus_to_idx = {bus_n: i for i, bus_n in enumerate(index_set_bus)}

    _len_branch = len(index_set_branch)

    row_from = []
    row_to = []
    col = []
    data = []

    for idx_col, branch_name in enumerate(index_set_branch):
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        tau = 1.0
        shift = 0.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])

        g = 0.
        if approximation_type == ApproximationType.PTDF:
            r = branch['resistance']
            g = (1/r)*(1/tau)*shift**2
        elif approximation_type == ApproximationType.PTDF_LOSSES:
            g = calculate_conductance(branch)*(1/tau)*shift**2

        row_from.append(mapping_bus_to_idx[from_bus])
        row_to.append(mapping_bus_to_idx[to_bus])
        col.append(idx_col)
        data.append(g)

    phi_loss_from = sp.coo_matrix((data,(row_from,col)),shape=(_len_bus,_len_branch))
    phi_loss_to = sp.coo_matrix((data,(row_to,col)),shape=(_len_bus,_len_branch))

    return phi_loss_from.tocsr(), phi_loss_to.tocsr()

def _calculate_pf_constant(branches,buses,index_set_branch,base_point=BasePointType.FLATSTART):
    """
    Compute the power flow constant for the taylor series expansion of real power flow as
    a convex combination of the from/to directions, i.e.,
    pf = 0.5*g*((tau*vn)^2 - vm^2) - tau*vn*vm*b*sin(tn-tm-shift)
    """

    _len_branch = len(index_set_branch)
    ## this will be fully dense
    pf_constant = np.zeros(_len_branch)

    for idx_row, branch_name in enumerate(index_set_branch):
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        tau = 1.0
        shift = 0.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])
        g = calculate_conductance(branch)
        b = calculate_susceptance(branch)/tau

        if base_point == BasePointType.FLATSTART:
            vn = 1.
            vm = 1.
            tn = 0.
            tm = 0.
        elif base_point == BasePointType.SOLUTION: # TODO: check that we are loading the correct values (or results)
            vn = buses[from_bus]['vm']
            vm = buses[to_bus]['vm']
            tn = buses[from_bus]['va']
            tm = buses[to_bus]['va']

        pf_constant[idx_row] = 0.5 * g * ((vn/tau) ** 2 - vm ** 2) \
                               - b * vn * vm * (sin(tn - tm + shift) - cos(tn - tm + shift)*(tn - tm))

    return pf_constant


def _calculate_pfl_constant(branches,buses,index_set_branch,base_point=BasePointType.FLATSTART):
    """
    Compute the power losses constant for the taylor series expansion of real power losses as
    a convex combination of the from/to directions, i.e.,
    pfl = g*((tau*vn)^2 + vm^2) - 2*tau*vn*vm*g*cos(tn-tm-shift)
    """

    _len_branch = len(index_set_branch)

    ## this will be fully dense
    pfl_constant = np.zeros(_len_branch)

    for idx_row, branch_name in enumerate(index_set_branch):
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        tau = 1.0
        shift = 0.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])
        _g = calculate_conductance(branch)
        g = _g/tau
        g2 = _g/tau**2

        if base_point == BasePointType.FLATSTART:
            vn = 1.
            vm = 1.
            tn = 0.
            tm = 0.
        elif base_point == BasePointType.SOLUTION: # TODO: check that we are loading the correct values (or results)
            vn = buses[from_bus]['vm']
            vm = buses[to_bus]['vm']
            tn = buses[from_bus]['va']
            tm = buses[to_bus]['va']

        pfl_constant[idx_row] = g2 * (vn ** 2) + _g * (vm ** 2) \
                              - 2 * g * vn * vm * (sin(tn - tm + shift) * (tn - tm) + cos(tn - tm + shift))

    return pfl_constant

def calculate_ptdf_factorization(branches,buses,index_set_branch,index_set_bus,reference_bus,
                                 base_point=BasePointType.FLATSTART,
                                 contingencies=None,
                                 mapping_bus_to_idx=None,
                                 mapping_branch_to_idx=None,
                                 interfaces=None,
                                 index_set_interface=None):

    if interfaces is None:
        assert index_set_interface is None
    if index_set_interface is None:
        assert interfaces is None

    _len_bus = len(index_set_bus)

    if mapping_bus_to_idx is None:
        mapping_bus_to_idx = {bus_n: i for i, bus_n in enumerate(index_set_bus)}

    _len_branch = len(index_set_branch)

    _ref_bus_idx = mapping_bus_to_idx[reference_bus]

    ## check if the network is connected
    graph = construct_connection_graph(branches, mapping_bus_to_idx)
    connected = check_network_connection(graph, index_set_bus)

    if not connected:
        raise RuntimeError("Network is not connected, cannot use PTDF formulation")

    #(A^T)
    At = calculate_adjacency_matrix_transpose(branches,index_set_branch,index_set_bus,mapping_bus_to_idx)

    ref_bus_mask = np.ones(_len_bus, dtype=bool)
    ref_bus_mask[_ref_bus_idx] = False

    At_masked = At[ref_bus_mask]

    Bd = _calculate_Bd(branches, index_set_branch)
    B_dA = Bd@(At_masked.T)

    # M is now (A^T B_d A) with
    # row and column of reference
    # bus removed
    M = At_masked@B_dA

    ## LU factorization
    MLU_MP = scipy.sparse.linalg.splu(M)

    if contingencies:
        contingency_compensators = \
            precompute_contingency_matricies( graph, MLU_MP, At_masked.T, Bd,
                                              mapping_bus_to_idx, mapping_branch_to_idx, 
                                              ref_bus_mask,
                                              branches, contingencies )
    else:
        contingency_compensators = {}

    if interfaces is None:
        return MLU_MP, B_dA, ref_bus_mask, contingency_compensators
    else:
        if mapping_bus_to_idx is None:
            mapping_branch_to_idx = {branch_n: i for i, branch_n in enumerate(index_set_branch)}
        I = _calculate_interface_matrix(interfaces, index_set_interface, mapping_branch_to_idx)
        B_dA_I = I@B_dA

        return MLU_MP, B_dA, ref_bus_mask, contingency_compensators, B_dA_I, I

class _ContingencyCompensator:
    def __init__(self, M, c, W, Wbar, phi_compensator, VA_compensator, branch_out):
        self._M = M
        self._c = c
        self._W = W
        self._Wbar = Wbar
        self._phi_compensator = phi_compensator
        self._VA_compensator = VA_compensator
        self._branch_out = branch_out
        self._global = None

    @property
    def M(self):
        return self._M
    @property
    def c(self):
        return self._c
    @property
    def W(self): 
        return self._W
    @property
    def Wbar(self): 
        return self._Wbar
    @property
    def phi_compensator(self):
        return self._phi_compensator
    @property
    def VA_compensator(self):
        return self._VA_compensator
    @property
    def branch_out(self):
        return self._branch_out

    @property
    def L(self):
        return self._global()._L
    @property
    def U(self):
        return self._global()._U
    @property
    def Pr(self):
        return self._global()._Pr
    @property
    def Pc(self):
        return self._global()._Pc

class _ContingencyCompensators(abc.Mapping):
    def __init__(self, compensators, L, U, Pr, Pc):
        self._compensators = compensators
        for c in compensators.values():
            c._global = weakref.ref(self)
        self._L = L
        self._U = U
        self._Pr = Pr
        self._Pc = Pc

    def __getitem__(self, key):
        return self._compensators[key]

    def __iter__(self):
        return iter(self._compensators)

    def __len__(self):
        return len(self._compensators)

    @property
    def L(self):
        return self._L
    @property
    def U(self):
        return self._U
    @property
    def Pr(self): 
        return self._Pr
    @property
    def Pc(self): 
        return self._Pc


def precompute_contingency_matricies( graph, MLU_MP, A, Bd,\
                                      mapping_bus_to_idx,  mapping_branch_to_idx, 
                                      ref_bus_mask,
                                      branches, contingencies):

    contingencies_monitored = {}
    for c, cdict in contingencies.items():
        if 'branch_contingency' not in cdict:
            logger.warning(f"Contingency {c} does not have a branch specified; ignoring")
            continue
        branches_out = cdict['branch_contingency'] 
        if isinstance( branches_out, list ):
            if len(branches_out) == 0:
                logger.warning(f"Contingency {c} does not have a branch specified; ignoring")
                continue
            if len(branches_out) > 1:
                raise RuntimeError(f"Contingency {c} has multiple branches. This is not currently supported")
            branch_out = branches_out[0]
            if branch_out not in branches:
                raise RuntimeError(f"Contingency {c} is already out!")
        elif branches_out in mapping_branch_to_idx:
            branch_out = branches_out
        else:
            raise RuntimeError(f"Contingencies must be specified as a list of branches or single branch")

        contingencies_monitored[c] = branch_out

    _check_contingencies_not_disconnecting(graph, branches, mapping_bus_to_idx, contingencies_monitored.values()) 
    
    ## things for every possible modification
    _bus_len = A.shape[1]
    Pr = sp.csc_matrix((np.ones(_bus_len), (MLU_MP.perm_r, np.arange(_bus_len))))
    Pc = sp.csc_matrix((np.ones(_bus_len), (np.arange(_bus_len), MLU_MP.perm_c)))

    ## shouldn't need to re-order
    splu_options = {
                     "Equil":False,
                     "ColPerm":"NATURAL",
                     #"DiagPivotThresh":0.0,
                   }
    L_factor = sp.linalg.splu(MLU_MP.L,options=splu_options)
    U_factor = sp.linalg.splu(MLU_MP.U,options=splu_options)

    buff = np.zeros((_bus_len,1))

    compensators = {}

    for cn, branch_out in contingencies_monitored.items():
        branch_out_idx = mapping_branch_to_idx[branch_out]

        M = A[branch_out_idx].T
        dely = -Bd[branch_out_idx, branch_out_idx]

        # NOTE: The conversions involved here are a bottleneck. 
        #       Egret should probably implement its own sparse 
        #       triangular solver. Batching (collecting Pr@M,
        #       Pc.T@M for every branch_out) could also be tried.
        W = sp.csc_matrix( L_factor.solve((Pr@M).toarray(out=buff)) )
        Wbar = sp.csc_matrix( U_factor.solve((Pc.T@M).toarray(out=buff), 'T') )

        # NOTE: With a single change, these are simple inverses.
        #       If we go to multiple contingencies, this needs to 
        #       use matrix inverses and the code should be re-visited.
        z = (Wbar.T@W)[0,0]
        c = 1./((1./dely) + z)

        # Compute phi_compensator
        branch = branches[branch_out]

        if branch['branch_type'] == 'transformer' and branch['transformer_phase_shift'] != 0.:
            shift = math.radians(branch['transformer_phase_shift'])

            neg_b = shift*dely

            row = [mapping_bus_to_idx[branch['from_bus']], mapping_bus_to_idx[branch['to_bus']]]
            col = [0, 0]
            data = [neg_b, -neg_b]

            phi_comp = sp.coo_matrix((data,(row,col)), shape=(_bus_len+1,1)).tocsc()[ref_bus_mask]
            VA_comp = MLU_MP.solve(phi_comp.toarray(out=buff).T[0])

        else:
            phi_comp = sp.coo_matrix(([],([],[])), shape=(_bus_len,1)).tocsc()
            VA_comp = None

        comp = _ContingencyCompensator(M=M, c=c, W=W, Wbar=Wbar, phi_compensator=phi_comp,\
                                        VA_compensator=VA_comp, branch_out=branch_out)

        compensators[cn] = comp

    contingency_compensators = _ContingencyCompensators(compensators=compensators, L=L_factor, U=U_factor, Pr=Pr, Pc=Pc)

    return contingency_compensators

def _check_contingencies_not_disconnecting( graph, branches, mapping_bus_to_idx, contingency_branches ):
    if len(contingency_branches) < 10:
        for bn in contingency_branches:
            connected = check_contingency_connection(graph, branches, [bn], mapping_bus_to_idx)
            if not connected:
                raise RuntimeError(f"Contingency {bn} disconnects the network!")
        return
    else:
        all_connected_contigencies = get_N_minus_1_branches(graph, branches, mapping_bus_to_idx)
        bad_contingencies = set(contingency_branches).difference(all_connected_contigencies)

        if len(bad_contingencies) == 0:
            return
        for bn in bad_contingencies:
            raise RuntimeError(f"Contingency {bn} disconnects the network!")

def _calculate_interface_matrix(interfaces, index_set_interface, mapping_branch_to_idx):
    """
    calculate an interface matrix, where a the rows correspond to each interface, and
    the columns correspond to each branch
    """
    _len_interface = len(index_set_interface)
    _len_branch = len(mapping_branch_to_idx)

    row = []
    col = []
    data = []

    for idx, i_n in enumerate(index_set_interface):
        interface = interfaces[i_n]
        for l, val in zip(interface['lines'], interface['line_orientation']):
            if val == 0:
                continue
            elif val == -1 or val == 1:
                branch_idx = mapping_branch_to_idx[l]
                row.append(idx)
                col.append(branch_idx)
                data.append(val)
            else:
                ## TODO: do we need enforce this requirement?
                raise Exception("Interface {} has line {} with line_orientation {} "
                        "not in [-1,0,1].".format(i_n, l, val))
    I = sp.coo_matrix((data,(row,col)),shape=(_len_interface,_len_branch))
    return I


def calculate_ptdf(branches,buses,index_set_branch,index_set_bus,reference_bus,base_point=BasePointType.FLATSTART,sparse_index_set_branch=None,mapping_bus_to_idx=None):
    """
    Calculates the sensitivity of the voltage angle to real power injections
    Parameters
    ----------
    branches: dict{}
        The dictionary of branches for the test case
    buses: dict{}
        The dictionary of buses for the test case
    index_set_branch: list
        The list of keys for branches for the test case
    index_set_bus: list
        The list of keys for buses for the test case
    reference_bus: key value
        The reference bus key value
    base_point: egret.model_library_defn.BasePointType
        The base-point type for calculating the PTDF matrix
    sparse_index_set_branch: list
        The list of keys for branches needed to compute a sparse PTDF matrix
        If this is None, a dense PTDF matrix is returned
    mapping_bus_to_idx: dict
        A map from bus names to indices for matrix construction. If None,
        will be inferred from index_set_bus.
    """
    _len_bus = len(index_set_bus)

    if mapping_bus_to_idx is None:
        mapping_bus_to_idx = {bus_n: i for i, bus_n in enumerate(index_set_bus)}

    _len_branch = len(index_set_branch)

    _ref_bus_idx = mapping_bus_to_idx[reference_bus]

    ## check if the network is connected
    graph = construct_connection_graph(branches, mapping_bus_to_idx)
    connected = check_network_connection(graph, index_set_bus)

    J = _calculate_J11(branches,buses,index_set_branch,index_set_bus,mapping_bus_to_idx,base_point,approximation_type=ApproximationType.PTDF)
    A = calculate_adjacency_matrix_transpose(branches,index_set_branch,index_set_bus,mapping_bus_to_idx)
    M = A@J

    if sparse_index_set_branch is None or len(sparse_index_set_branch) == _len_branch:
        ## the resulting matrix after inversion will be fairly dense,
        ## the scipy documenation recommends using dense for the inversion
        ## as well

        ref_bus_mask = np.ones(_len_bus, dtype=bool)
        ref_bus_mask[_ref_bus_idx] = False

        # M is now (A^T B_d A) with
        # row and column of reference
        # bus removed
        J0 = M[ref_bus_mask,:][:,ref_bus_mask]

        # (B_d A) with reference bus column removed
        B_dA = J[:,ref_bus_mask].A

        if connected:
            try:
                PTDF = np.linalg.solve(J0.T.A, B_dA.T).T
            except np.linalg.LinAlgError:
                logger.warning("Matrix not invertible. Calculating pseudo-inverse instead.")
                SENSI = np.linalg.pinv(J0.A,rcond=1e-7)
                PTDF = np.matmul(B_dA,SENSI)
        else:
            logger.warning("Using pseudo-inverse method as network is disconnected")
            SENSI = np.linalg.pinv(J0.A,rcond=1e-7)
            PTDF = np.matmul(B_dA,SENSI)

        # insert 0 column for reference bus
        PTDF = np.insert(PTDF, _ref_bus_idx, np.zeros(_len_branch), axis=1)

    elif len(sparse_index_set_branch) < _len_branch:
        ref_bus_row = sp.coo_matrix(([1],([0],[_ref_bus_idx])), shape=(1,_len_bus))
        ref_bus_col = sp.coo_matrix(([1],([_ref_bus_idx],[0])), shape=(_len_bus,1))
 
        J0 = sp.bmat([[M,ref_bus_col],[ref_bus_row,0]], format='coo')

        B = np.array([], dtype=np.int64).reshape(_len_bus + 1,0)
        _sparse_mapping_branch = {i: branch_n for i, branch_n in enumerate(index_set_branch) if branch_n in sparse_index_set_branch}

        ## TODO: Maybe just keep the sparse PTDFs as a dict of ndarrays?
        ## Right now the return type depends on the options 
        ## passed in
        for idx, branch_name in _sparse_mapping_branch.items():
            b = np.zeros((_len_branch,1))
            b[idx] = 1
            _tmp = J.transpose()@b
            _tmp = np.vstack([_tmp,0])
            B = np.concatenate((B,_tmp), axis=1)
        row_idx = list(_sparse_mapping_branch.keys())
        PTDF = sp.lil_matrix((_len_branch,_len_bus))
        _ptdf = sp.linalg.spsolve(J0.transpose().tocsr(), B).T
        PTDF[row_idx] = _ptdf[:,:-1]

    return PTDF


def calculate_ptdf_ldf(branches,buses,index_set_branch,index_set_bus,reference_bus,base_point=BasePointType.SOLUTION,sparse_index_set_branch=None,mapping_bus_to_idx=None):
    """
    Calculates the sensitivity of the voltage angle to real power injections and losses on the lines. Includes the
    calculation of the constant term for the quadratic losses on the lines.
    Parameters
    ----------
    branches: dict{}
        The dictionary of branches for the test case
    buses: dict{}
        The dictionary of buses for the test case
    index_set_branch: list
        The list of keys for branches for the test case
    index_set_bus: list
        The list of keys for buses for the test case
    reference_bus: key value
        The reference bus key value
    base_point: egret.model_library_defn.BasePointType
        The base-point type for calculating the PTDF and LDF matrix
    sparse_index_set_branch: list
        The list of keys for branches needed to compute a sparse PTDF matrix
    mapping_bus_to_idx: dict
        A map from bus names to indices for matrix construction. If None,
        will be inferred from index_set_bus.
    """
    _len_bus = len(index_set_bus)

    if mapping_bus_to_idx is None:
        mapping_bus_to_idx = {bus_n: i for i, bus_n in enumerate(index_set_bus)}

    _len_branch = len(index_set_branch)

    _ref_bus_idx = mapping_bus_to_idx[reference_bus]

    J = _calculate_J11(branches,buses,index_set_branch,index_set_bus,mapping_bus_to_idx,base_point,approximation_type=ApproximationType.PTDF_LOSSES)
    L = _calculate_L11(branches,buses,index_set_branch,index_set_bus,mapping_bus_to_idx,base_point)
    Jc = _calculate_pf_constant(branches,buses,index_set_branch,base_point)
    Lc = _calculate_pfl_constant(branches,buses,index_set_branch,base_point)

    if np.all(Jc == 0) and np.all(Lc == 0):
        return np.zeros((_len_branch, _len_bus)), np.zeros((_len_branch, _len_bus)), np.zeros((1,_len_branch))

    ## check if the network is connected
    graph = construct_connection_graph(branches, mapping_bus_to_idx)
    connected = check_network_connection(graph, index_set_bus)

    A = calculate_adjacency_matrix_transpose(branches,index_set_branch,index_set_bus, mapping_bus_to_idx)
    AA = calculate_absolute_adjacency_matrix(A)
    M1 = A@J
    M2 = AA@L
    M = M1 + 0.5 * M2

    ref_bus_row = sp.coo_matrix(([1],([0],[_ref_bus_idx])), shape=(1,_len_bus))
    ref_bus_col = sp.coo_matrix(([1],([_ref_bus_idx],[0])), shape=(_len_bus,1))

    J0 = sp.bmat([[M,ref_bus_col],[ref_bus_row,0]], format='coo')

    if sparse_index_set_branch is None or len(sparse_index_set_branch) == _len_branch:
        ## the resulting matrix after inversion will be fairly dense,
        ## the scipy documenation recommends using dense for the inversion
        ## as well
        if connected:
            try:
                SENSI = np.linalg.inv(J0.A)
            except np.linalg.LinAlgError:
                logger.warning("Matrix not invertible. Calculating pseudo-inverse instead.")
                SENSI = np.linalg.pinv(J0.A,rcond=1e-7)
        else:
            logger.warning("Using pseudo-inverse method as network is disconnected")
            SENSI = np.linalg.pinv(J0.A,rcond=1e-7)
        SENSI = SENSI[:-1,:-1]

        PTDF = np.matmul(J.A, SENSI)
        LDF = np.matmul(L.A, SENSI)
    elif len(sparse_index_set_branch) < _len_branch:
        B_J = np.array([], dtype=np.int64).reshape(_len_bus + 1, 0)
        B_L = np.array([], dtype=np.int64).reshape(_len_bus + 1, 0)
        _sparse_mapping_branch = {i: branch_n for i, branch_n in enumerate(index_set_branch) if branch_n in sparse_index_set_branch}

        for idx, branch_name in _sparse_mapping_branch.items():
            b = np.zeros((_len_branch, 1))
            b[idx] = 1

            _tmp_J = np.matmul(J.transpose(), b)
            _tmp_J = np.vstack([_tmp_J, 0])
            B_J = np.concatenate((B_J, _tmp_J), axis=1)

            _tmp_L = np.matmul(L.transpose(), b)
            _tmp_L = np.vstack([_tmp_L, 0])
            B_L = np.concatenate((B_L, _tmp_L), axis=1)

        row_idx = list(_sparse_mapping_branch.keys())
        PTDF = sp.lil_matrix((_len_branch, _len_bus))
        _ptdf = sp.linalg.spsolve(J0.transpose().tocsr(), B_J).T
        PTDF[row_idx] = _ptdf[:, :-1]

        LDF = sp.lil_matrix((_len_branch, _len_bus))
        _ldf = sp.linalg.spsolve(J0.transpose().tocsr(), B_L).T
        LDF[row_idx] = _ldf[:, :-1]

    M1 = A@Jc
    M2 = AA@Lc
    M = M1 + 0.5 * M2
    LDF_constant = -LDF@M + Lc

    return PTDF, LDF, LDF_constant


def calculate_interface_sensitivities(interfaces,index_set_interface,PTDFM,phase_shift_array,phi_adjust_array,mapping_branch_to_idx):
    """
    Calculates the sensitivity of interface flows to real power injections from a PTDF matrix
    Parameters
    ----------
    interfaces: dict{}
        The dictionary of interfaces for the test case
    index_set_interface : tuple
        The tuple of keys for interfaces for the test case
    PTDFM : numpy.array
        The PTDF matrix
    phase_shift_array : numpy.array
        The array of phase shifts per branch
    phi_adjust_array : numpy.array
        The array of phi adjusts per bus
    mapping_branch_to_idx: dict
        A map from branch names to indices for the PTDF matrix. If None,
        will be inferred from index_set_branch.
    """

    ## pre-allocate a matrix for the interface sensitivities
    ## size is number of interfaces (rows) by number of buses (cols)
    PTDF_I = np.zeros((len(index_set_interface),PTDFM.shape[1]))
    PTDF_I_const = np.zeros(len(index_set_interface))

    for idx, i_n in enumerate(index_set_interface):
        interface = interfaces[i_n]
        PTDF_I_row = PTDF_I[idx]
        const = 0
        for l, val in zip(interface['lines'], interface['line_orientation']):
            if val == 0:
                continue
            else:
                branch_idx = mapping_branch_to_idx[l]
                PTDF_row = PTDFM[branch_idx]
                PTDF_row_const = phase_shift_array[branch_idx] + PTDF_row.dot(phi_adjust_array)
                if val == 1:
                    PTDF_I_row += PTDF_row
                    const += PTDF_row_const
                elif val == -1:
                    PTDF_I_row -= PTDF_row
                    const -= PTDF_row_const
                else:
                    raise Exception("Interface {} has line {} with line_orientation {} "
                            "not in [-1,0,1].".format(i_n, l, val))
        PTDF_I_const[idx] = const

    return PTDF_I, PTDF_I_const


def calculate_adjacency_matrix_transpose(branches,index_set_branch,index_set_bus, mapping_bus_to_idx):
    """
    Calculates the adjacency matrix where (-1) represents flow from the bus and (1) represents flow to the bus
    for a given branch
    """
    _len_bus = len(index_set_bus)

    _len_branch = len(index_set_branch)

    row = []
    col = []
    data = []

    for idx_col, branch_name in enumerate(index_set_branch):
        branch = branches[branch_name]

        from_bus = branch['from_bus']
        row.append(mapping_bus_to_idx[from_bus])
        col.append(idx_col)
        data.append(-1)

        to_bus = branch['to_bus']
        row.append(mapping_bus_to_idx[to_bus])
        col.append(idx_col)
        data.append(1)

    adjacency_matrix = sp.coo_matrix((data,(row,col)), shape=(_len_bus, _len_branch))
    return adjacency_matrix.tocsc()


def calculate_absolute_adjacency_matrix(adjacency_matrix):
    """
    Calculates the absolute value of the adjacency matrix
    """
    return np.absolute(adjacency_matrix)

def construct_connection_graph(branches, mapping_bus_to_idx):
    _len_bus = len(mapping_bus_to_idx)

    row = []
    col = []
    data = []

    for branch in branches.values():
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        row.append(mapping_bus_to_idx[from_bus])
        col.append(mapping_bus_to_idx[to_bus])

    data = np.ones((len(branches),), dtype=np.uint8)

    graph = sp.coo_matrix((data,(row,col)), shape=(_len_bus, _len_bus)).tocsr()

    return graph

def check_contingency_connection(graph, branches, branches_removed, mapping_bus_to_idx):
    """
    Checks the connectivity after removing the branches in branches_removed

    Parameters
    ----------
    graph : output from construct_connection_graph
    branches : dictionary of branches
    branches_removed : list of names of branches removed by this contingency
    mapping_bus_to_idx : bus name to index dictionary
    """
    _len_bus = len(mapping_bus_to_idx)

    row = []
    col = []
    data = []

    for branch_name in branches_removed:
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        row.append(mapping_bus_to_idx[from_bus])
        col.append(mapping_bus_to_idx[to_bus])

    data = np.ones((len(branches_removed),), dtype=np.uint8)

    graph_delta = sp.coo_matrix((data,(row,col)), shape=(_len_bus, _len_bus)).tocsr()

    graph -= graph_delta
    n_components = sp.csgraph.connected_components(csgraph=graph, directed=False, return_labels=False)
    graph += graph_delta

    return (n_components == 1)

def get_N_minus_1_branches(graph, branches, mapping_bus_to_idx):
    """
    Gets a list of branches which can be monitored using N-1 tools

    Returns
    -------
    List of branches for which removing one does not disconnect the network
    """
    potential_bridges = set(nx.bridges(nx.Graph(graph)))

    ## networkx.bridges does not handle multi-edges cleanly, 
    ## so we need to remove some it found if a redundant
    ## branch connects two nodes
    bridges_to_remove = set()
    for bridge in potential_bridges:
        swap = (bridge[1], bridge[0])
        if graph[bridge] + graph[swap] > 1:
            bridges_to_remove.add(bridge)

    bridges = potential_bridges - bridges_to_remove

    branches_not_disconnecting = []
    for bn, branch in branches.items():
        idx_from, idx_to = mapping_bus_to_idx[branch['from_bus']], \
                            mapping_bus_to_idx[branch['to_bus']]
        if (idx_from, idx_to) in bridges:
            continue
        elif (idx_to, idx_from) in bridges:
            continue
        else:
            branches_not_disconnecting.append(bn)

    return branches_not_disconnecting 

def check_network_connection(graph, index_set_bus):
    """
    Checks for the connectivity of the network and prints some helpful information to the
    logger if the network is disconnected

    Parameters
    ----------
    graph : output from construct_connection_graph
    index_set_bus : list mapping bus indices to bus names (only used to generate warnings)
    """
    
    n_components, labels = sp.csgraph.connected_components(csgraph=graph, directed=False, return_labels=True)

    if n_components > 1:
        logger.warning("Network is disconnected. Number of components: {}".format(n_components))
        ### get the counts to eliminate the largest connected component
        unique, counts = np.unique(labels, return_counts=True)

        largest_component_label = unique[counts.argmax()]

        ## These are the indicies of the small connected components
        small_connected_components = np.nonzero(labels != largest_component_label)[0]

        components = { comp: [] for comp in unique if comp != largest_component_label }

        for idx, comp in zip(small_connected_components, labels[small_connected_components]):
            components[comp].append(index_set_bus[idx])

        logger.warning("Buses not in largest component:")
        for comp, buses in components.items():
            logger.warning("{} : {}".format(comp, buses))

    return (n_components == 1)
