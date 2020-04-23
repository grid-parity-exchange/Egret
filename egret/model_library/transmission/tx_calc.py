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
import numpy as np
import scipy.sparse as sp
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

        tau = 1.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']

        if approximation_type == ApproximationType.PTDF:
            x = branch['reactance']
            b = -1/(tau*x)
        elif approximation_type == ApproximationType.PTDF_LOSSES:
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
    return J11.tocsr()


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

        tau = 1.0
        shift = 0.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])

        b = 0.
        if approximation_type == ApproximationType.PTDF:
            x = branch['reactance']
            b = -(1/x)*(shift/tau)
        elif approximation_type == ApproximationType.PTDF_LOSSES:
            b = calculate_susceptance(branch)*(shift/tau)

        row_from.append(mapping_bus_to_idx[from_bus])
        row_to.append(mapping_bus_to_idx[to_bus])
        col.append(idx_col)
        data.append(b)

    phi_from = sp.coo_matrix((data,(row_from,col)), shape=(_len_bus,_len_branch))
    phi_to = sp.coo_matrix((data,(row_to,col)), shape=(_len_bus,_len_branch))

    return phi_from.tocsr(), phi_to.tocsr()


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
    connected = check_network_connection(branches, index_set_branch, index_set_bus, mapping_bus_to_idx)

    J = _calculate_J11(branches,buses,index_set_branch,index_set_bus,mapping_bus_to_idx,base_point,approximation_type=ApproximationType.PTDF)
    A = calculate_adjacency_matrix_transpose(branches,index_set_branch,index_set_bus,mapping_bus_to_idx)
    M = A@J

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
        PTDF = np.matmul(J.A,SENSI)
    elif len(sparse_index_set_branch) < _len_branch:
        B = np.array([], dtype=np.int64).reshape(_len_bus + 1,0)
        _sparse_mapping_branch = {i: branch_n for i, branch_n in enumerate(index_set_branch) if branch_n in sparse_index_set_branch}

        ## TODO: Maybe just keep the sparse PTDFs as a dict of ndarrays?
        ## Right now the return type depends on the options 
        ## passed in
        for idx, branch_name in _sparse_mapping_branch.items():
            b = np.zeros((_len_branch,1))
            b[idx] = 1
            _tmp = np.matmul(J.transpose(),b)
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
    connected = check_network_connection(branches, index_set_branch, index_set_bus, mapping_bus_to_idx)

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
    return adjacency_matrix.tocsr()


def calculate_absolute_adjacency_matrix(adjacency_matrix):
    """
    Calculates the absolute value of the adjacency matrix
    """
    return np.absolute(adjacency_matrix)

def check_network_connection(branches, index_set_branch, index_set_bus, mapping_bus_to_idx):
    """
    Checks for the connectivity of the network and prints some helpful information to the
    logger if the network is disconnected
    """
    _len_bus = len(index_set_bus)

    row = []
    col = []
    data = []

    for branch in branches.values():
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        row.append(mapping_bus_to_idx[from_bus])
        col.append(mapping_bus_to_idx[to_bus])

    data = np.ones((len(branches),), dtype=int)

    graph = sp.coo_matrix((data,(row,col)), shape=(_len_bus, _len_bus)).tocsr()

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
