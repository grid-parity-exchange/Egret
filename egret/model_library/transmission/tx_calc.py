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
from math import cos, sin
from egret.model_library.defn import BasePointType, ApproximationType

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


def _calculate_J11(branches,buses,index_set_branch,index_set_bus,base_point=BasePointType.FLATSTART,approximation_type=ApproximationType.PTDF):
    """
    Compute the power flow Jacobian for partial derivative of real power flow to voltage angle
    """
    _len_bus = len(index_set_bus)
    _mapping_bus = {i: index_set_bus[i] for i in list(range(0,_len_bus))}

    _len_branch = len(index_set_branch)
    _mapping_branch = {i: index_set_branch[i] for i in list(range(0,_len_branch))}

    J11 = np.zeros((_len_branch,_len_bus))

    for idx_row, branch_name in _mapping_branch.items():
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        tau = 1.0
        shift = 0.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])

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
        idx_col = [key for key, value in _mapping_bus.items() if value == from_bus][0]
        J11[idx_row][idx_col] = -b * vn * vm * cos(tn - tm - shift)

        idx_col = [key for key, value in _mapping_bus.items() if value == to_bus][0]
        J11[idx_row][idx_col] = b * vn * vm * cos(tn - tm - shift)

    return J11


def _calculate_L11(branches,buses,index_set_branch,index_set_bus,base_point=BasePointType.FLATSTART):
    """
    Compute the power flow Jacobian for partial derivative of real power losses to voltage angle
    """
    _len_bus = len(index_set_bus)
    _mapping_bus = {i: index_set_bus[i] for i in list(range(0,_len_bus))}

    _len_branch = len(index_set_branch)
    _mapping_branch = {i: index_set_branch[i] for i in list(range(0,_len_branch))}

    L11 = np.zeros((_len_branch,_len_bus))

    for idx_row, branch_name in _mapping_branch.items():
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        tau = 1.0
        shift = 0.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])
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

        idx_col = [key for key, value in _mapping_bus.items() if value == from_bus][0]
        L11[idx_row][idx_col] = 2 * g * vn * vm * sin(tn - tm - shift)

        idx_col = [key for key, value in _mapping_bus.items() if value == to_bus][0]
        L11[idx_row][idx_col] = -2 * g * vn * vm * sin(tn - tm - shift)

    return L11


def _calculate_pf_constant(branches,buses,index_set_branch,base_point=BasePointType.FLATSTART):
    """
    Compute the power flow constant for the taylor series expansion of real power flow as
    a convex combination of the from/to directions, i.e.,
    pf = 0.5*g*((tau*vn)^2 - vm^2) - tau*vn*vm*b*sin(tn-tm-shift)
    """

    _len_branch = len(index_set_branch)
    _mapping_branch = {i: index_set_branch[i] for i in list(range(0,_len_branch))}

    pf_constant = np.zeros(_len_branch)

    for idx_row, branch_name in _mapping_branch.items():
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
                               - b * vn * vm * (sin(tn - tm - shift) - cos(tn - tm - shift)*(tn - tm))

    return pf_constant


def _calculate_pfl_constant(branches,buses,index_set_branch,base_point=BasePointType.FLATSTART):
    """
    Compute the power losses constant for the taylor series expansion of real power losses as
    a convex combination of the from/to directions, i.e.,
    pfl = g*((tau*vn)^2 + vm^2) - 2*tau*vn*vm*g*cos(tn-tm-shift)
    """

    _len_branch = len(index_set_branch)
    _mapping_branch = {i: index_set_branch[i] for i in list(range(0,_len_branch))}

    pfl_constant = np.zeros(_len_branch)

    for idx_row, branch_name in _mapping_branch.items():
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
                              - 2 * g * vn * vm * (sin(tn - tm - shift) * (tn - tm) + cos(tn - tm - shift))

    return pfl_constant

def calculate_ptdf(branches,buses,index_set_branch,index_set_bus,reference_bus,base_point=BasePointType.FLATSTART):
    """
    Calculates the sensitivity of the voltage angle to real power injections
    """

    _len_bus = len(index_set_bus)
    _mapping_bus = {i: index_set_bus[i] for i in list(range(0,_len_bus))}

    _len_branch = len(index_set_branch)
    _mapping_branch = {i: index_set_branch[i] for i in list(range(0,_len_branch))}

    _ref_bus_idx = [key for key, value in _mapping_bus.items() if value == reference_bus][0]

    J = _calculate_J11(branches,buses,index_set_branch,index_set_bus,base_point,approximation_type=ApproximationType.PTDF)
    A = calculate_adjacency_matrix(branches,index_set_branch,index_set_bus)
    M = np.matmul(A.transpose(),J)

    J0 = np.zeros((_len_bus+1,_len_bus+1))
    J0[:-1,:-1] = M
    J0[-1][_ref_bus_idx] = 1
    J0[_ref_bus_idx][-1] = 1

    try:
        SENSI = np.linalg.inv(J0)
    except np.linalg.LinAlgError:
        print("Matrix not invertible. Calculating pseudo-inverse instead.")
        SENSI = np.linalg.pinv(J0,rcond=1e-7)
        pass
    SENSI = SENSI[:-1,:-1]

    PTDF = np.matmul(J,SENSI)

    return PTDF


def calculate_ptdf_ldf(branches,buses,index_set_branch,index_set_bus,reference_bus,base_point=BasePointType.SOLUTION):
    """
    Calculates the sensitivity of the voltage angle to real power injections and losses on the lines. Includes the
    calculation of the constant term for the quadratic losses on the lines.
    """

    _len_bus = len(index_set_bus)
    _mapping_bus = {i: index_set_bus[i] for i in list(range(0,_len_bus))}

    _len_branch = len(index_set_branch)
    _mapping_branch = {i: index_set_branch[i] for i in list(range(0,_len_branch))}

    _ref_bus_idx = [key for key, value in _mapping_bus.items() if value == reference_bus][0]

    J = _calculate_J11(branches,buses,index_set_branch,index_set_bus,base_point,approximation_type=ApproximationType.PTDF_LOSSES)
    L = _calculate_L11(branches,buses,index_set_branch,index_set_bus,base_point)
    Jc = _calculate_pf_constant(branches,buses,index_set_branch,base_point)
    Lc = _calculate_pfl_constant(branches,buses,index_set_branch,base_point)

    A = calculate_adjacency_matrix(branches,index_set_branch,index_set_bus)
    AA = calculate_absolute_adjacency_matrix(A)
    M1 = np.matmul(A.transpose(),J)
    M2 = np.matmul(AA.transpose(),L)
    M = M1 + 0.5 * M2

    J0 = np.zeros((_len_bus+1,_len_bus+1))
    J0[:-1,:-1] = M
    J0[-1][_ref_bus_idx] = 1
    J0[_ref_bus_idx][-1] = 1

    try:
        SENSI = np.linalg.inv(J0)
    except np.linalg.LinAlgError:
        print("Matrix not invertible. Calculating pseudo-inverse instead.")
        SENSI = np.linalg.pinv(J0,rcond=1e-7)
        pass
    SENSI = SENSI[:-1,:-1]

    PTDF = np.matmul(J, SENSI)
    LDF = np.matmul(L,SENSI)

    M1 = np.matmul(A.transpose(), Jc)
    M2 = np.matmul(AA.transpose(), Lc)
    M = M1 + 0.5 * M2
    LDF_constant = -np.matmul(LDF,M) + Lc

    return PTDF, LDF, LDF_constant


def calculate_adjacency_matrix(branches,index_set_branch,index_set_bus):
    """
    Calculates the adjacency matrix where (-1) represents flow from the bus and (1) represents flow to the bus
    for a given branch
    """
    _len_bus = len(index_set_bus)
    _mapping_bus = {i: index_set_bus[i] for i in list(range(0,_len_bus))}

    _len_branch = len(index_set_branch)
    _mapping_branch = {i: index_set_branch[i] for i in list(range(0,_len_branch))}

    adjacency_matrix = np.zeros((_len_branch,_len_bus))

    for idx_row, branch_name in _mapping_branch.items():
        branch = branches[branch_name]

        from_bus = branch['from_bus']
        idx_col = [key for key, value in _mapping_bus.items() if value == from_bus][0]
        adjacency_matrix[idx_row,idx_col] = -1

        to_bus = branch['to_bus']
        idx_col = [key for key, value in _mapping_bus.items() if value == to_bus][0]
        adjacency_matrix[idx_row,idx_col] = 1

    return adjacency_matrix


def calculate_absolute_adjacency_matrix(adjacency_matrix):
    """
    Calculates the absolute value of the adjacency matrix
    """
    return np.absolute(adjacency_matrix)
