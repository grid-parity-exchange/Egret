#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## helpers for flow verification across dcopf and unit commitment models
from pyomo.solvers.plugins.solvers.persistent_solver import PersistentSolver
import egret.model_library.transmission.branch as libbranch
from egret.model_library.defn import ApproximationType
import pyomo.environ as pe
import numpy as np

from enum import Enum

class LazyPTDFTerminationCondition(Enum):
    NORMAL = 1
    ITERATION_LIMIT = 2
    FLOW_VIOLATION = 3

def populate_default_ptdf_options(ptdf_options_dict):
    if 'rel_ptdf_tol' not in ptdf_options_dict:
        ptdf_options_dict['rel_ptdf_tol'] = 1.e-6
    if 'abs_ptdf_tol' not in ptdf_options_dict:
        ptdf_options_dict['abs_ptdf_tol'] = 1.e-10
    if 'abs_flow_tol' not in ptdf_options_dict:
        ptdf_options_dict['abs_flow_tol'] = 1.e-3
    if 'rel_flow_tol' not in ptdf_options_dict:
        ptdf_options_dict['rel_flow_tol'] = 1.e-5
    if 'lazy_rel_flow_tol' not in ptdf_options_dict:
        ptdf_options_dict['lazy_rel_flow_tol'] = -0.15
    if 'iteration_limit' not in ptdf_options_dict:
        ptdf_options_dict['iteration_limit'] = 100000
    if 'lp_iteration_limit' not in ptdf_options_dict:
        ptdf_options_dict['lp_iteration_limit'] = 100

def check_and_scale_ptdf_options(ptdf_options_dict, baseMVA):
    ## scale to base MVA
    ptdf_options_dict['abs_ptdf_tol'] /= baseMVA
    ptdf_options_dict['abs_flow_tol'] /= baseMVA

    rel_flow_tol = ptdf_options_dict['rel_flow_tol']
    abs_flow_tol = ptdf_options_dict['abs_flow_tol']

    rel_ptdf_tol = ptdf_options_dict['rel_ptdf_tol']
    abs_ptdf_tol = ptdf_options_dict['abs_ptdf_tol']

    lazy_rel_flow_tol = ptdf_options_dict['lazy_rel_flow_tol']

    if abs_flow_tol < lazy_rel_flow_tol:
        raise Exception("abs_flow_tol (when scaled by baseMVA) cannot be less than lazy_flow_tol"
                        " abs_flow_tol={0}, lazy_flow_tol={1}, baseMVA={2}".format(abs_flow_tol*baseMVA, lazy_flow_tol, baseMVA))
    if abs_flow_tol < 1e-6:
        print("WARNING: abs_flow_tol={0}, which is below the numeric threshold of most solvers.".format(abs_flow_tol*baseMVA))
    if abs_flow_tol < rel_ptdf_tol*10:
        print("WARNING: abs_flow_tol={0}, rel_ptdf_tol={1}, which will likely result in violations. Consider raising abs_flow_tol or lowering rel_ptdf_tol.".format(abs_flow_tol*baseMVA, rel_ptdf_tol))
    if rel_ptdf_tol < 1e-6:
        print("WARNING: rel_ptdf_tol={0}, which is low enough it may cause numerical issues in the solver. Consider rasing rel_ptdf_tol.".format(rel_ptdf_tol))
    if abs_ptdf_tol < 1e-12:
        print("WARNING: abs_ptdf_tol={0}, which is low enough it may cause numerical issues in the solver. Consider rasing abs_ptdf_tol.".format(abs_ptdf_tol*baseMVA))

## violation checker
def check_violations(PTDF_dict, bus_nw_exprs):

    PTDFM = PTDF_dict['PTDFM']
    enforced_branch_limits = PTDF_dict['enforced_branch_limits']
    lazy_branch_limits = PTDF_dict['lazy_branch_limits']

    NWV = np.array([pe.value(bus_nw_expr) for bus_nw_expr in bus_nw_exprs])

    PFV  = np.dot(PTDFM, NWV)

    ## get the indices of the violations, but do it in numpy
    gt_viol_lazy = np.nonzero(np.greater(PFV, lazy_branch_limits))[0]
    lt_viol_lazy = np.nonzero(np.less(PFV, -lazy_branch_limits))[0]

    gt_viol = np.nonzero(np.greater(PFV, enforced_branch_limits))[0]
    lt_viol = np.nonzero(np.less(PFV, -enforced_branch_limits))[0]

    viol_num = len(gt_viol)+len(lt_viol)

    return PFV, viol_num, (gt_viol, lt_viol, gt_viol_lazy, lt_viol_lazy)
    
def _generate_flow_viol_warning(sense, bn, flow, limit, baseMVA, time):
    ret_str = "WARNING: line {0} ({1}) is in the  monitored set".format(bn, sense)
    if time is not None:
        ret_str += " at time {}".format(time)
    ret_str += ", but flow exceeds limit!!\n\t flow={0}, limit={1}".format(flow*baseMVA, limit*baseMVA, sense)
    return ret_str

def _generate_flow_monitor_message(sense, bn, flow, limit, baseMVA, time): 
    ret_str = "Adding line {0} ({1}) to monitored set".format(bn, sense)
    if time is not None:
        ret_str += " at time {}".format(time)
    ret_str += ", flow={0}, limit={1}".format(flow*baseMVA, limit*baseMVA)
    return ret_str

## violation adder
def add_violations(viols_tup, PFV, mb, md, solver, ptdf_options_dict,
                    PTDF_dict, bus_nw_exprs, bus_p_loads,
                    time=None):

    model = mb.model()

    persistent_solver = isinstance(solver, PersistentSolver)
    baseMVA = md.data['system']['baseMVA']

    gt_viol, lt_viol, gt_viol_lazy, lt_viol_lazy = viols_tup
    ## static information between runs
    rel_ptdf_tol = ptdf_options_dict['rel_ptdf_tol']
    abs_ptdf_tol = ptdf_options_dict['abs_ptdf_tol']

    ## get needed data
    PTDFM = PTDF_dict['PTDFM']
    buses_idx = PTDF_dict['buses_idx']
    branches_idx = PTDF_dict['branches_idx']
    branch_limits = PTDF_dict['branch_limits']
    branches = model._branches

    ## helper for generating pf
    def _iter_over_viol_set(viol_set):
        for i in viol_set:
            bn = branches_idx[i]
            branch = branches[bn]
            if mb.pf[bn].expr is None:
                if 'ptdf' not in branch:
                    branch['ptdf'] = {bus : PTDFM[i,j] for j, bus in enumerate(buses_idx)}
                expr = libbranch.get_power_flow_expr_ptdf_approx_from_nwe(mb, branch, buses_idx, bus_nw_exprs, abs_ptdf_tol=abs_ptdf_tol, rel_ptdf_tol=rel_ptdf_tol, approximation_type=ApproximationType.PTDF)
                mb.pf[bn] = expr
            yield i, bn, branch

    lt_viol_in_constr = 0
    for i, bn, branch in _iter_over_viol_set(lt_viol_lazy):
        constr = mb.ineq_pf_branch_thermal_lb
        thermal_limit = branch['rating_long_term']
        if bn in constr and i in lt_viol:
            print(_generate_flow_viol_warning('LB', bn, PFV[i], -thermal_limit, baseMVA, time))
            lt_viol_in_constr += 1
        elif bn not in constr: 
            print(_generate_flow_monitor_message('LB', bn, PFV[i], -thermal_limit, baseMVA, time))
            constr[bn] = (-thermal_limit, mb.pf[bn], None)
            if persistent_solver:
                solver.add_constraint(constr[bn])

    gt_viol_in_constr = 0
    for i, bn, branch in _iter_over_viol_set(gt_viol_lazy):
        constr = mb.ineq_pf_branch_thermal_ub
        thermal_limit = branch['rating_long_term']
        if bn in constr and i in gt_viol:
            print(_generate_flow_viol_warning('UB', bn, PFV[i], thermal_limit, baseMVA, time))
            gt_viol_in_constr += 1
        elif bn not in constr:
            print(_generate_flow_monitor_message('UB', bn, PFV[i], thermal_limit, baseMVA, time))
            constr[bn] = (None, mb.pf[bn], thermal_limit)
            if persistent_solver:
                solver.add_constraint(constr[bn])

    all_viol_in_mb = (len(lt_viol) == lt_viol_in_constr) \
                      and (len(gt_viol) == gt_viol_in_constr)


def _binary_var_generator(instance):
    regulation =  hasattr(instance, 'regulation_service')
    if instance.status_vars in ['CA_1bin_vars', 'garver_3bin_vars', 'garver_2bin_vars', 'garver_3bin_relaxed_stop_vars']:
        yield instance.UnitOn
    if instance.status_vars in ['ALS_state_transition_vars']:
        yield instance.UnitStayOn
    if instance.status_vars in ['garver_3bin_vars', 'garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars']:
        yield instance.UnitStart
    if instance.status_vars in ['garver_3bin_vars', 'ALS_state_transition_vars']:
        yield instance.UnitStop
    if regulation:
        yield instance.RegulationOn

    yield instance.OutputStorage
    yield instance.InputStorage

    if instance.startup_costs in ['KOW_startup_costs']:
        yield instance.StartupIndicator
    elif instance.startup_costs in ['MLR_startup_costs', 'MLR_startup_costs2',]:
        yield instance.delta

def uc_instance_binary_relaxer(model, solver):
    persistent_solver = isinstance(solver, PersistentSolver)
    for ivar in _binary_var_generator(model):
        ivar.domain = pe.UnitInterval
        if persistent_solver:
            for var in ivar.itervalues():
                solver.update_var(var)

def uc_instance_binary_enforcer(model, solver):
    persistent_solver = isinstance(solver, PersistentSolver)
    for ivar in _binary_var_generator(model):
        ivar.domain = pe.Binary
        if persistent_solver:
            for var in ivar.itervalues():
                solver.update_var(var)

