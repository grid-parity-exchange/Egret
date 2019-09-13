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
from egret.model_library.defn import ApproximationType
from egret.common.log import logger
import egret.model_library.transmission.branch as libbranch
import pyomo.environ as pe
import numpy as np

from enum import Enum


class LazyPTDFTerminationCondition(Enum):
    NORMAL = 1
    ITERATION_LIMIT = 2
    FLOW_VIOLATION = 3

def populate_default_ptdf_options(ptdf_options):
    if 'rel_ptdf_tol' not in ptdf_options:
        ptdf_options['rel_ptdf_tol'] = 1.e-6
    if 'abs_ptdf_tol' not in ptdf_options:
        ptdf_options['abs_ptdf_tol'] = 1.e-10
    if 'abs_flow_tol' not in ptdf_options:
        ptdf_options['abs_flow_tol'] = 1.e-3
    if 'rel_flow_tol' not in ptdf_options:
        ptdf_options['rel_flow_tol'] = 1.e-5
    if 'iteration_limit' not in ptdf_options:
        ptdf_options['iteration_limit'] = 100000
    if 'lp_iteration_limit' not in ptdf_options:
        ptdf_options['lp_iteration_limit'] = 100
    if 'max_violations_per_iteration' not in ptdf_options:
        ptdf_options['max_violations_per_iteration'] = 1
    if 'lazy' not in ptdf_options:
        ptdf_options['lazy'] = True
    if 'load_from' not in ptdf_options:
        ptdf_options['load_from'] = None
    if 'save_to' not in ptdf_options:
        ptdf_options['save_to'] = None

def check_and_scale_ptdf_options(ptdf_options, baseMVA):
    ## scale to base MVA
    ptdf_options['abs_ptdf_tol'] /= baseMVA
    ptdf_options['abs_flow_tol'] /= baseMVA

    rel_flow_tol = ptdf_options['rel_flow_tol']
    abs_flow_tol = ptdf_options['abs_flow_tol']

    rel_ptdf_tol = ptdf_options['rel_ptdf_tol']
    abs_ptdf_tol = ptdf_options['abs_ptdf_tol']

    max_violations_per_iteration = ptdf_options['max_violations_per_iteration']

    if max_violations_per_iteration < 1 or (not isinstance(max_violations_per_iteration, int)):
        raise Exception("max_violations_per_iteration must be an integer least 1, max_violations_per_iteration={}".format(max_violations_per_iteration))

    if abs_flow_tol < 1e-6:
        logger.warning("WARNING: abs_flow_tol={0}, which is below the numeric threshold of most solvers.".format(abs_flow_tol*baseMVA))
    if abs_flow_tol < rel_ptdf_tol*10:
        logger.warning("WARNING: abs_flow_tol={0}, rel_ptdf_tol={1}, which will likely result in violations. Consider raising abs_flow_tol or lowering rel_ptdf_tol.".format(abs_flow_tol*baseMVA, rel_ptdf_tol))
    if rel_ptdf_tol < 1e-6:
        logger.warning("WARNING: rel_ptdf_tol={0}, which is low enough it may cause numerical issues in the solver. Consider rasing rel_ptdf_tol.".format(rel_ptdf_tol))
    if abs_ptdf_tol < 1e-12:
        logger.warning("WARNING: abs_ptdf_tol={0}, which is low enough it may cause numerical issues in the solver. Consider rasing abs_ptdf_tol.".format(abs_ptdf_tol*baseMVA))

## to hold the indicies of the violations
## in the model or block
def add_monitored_branch_tracker(mb):
    mb._lt_idx_monitored = list()
    mb._gt_idx_monitored = list()

## violation checker
def check_violations(mb, md, PTDF, max_viol_add, time=None):

    NWV = np.fromiter((pe.value(mb.p_nw[b]) for b in PTDF.bus_iterator()), float, count=len(PTDF.buses_keys))
    NWV += PTDF.phi_adjust_array

    PFV  = PTDF.PTDFM.dot(NWV)
    PFV += PTDF.phase_shift_array

    ## calculate the negative of the violations (for easy sorting)
    gt_viol_array = PTDF.enforced_branch_limits - PFV
    lt_viol_array = PFV + PTDF.enforced_branch_limits

    gt_viol = np.nonzero(gt_viol_array < 0)[0]
    lt_viol = np.nonzero(lt_viol_array < 0)[0]

    ## these will hold the violations 
    ## we found this iteration
    gt_viol = frozenset(gt_viol)
    lt_viol = frozenset(lt_viol)

    ## get the lines we're monitoring
    gt_idx_monitored = mb._gt_idx_monitored
    lt_idx_monitored = mb._lt_idx_monitored

    ## get the lines for which we've found a violation that's
    ## in the model
    gt_viol_in_mb = gt_viol.intersection(gt_idx_monitored)
    lt_viol_in_mb = lt_viol.intersection(lt_idx_monitored)

    ## print a warning for these lines
    ## check if the found violations are in the model and print warning
    baseMVA = md.data['system']['baseMVA']
    for i in lt_viol_in_mb:
        bn = PTDF.branches_keys[i]
        thermal_limit = PTDF.branch_limits_array[i]
        logger.warning(_generate_flow_viol_warning('LB', mb, bn, PFV[i], -thermal_limit, baseMVA, time))

    for i in gt_viol_in_mb:
        bn = PTDF.branches_keys[i]
        thermal_limit = PTDF.branch_limits_array[i]
        logger.warning(_generate_flow_viol_warning('UB', mb, bn, PFV[i], thermal_limit, baseMVA, time))

    ## *t_viol_lazy will hold the lines we're adding
    ## this iteration -- don't want to add lines
    ## that are already in the monitored set
    gt_viol_lazy = gt_viol.difference(gt_idx_monitored)
    lt_viol_lazy = lt_viol.difference(lt_idx_monitored)

    ## limit the number of lines we add in one iteration
    ## if we have too many violations, just take those largest
    ## in absolute value in either direction
    if len(gt_viol_lazy)+len(lt_viol_lazy) > max_viol_add:

        ## for those in the monitored set, assume they're feasible for
        ## the purposes of sorting the worst violations, which means
        ## resetting the values for these lines as computed above

        ## use what most solvers consider +infty
        LARGE_CONST = 1e+100
        gt_viol_array[gt_idx_monitored] = LARGE_CONST
        lt_viol_array[lt_idx_monitored] = LARGE_CONST

        ## give the order of the first max_viol_add violations
        measured_gt_viol = np.argpartition(gt_viol_array, range(max_viol_add))
        measured_lt_viol = np.argpartition(lt_viol_array, range(max_viol_add))

        measured_gt_viol_pos = 0
        measured_lt_viol_pos = 0
        gt_viol_lazy = set()
        lt_viol_lazy = set()
        for _ in range(max_viol_add):
            gt_v = gt_viol_array[measured_gt_viol[measured_gt_viol_pos]]
            lt_v = lt_viol_array[measured_lt_viol[measured_lt_viol_pos]]

            ## because we negated for sorting, this means the
            ## overall violation is more for the gt side
            ## dont have any more actual violations
            if gt_v > 0 and lt_v > 0:
                break
            elif gt_v < lt_v:
                gt_viol_lazy.add(measured_gt_viol[measured_gt_viol_pos])
                measured_gt_viol_pos += 1
            else:
                lt_viol_lazy.add(measured_lt_viol[measured_lt_viol_pos])
                measured_lt_viol_pos += 1

    viol_num = len(gt_viol)+len(lt_viol)
    monitored_viol_num = len(lt_viol_in_mb)+len(gt_viol_in_mb)

    return PFV, viol_num, monitored_viol_num, gt_viol_lazy, lt_viol_lazy
    
def _generate_flow_viol_warning(sense, mb, bn, flow, limit, baseMVA, time):
    ret_str = "WARNING: line {0} ({1}) is in the  monitored set".format(bn, sense)
    if time is not None:
        ret_str += " at time {}".format(time)
    ret_str += ", but flow exceeds limit!!\n\t flow={0}, limit={1}".format(flow*baseMVA, limit*baseMVA, sense)
    ret_str += ", model_flow={}".format(pe.value(mb.pf[bn])*baseMVA)
    return ret_str

def _generate_flow_monitor_message(sense, bn, flow, limit, baseMVA, time): 
    ret_str = "Adding line {0} ({1}) to monitored set".format(bn, sense)
    if time is not None:
        ret_str += " at time {}".format(time)
    ret_str += ", flow={0}, limit={1}".format(flow*baseMVA, limit*baseMVA)
    return ret_str

## violation adder
def add_violations(gt_viol_lazy, lt_viol_lazy, PFV, mb, md, solver, ptdf_options,
                    PTDF, time=None):

    model = mb.model()

    baseMVA = md.data['system']['baseMVA']

    persistent_solver = isinstance(solver, PersistentSolver)

    ## static information between runs
    rel_ptdf_tol = ptdf_options['rel_ptdf_tol']
    abs_ptdf_tol = ptdf_options['abs_ptdf_tol']

    ## helper for generating pf
    def _iter_over_viol_set(viol_set):
        for i in viol_set:
            bn = PTDF.branches_keys[i]
            if mb.pf[bn].expr is None:
                expr = libbranch.get_power_flow_expr_ptdf_approx(mb, bn, PTDF, abs_ptdf_tol=abs_ptdf_tol, rel_ptdf_tol=rel_ptdf_tol)
                mb.pf[bn] = expr
            yield i, bn

    constr = mb.ineq_pf_branch_thermal_lb
    lt_viol_in_mb = mb._lt_idx_monitored
    for i, bn in _iter_over_viol_set(lt_viol_lazy):
        thermal_limit = PTDF.branch_limits_array[i]
        logger.info(_generate_flow_monitor_message('LB', bn, PFV[i], -thermal_limit, baseMVA, time))
        constr[bn] = (-thermal_limit, mb.pf[bn], None)
        lt_viol_in_mb.append(i)
        if persistent_solver:
            solver.add_constraint(constr[bn])

    constr = mb.ineq_pf_branch_thermal_ub
    gt_viol_in_mb = mb._gt_idx_monitored
    for i, bn in _iter_over_viol_set(gt_viol_lazy):
        thermal_limit = PTDF.branch_limits_array[i]
        logger.info(_generate_flow_monitor_message('UB', bn, PFV[i], thermal_limit, baseMVA, time))
        constr[bn] = (None, mb.pf[bn], thermal_limit)
        gt_viol_in_mb.append(i)
        if persistent_solver:
            solver.add_constraint(constr[bn])


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

