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
import copy as cp

from enum import Enum


class LazyPTDFTerminationCondition(Enum):
    NORMAL = 1
    ITERATION_LIMIT = 2
    FLOW_VIOLATION = 3

def populate_default_ptdf_options(ptdf_options):
    if ptdf_options is None:
        ptdf_options = dict()
    else:
        ## get a copy
        ptdf_options = cp.deepcopy(ptdf_options)
    if 'rel_ptdf_tol' not in ptdf_options:
        ptdf_options['rel_ptdf_tol'] = 1.e-6
    if 'abs_ptdf_tol' not in ptdf_options:
        ptdf_options['abs_ptdf_tol'] = 1.e-10
    if 'abs_flow_tol' not in ptdf_options:
        ptdf_options['abs_flow_tol'] = 1.e-3
    if 'rel_flow_tol' not in ptdf_options:
        ptdf_options['rel_flow_tol'] = 1.e-5
    if 'lazy_rel_flow_tol' not in ptdf_options:
        ptdf_options['lazy_rel_flow_tol'] = -0.01
    if 'iteration_limit' not in ptdf_options:
        ptdf_options['iteration_limit'] = 100000
    if 'lp_iteration_limit' not in ptdf_options:
        ptdf_options['lp_iteration_limit'] = 100
    if 'max_violations_per_iteration' not in ptdf_options:
        ptdf_options['max_violations_per_iteration'] = 5
    if 'lazy' not in ptdf_options:
        ptdf_options['lazy'] = True
    if 'load_from' not in ptdf_options:
        ptdf_options['load_from'] = None
    if 'save_to' not in ptdf_options:
        ptdf_options['save_to'] = None
    if 'branch_kv_threshold' not in ptdf_options:
        ptdf_options['branch_kv_threshold'] = None
    if 'kv_threshold_type' not in ptdf_options:
        ptdf_options['kv_threshold_type'] = 'one'
    if 'pre_lp_iteration_limit' not in ptdf_options:
        ptdf_options['pre_lp_iteration_limit'] = 100
    if 'active_flow_tol' not in ptdf_options:
        ptdf_options['active_flow_tol'] = 10.
    if 'lp_cleanup_phase' not in ptdf_options:
        ptdf_options['lp_cleanup_phase'] = True
    return ptdf_options

def check_and_scale_ptdf_options(ptdf_options, baseMVA):
    ## scale to base MVA
    ptdf_options['abs_ptdf_tol'] /= baseMVA
    ptdf_options['abs_flow_tol'] /= baseMVA
    ptdf_options['active_flow_tol'] /= baseMVA

    ## lowercase keyword options
    ptdf_options['kv_threshold_type'] = ptdf_options['kv_threshold_type'].lower()

    rel_flow_tol = ptdf_options['rel_flow_tol']
    abs_flow_tol = ptdf_options['abs_flow_tol']

    rel_ptdf_tol = ptdf_options['rel_ptdf_tol']
    abs_ptdf_tol = ptdf_options['abs_ptdf_tol']

    lazy_rel_flow_tol = ptdf_options['lazy_rel_flow_tol']

    max_violations_per_iteration = ptdf_options['max_violations_per_iteration']

    if max_violations_per_iteration < 1 or (not isinstance(max_violations_per_iteration, int)):
        raise Exception("max_violations_per_iteration must be an integer least 1, max_violations_per_iteration={}".format(max_violations_per_iteration))

    if abs_flow_tol < lazy_rel_flow_tol:
        raise Exception("abs_flow_tol (when scaled by baseMVA) cannot be less than lazy_flow_tol"
                        " abs_flow_tol={0}, lazy_flow_tol={1}, baseMVA={2}".format(abs_flow_tol*baseMVA, lazy_flow_tol, baseMVA))

    if ptdf_options['kv_threshold_type'] not in ['one', 'both']:
        raise Exception("kv_threshold_type must be either 'one' (for at least one end of the line"
                        " above branch_kv_threshold) or 'both' (for both end of the line above"
                        " branch_kv_threshold), kv_threshold_type={}".format(ptdf_options['kv_threshold_type']))

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

def calculate_PFV(mb, PTDF):
    NWV = np.fromiter((pe.value(mb.p_nw[b]) for b in PTDF.bus_iterator()), float, count=len(PTDF.buses_keys))
    NWV += PTDF.phi_adjust_array

    PFV  = PTDF.PTDFM_masked.dot(NWV)
    PFV += PTDF.phase_shift_array_masked

    return PFV


## violation checker
def check_violations(mb, md, PTDF, max_viol_add, time=None, prepend_str=""):

    PFV = calculate_PFV(mb, PTDF)

    ## calculate the lazy violations
    gt_viol_lazy_array = PFV - PTDF.lazy_branch_limits
    lt_viol_lazy_array = -PFV - PTDF.lazy_branch_limits

    ## *_viol_lazy has the indices of the violations at
    ## the lazy limit
    gt_viol_lazy = np.nonzero(gt_viol_lazy_array > 0)[0]
    lt_viol_lazy = np.nonzero(lt_viol_lazy_array > 0)[0]

    ## calculate the violations
    ## these will be just a subset
    gt_viol_array = PFV[gt_viol_lazy] - PTDF.enforced_branch_limits[gt_viol_lazy]
    lt_viol_array = -PFV[lt_viol_lazy]- PTDF.enforced_branch_limits[lt_viol_lazy]

    ## *_viol will be indexed by *_viol_lazy
    gt_viol = np.nonzero(gt_viol_array > 0)[0]
    lt_viol = np.nonzero(lt_viol_array > 0)[0]

    ## these will hold the violations 
    ## we found this iteration
    gt_viol = frozenset(gt_viol_lazy[gt_viol])
    lt_viol = frozenset(lt_viol_lazy[lt_viol])

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
        bn = PTDF.branches_keys_masked[i]
        thermal_limit = PTDF.branch_limits_array_masked[i]
        logger.warning(prepend_str+_generate_flow_viol_warning('LB', mb, bn, PFV[i], -thermal_limit, baseMVA, time))

    for i in gt_viol_in_mb:
        bn = PTDF.branches_keys_masked[i]
        thermal_limit = PTDF.branch_limits_array_masked[i]
        logger.warning(prepend_str+_generate_flow_viol_warning('UB', mb, bn, PFV[i], thermal_limit, baseMVA, time))

    ## *t_viol_lazy will hold the lines we're adding
    ## this iteration -- don't want to add lines
    ## that are already in the monitored set

    # eliminate lines in the monitored set
    gt_viol_lazy = set(gt_viol_lazy).difference(gt_idx_monitored)
    lt_viol_lazy = set(lt_viol_lazy).difference(lt_idx_monitored)

    ## limit the number of lines we add in one iteration
    ## if we have too many violations, just take those largest
    ## in absolute value in either direction
    if len(gt_viol_lazy)+len(lt_viol_lazy) > max_viol_add:

        tracking_gt_viol_lazy = list(gt_viol_lazy)
        tracking_lt_viol_lazy = list(lt_viol_lazy)

        gt_viol_lazy = list()
        lt_viol_lazy = list()

        ## one of the tracking_*t_viol_lazy could be empty

        if not tracking_gt_viol_lazy: 
            idx = np.argmax(lt_viol_lazy_array[tracking_lt_viol_lazy])
            ptdf_idx = tracking_lt_viol_lazy.pop(idx)
            lt_viol_lazy.append(ptdf_idx)

        elif not tracking_lt_viol_lazy:
            idx = np.argmax(gt_viol_lazy_array[tracking_gt_viol_lazy])
            ptdf_idx = tracking_gt_viol_lazy.pop(idx)
            gt_viol_lazy.append(ptdf_idx)

        else: ## get the worst of both
            gt_idx = np.argmax(gt_viol_lazy_array[tracking_gt_viol_lazy])
            lt_idx = np.argmax(lt_viol_lazy_array[tracking_lt_viol_lazy])
            gt_branch_idx = tracking_gt_viol_lazy[gt_idx]
            lt_branch_idx = tracking_lt_viol_lazy[lt_idx]

            if gt_viol_lazy_array[gt_branch_idx] > lt_viol_lazy_array[lt_branch_idx]:
                ptdf_idx = gt_branch_idx
                gt_viol_lazy.append(ptdf_idx)
                del tracking_gt_viol_lazy[gt_idx]
            else:
                ptdf_idx = lt_branch_idx
                lt_viol_lazy.append(ptdf_idx)
                del tracking_lt_viol_lazy[lt_idx]

        if max_viol_add > 1:
            ptdf_lin = np.zeros(len(PTDF.buses_keys))

        ## for those in the monitored set, assume they're feasible for
        ## the purposes of sorting the worst violations, which means
        ## resetting the values for these lines as computed above
        for _ in range(max_viol_add-1):

            ptdf_lin += PTDF.PTDFM_masked[ptdf_idx]

            all_other_violations = list(tracking_gt_viol_lazy + tracking_lt_viol_lazy)

            other_gt_viols = gt_viol_lazy_array[all_other_violations]
            other_lt_viols = lt_viol_lazy_array[all_other_violations]

            other_viols = np.maximum(other_gt_viols, other_lt_viols)

            ## put this in baseMVA
            other_viols *= baseMVA

            other_viol_rows = PTDF.PTDFM_masked[all_other_violations]

            orthogonality = np.absolute(np.dot(other_viol_rows, ptdf_lin))

            ## divide by transmission limits to give higher
            ## priority to those lines with larger violations

            ## larger values emphasize violation
            ## smaller emphasize orthogonality
            ## TODO: try weighting by number of nonzeros
            orthogonality /= other_viols

            ## this is the index into the orthogonality matrix,
            ## which is indexed by all_other_violations
            all_other_idx = np.argmin(orthogonality)

            ptdf_idx = all_other_violations[all_other_idx]

            if ptdf_idx in tracking_gt_viol_lazy:
                tracking_gt_viol_lazy.remove(ptdf_idx)
                gt_viol_lazy.append(ptdf_idx)
            elif ptdf_idx in tracking_lt_viol_lazy:
                tracking_lt_viol_lazy.remove(ptdf_idx)
                lt_viol_lazy.append(ptdf_idx)
            else:
                raise Exception("Unexpected case")


    viol_num = len(gt_viol)+len(lt_viol)
    monitored_viol_num = len(lt_viol_in_mb)+len(gt_viol_in_mb)

    return PFV, viol_num, monitored_viol_num, gt_viol_lazy, lt_viol_lazy


def _generate_branch_remove_message(sense, bn, slack, baseMVA, time):
    ret_str = "removing line {0} ({1}) from monitored set".format(bn, sense)
    if time is not None:
        ret_str += " at time {}".format(time)
    ret_str += ", flow slack={0}".format(slack*baseMVA)
    return ret_str

## flow constraint remover
def remove_inactive(mb, solver, time=None, prepend_str=""):
    model = mb.model()
    PTDF = mb._PTDF
    ptdf_options = model._ptdf_options
    baseMVA = model.model_data.data['system']['baseMVA']

    slack_tol = ptdf_options['active_flow_tol']

    persistent_solver = isinstance(solver, PersistentSolver)

    ## get the lines we're monitoring
    gt_idx_monitored = mb._gt_idx_monitored
    lt_idx_monitored = mb._lt_idx_monitored

    ## get the branchnname to index map
    branchname_index_map = PTDF.branchname_to_index_masked_map

    constr_to_remove = list()

    for bn, constr in mb.ineq_pf_branch_thermal_lb.items():
        slack = constr.slack()
        if slack_tol <= abs(slack):
            logger.debug(prepend_str+_generate_branch_remove_message('LB', bn, abs(slack), baseMVA, time))
            constr_to_remove.append(constr)
            ## remove the index from the lines we're monitoring
            lt_idx_monitored.remove(branchname_index_map[bn])

    for bn, constr in mb.ineq_pf_branch_thermal_ub.items():
        slack = constr.slack()
        if slack_tol <= abs(slack):
            logger.debug(prepend_str+_generate_branch_remove_message('UB', bn, abs(slack), baseMVA, time))
            constr_to_remove.append(constr)
            ## remove the index from the lines we're monitoring
            gt_idx_monitored.remove(branchname_index_map[bn])

    msg = prepend_str+"removing {} inactive transmission constraint(s)".format(len(constr_to_remove))
    if time is not None:
        msg += " at time {}".format(time)
    logger.debug(msg)

    for constr in constr_to_remove:
        if persistent_solver:
            solver.remove_constraint(constr)
        del constr
    return len(constr_to_remove)


def _generate_flow_viol_warning(sense, mb, bn, flow, limit, baseMVA, time):
    ret_str = "WARNING: line {0} ({1}) is in the  monitored set".format(bn, sense)
    if time is not None:
        ret_str += " at time {}".format(time)
    ret_str += ", but flow exceeds limit!!\n\t flow={0}, limit={1}".format(flow*baseMVA, limit*baseMVA, sense)
    ret_str += ", model_flow={}".format(pe.value(mb.pf[bn])*baseMVA)
    return ret_str

def _generate_flow_monitor_message(sense, bn, flow=None, limit=None, baseMVA=None, time=None):
    ret_str = "adding line {0} ({1}) to monitored set".format(bn, sense)
    if time is not None:
        ret_str += " at time {}".format(time)
    if flow is not None:
        ret_str += ", flow={0}, limit={1}".format(flow*baseMVA, limit*baseMVA)
    return ret_str

## violation adder
def add_violations(gt_viol_lazy, lt_viol_lazy, PFV, mb, md, solver, ptdf_options,
                    PTDF, time=None, prepend_str=""):

    model = mb.model()

    baseMVA = md.data['system']['baseMVA']

    persistent_solver = isinstance(solver, PersistentSolver)

    ## static information between runs
    rel_ptdf_tol = ptdf_options['rel_ptdf_tol']
    abs_ptdf_tol = ptdf_options['abs_ptdf_tol']

    ## helper for generating pf
    def _iter_over_viol_set(viol_set):
        for i in viol_set:
            bn = PTDF.branches_keys_masked[i]
            if mb.pf[bn].expr is None:
                expr = libbranch.get_power_flow_expr_ptdf_approx(mb, bn, PTDF, abs_ptdf_tol=abs_ptdf_tol, rel_ptdf_tol=rel_ptdf_tol)
                mb.pf[bn] = expr
            yield i, bn

    constr = mb.ineq_pf_branch_thermal_lb
    lt_viol_in_mb = mb._lt_idx_monitored
    for i, bn in _iter_over_viol_set(lt_viol_lazy):
        thermal_limit = PTDF.branch_limits_array_masked[i]
        if PFV is None:
            logger.debug(prepend_str+_generate_flow_monitor_message('LB', bn, time=time))
        else:
            logger.debug(prepend_str+_generate_flow_monitor_message('LB', bn, PFV[i], -thermal_limit, baseMVA, time))
        constr[bn] = (-thermal_limit, mb.pf[bn], None)
        lt_viol_in_mb.append(i)
        if persistent_solver:
            solver.add_constraint(constr[bn])

    constr = mb.ineq_pf_branch_thermal_ub
    gt_viol_in_mb = mb._gt_idx_monitored
    for i, bn in _iter_over_viol_set(gt_viol_lazy):
        thermal_limit = PTDF.branch_limits_array_masked[i]
        if PFV is None:
            logger.debug(prepend_str+_generate_flow_monitor_message('UB', bn, time=time))
        else:
            logger.debug(prepend_str+_generate_flow_monitor_message('UB', bn, PFV[i], thermal_limit, baseMVA, time))
        constr[bn] = (None, mb.pf[bn], thermal_limit)
        gt_viol_in_mb.append(i)
        if persistent_solver:
            solver.add_constraint(constr[bn])


def copy_active_to_next_time(m, b_next, PTDF_next, slacks_ub, slacks_lb):
    active_slack_tol = m._ptdf_options['active_flow_tol']

    branchname_index_map= PTDF_next.branchname_to_index_masked_map

    lt_viol_lazy = set()
    gt_viol_lazy = set()

    for (bn, constr), slack in slacks_lb.items():
        if abs(slack) <= active_slack_tol:
            ## in case the topology has changed
            if bn in branchname_index_map:
                lt_viol_lazy.add(branchname_index_map[bn])
    for (bn, constr), slack in slacks_ub.items():
        if abs(slack) <= active_slack_tol:
            ## in case the topology has changed
            if bn in branchname_index_map:
                gt_viol_lazy.add(branchname_index_map[bn])

    return None, gt_viol_lazy, lt_viol_lazy


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

