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
import collections.abc as abc
import egret.model_library.transmission.branch as libbranch
import pyomo.environ as pyo
import numpy as np
import copy as cp

from enum import Enum

class LazyPTDFTerminationCondition(Enum):
    NORMAL = 1
    ITERATION_LIMIT = 2
    FLOW_VIOLATION = 3

class _LazyViolations(abc.Sized):
    def __init__(self, branch_lazy_violations,
                       interface_lazy_violations=None,
                       contingency_lazy_violations=None):
        self._branch_lazy_violations = branch_lazy_violations
        if interface_lazy_violations is None:
            self._interface_lazy_violations = set()
        else:
            self._interface_lazy_violations = interface_lazy_violations
        if contingency_lazy_violations is None:
            self._contingency_lazy_violations = set()
        else:
            self._contingency_lazy_violations = contingency_lazy_violations

    def __len__(self):
        return len(self._branch_lazy_violations) \
                + len(self._interface_lazy_violations) \
                + len(self._contingency_lazy_violations)

    @property
    def branch_lazy_violations(self):
        return self._branch_lazy_violations
    @property
    def interface_lazy_violations(self):
        return self._interface_lazy_violations
    @property
    def contingency_lazy_violations(self):
        return self._contingency_lazy_violations

class _CalculatedFlows:
    def __init__(self, PFV=None, PFV_I=None):
        self._PFV = PFV
        self._PFV_I = PFV_I

    @property
    def PFV(self):
        return self._PFV
    @property
    def PFV_I(self):
        return self._PFV_I

class _MaximalViolationsStore:
    def __init__(self, max_viol_add, mb, md, time=None):
        self.max_viol_add = max_viol_add
        self.mb = mb
        self.md = md
        self.time = time
        self.violations_store = {}
        self.total_violations = 0
        self.monitored_violations = 0

    def get_violations_named(self, name):
        for key in self.violations_store:
            if key[0] == name:
                yield key[1]

    def _min_flow_violation(self):
        if self.violations_store:
            return min(self.violations_store.values())
        else:
            return 0.

    def _min_violation_key(self):
        d = self.violations_store
        return min(d, key=d.get)

    def _add_violation(self, name, other_name, index, val):
        if other_name:
            key = ( name, (other_name, index) )
        else:
            key = ( name, index )
        self.violations_store[key] = val

        # keep the violations_store <= self.max_viol_add
        if len(self.violations_store) > self.max_viol_add:
            min_key = self._min_violation_key()
            del self.violations_store[min_key]

            if min_key == key:
                raise RuntimeError(f"Circular condition: added {new_key} to violations_store "
                        f"only to delete it. violations_store: {violations_store}")
    
    def _add_violations( self, name, other_name, viol_array, viol_indices):
        val = np.inf
        while viol_indices and val > self._min_flow_violation():

            idx = np.argmax(viol_array[viol_indices])
            val = viol_array[viol_indices[idx]]

            self._add_violation( name, other_name, viol_indices[idx], val )
            del viol_indices[idx]

    def check_and_add_violations(self, name, flow_array,
                upper_lazy_limits, upper_enforced_limits,
                lower_lazy_limits, lower_enforced_limits,
                monitored_indices, outer_name=None):

        ## check upper bound
        upper_viol_lazy_array = flow_array - upper_lazy_limits

        ## get the indices of the violation
        ## here filter by least violation in violations_store
        ## in the limit, this will become 0 eventually --
        upper_viol_lazy_idx = np.nonzero(upper_viol_lazy_array > self._min_flow_violation())[0]

        upper_viol_array = flow_array[upper_viol_lazy_idx] - upper_enforced_limits[upper_viol_lazy_idx]
        self._calculate_total_and_monitored_violations(upper_viol_array, upper_viol_lazy_idx, monitored_indices)

        ## viol_lazy_idx will hold the lines we're adding
        ## this iteration -- don't want to add lines
        ## that are already in the monitored set

        # eliminate lines in the monitored set
        upper_viol_lazy_idx = list(set(upper_viol_lazy_idx).difference(monitored_indices))

        self._add_violations( name, outer_name, upper_viol_lazy_array, upper_viol_lazy_idx )

        ## check lower bound
        lower_viol_lazy_array = lower_lazy_limits - flow_array

        ## get the indices of the violation
        ## here filter by least violation in violations_store
        ## in the limit, this will become 0 eventually --
        lower_viol_lazy_idx = np.nonzero(lower_viol_lazy_array > self._min_flow_violation())[0]

        lower_viol_array =  lower_enforced_limits[lower_viol_lazy_idx] - flow_array[lower_viol_lazy_idx]
        self._calculate_total_and_monitored_violations(lower_viol_array, lower_viol_lazy_idx, monitored_indices) 

        ## viol_lazy_idx will hold the lines we're adding
        ## this iteration -- don't want to add lines
        ## that are already in the monitored set

        # eliminate lines in the monitored set
        lower_viol_lazy_idx = list(set(lower_viol_lazy_idx).difference(monitored_indices))

        self._add_violations( name, outer_name, lower_viol_lazy_array, lower_viol_lazy_idx )


    def _calculate_total_and_monitored_violations(self, viol_array, viol_lazy_idx, monitored_indices): 
        ## viol_idx_idx will be indexed by viol_lazy_idx
        viol_idx_idx = np.nonzero(viol_array > 0)[0]
        viol_idx = frozenset(viol_lazy_idx[viol_idx_idx])

        self.total_violations += len(viol_idx)

        viol_in_mb = viol_idx.intersection(monitored_indices)
        self.monitored_violations += len(viol_in_mb)

        '''
        ## TODO: GENERALIZE CODE BELOW
        for i in viol_in_mb:
            bn = PTDF.branches_keys_masked[i]
            thermal_limit = PTDF.branch_limits_array_masked[i]
            if 'violation_penalty' in branches[bn] \
                    and branches[bn]['violation_penalty'] is not None:
                log_func = logger.info
            else:
                log_func = logger.warning
            log_func(prepend_str+_generate_flow_viol_warning(mb.pf, 'branch', bn, PFV[i], -thermal_limit, baseMVA, time))
        '''


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
    if 'branch_kv_threshold' not in ptdf_options:
        ptdf_options['branch_kv_threshold'] = None
    if 'kv_threshold_type' not in ptdf_options:
        ptdf_options['kv_threshold_type'] = 'one'
    if 'pre_lp_iteration_limit' not in ptdf_options:
        ptdf_options['pre_lp_iteration_limit'] = 100
    if 'active_flow_tol' not in ptdf_options:
        ptdf_options['active_flow_tol'] = 50.
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
                        " abs_flow_tol={0}, lazy_rel_flow_tol={1}, baseMVA={2}".format(abs_flow_tol*baseMVA, lazy_rel_flow_tol, baseMVA))

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
def add_monitored_flow_tracker(mb):
    mb._idx_monitored = list()
    mb._interfaces_monitored = list()

    # add these if there are no slacks
    # so we don't have to check later
    # for these attributes
    if not hasattr(mb, 'pf_slack_pos'):
        mb.pf_slack_pos = pyo.Var([], dense=False)
    if not hasattr(mb, 'pfi_slack_pos'):
        mb.pfi_slack_pos = pyo.Var([], dense=False)
    if not hasattr(mb, 'pfc_slack_pos'):
        mb.pfc_slack_pos = pyo.Var([], dense=False)

## violation checker
def check_violations(mb, md, PTDF, max_viol_add, time=None, prepend_str=""):

    ## PFV -- power flow vector
    ## PFV_I -- interface power flow vector
    PFV, PFV_I, _ = PTDF.calculate_masked_PFV(mb)

    ## get the lines we're monitoring
    idx_monitored = mb._idx_monitored
    interfaces_monitored = mb._interfaces_monitored

    violations_store = _MaximalViolationsStore(max_viol_add=max_viol_add, mb=mb, md=md, time=time)

    violations_store.check_and_add_violations('branch', PFV,
                                        PTDF.lazy_branch_limits, PTDF.enforced_branch_limits,
                                       -PTDF.lazy_branch_limits, -PTDF.enforced_branch_limits,
                                        idx_monitored)
    
    violations_store.check_and_add_violations('interface', PFV_I,
                                        PTDF.lazy_interface_max_limits, PTDF.enforced_interface_max_limits,
                                        PTDF.lazy_interface_min_limits, PTDF.enforced_interface_min_limits,
                                        interfaces_monitored)

    viol_lazy = _LazyViolations(branch_lazy_violations=set(violations_store.get_violations_named('branch')),
                                interface_lazy_violations=set(violations_store.get_violations_named('interface')))
    flows = _CalculatedFlows(PFV=PFV, PFV_I=PFV_I)

    return flows, violations_store.total_violations, violations_store.monitored_violations, viol_lazy

def _generate_flow_monitor_remove_message(flow_type, bn, slack, baseMVA, time):
    ret_str = "removing {0} {1} from monitored set".format(flow_type, bn)
    if time is not None:
        ret_str += " at time {}".format(time)
    ret_str += ", flow slack={0}".format(slack*baseMVA)
    return ret_str

## flow constraint remover
def remove_inactive(mb, solver, time=None, prepend_str=""):
    if time is None: # DCOPF
        model = mb
    else: # UC
        model = mb.parent_block()
    PTDF = mb._PTDF
    ptdf_options = model._ptdf_options
    baseMVA = model.model_data.data['system']['baseMVA']

    slack_tol = ptdf_options['active_flow_tol']

    persistent_solver = isinstance(solver, PersistentSolver)

    if persistent_solver:
        _load_pf_slacks(solver, model, [time])

    ## get the lines we're monitoring
    idx_monitored = mb._idx_monitored
    interfaces_monitored = mb._interfaces_monitored

    ## get the branchnname to index map
    branchname_index_map = PTDF.branchname_to_index_masked_map
    interfacename_index_map = PTDF.interfacename_to_index_map

    ## branches
    branches = model.model_data.data['elements']['branch']
    interfaces = model.model_data.data['elements']['interface']

    constr_to_remove = list()

    for bn, constr in mb.ineq_pf_branch_thermal_bounds.items():
        ## don't take out branches we were told to monitor
        if 'lazy' in branches[bn] and not branches[bn]['lazy']:
            continue
        slack = constr.slack()
        if slack_tol <= abs(slack):
            logger.debug(prepend_str+_generate_flow_monitor_remove_message('branch', bn, abs(slack), baseMVA, time))
            constr_to_remove.append(constr)
            ## remove the index from the lines we're monitoring
            idx_monitored.remove(branchname_index_map[bn])

    for i_n, constr in mb.ineq_pf_interface_bounds.items():
        ## don't take out branches we were told to monitor
        if 'lazy' in interfaces[i_n] and not interfaces[i_n]['lazy']:
            continue
        slack = constr.slack()
        if slack_tol <= abs(slack):
            logger.debug(prepend_str+_generate_flow_monitor_remove_message('interface', i_n, abs(slack), baseMVA, time))
            constr_to_remove.append(constr)
            ## remove the index from the lines we're monitoring
            interfaces_monitored.remove(interfacename_index_map[i_n])

    msg = prepend_str+"removing {} inactive transmission constraint(s)".format(len(constr_to_remove))
    if time is not None:
        msg += " at time {}".format(time)
    logger.debug(msg)

    for constr in constr_to_remove:
        if persistent_solver:
            solver.remove_constraint(constr)
        del constr
    return len(constr_to_remove)

def _check_and_generate_flow_viol_warnings(mb, md, PTDF, PFV, PFV_I, prepend_str, \
        lt_viol, gt_viol, min_viol_int, max_viol_int, time):

    ## get the lines we're monitoring
    idx_monitored = mb._idx_monitored
    interfaces_monitored = mb._interfaces_monitored

    gt_viol_in_mb = gt_viol.intersection(idx_monitored)
    lt_viol_in_mb = lt_viol.intersection(idx_monitored)

    max_viol_int_in_mb = max_viol_int.intersection(interfaces_monitored)
    min_viol_int_in_mb = min_viol_int.intersection(interfaces_monitored)

    baseMVA = md.data['system']['baseMVA']
    ## print a warning for these lines
    ## check if the found violations are in the model and print warning
    branches = md.data['elements']['branch']
    branch_hard_violations = 0
    branch_soft_violations = 0
    for i in lt_viol_in_mb:
        bn = PTDF.branches_keys_masked[i]
        thermal_limit = PTDF.branch_limits_array_masked[i]
        if 'violation_penalty' in branches[bn] \
                and branches[bn]['violation_penalty'] is not None:
            branch_soft_violations += 1
            log_func = logger.info
        else:
            branch_hard_violations += 1
            log_func = logger.warning
        log_func(prepend_str+_generate_flow_viol_warning(mb.pf, 'branch', bn, PFV[i], -thermal_limit, baseMVA, time))

    for i in gt_viol_in_mb:
        bn = PTDF.branches_keys_masked[i]
        thermal_limit = PTDF.branch_limits_array_masked[i]
        if 'violation_penalty' in branches[bn] \
                and branches[bn]['violation_penalty'] is not None:
            branch_soft_violations += 1
            log_func = logger.info
        else:
            branch_hard_violations += 1
            log_func = logger.warning
        log_func(prepend_str+_generate_flow_viol_warning(mb.pf, 'branch', bn, PFV[i], thermal_limit, baseMVA, time))

    ## break here if no interfaces
    if 'interface' not in md.data['elements']:
        return len(gt_viol_in_mb)+len(lt_viol_in_mb), 0
    ## print a warning for these interfaces if they don't have slack
    ## check if the found violations are in the model and print warning
    interfaces = md.data['elements']['interface']
    interface_hard_violations = 0
    interface_soft_violations = 0
    for i in min_viol_int_in_mb:
        i_n = PTDF.interface_keys[i]
        limit = PTDF.interface_min_limits[i]
        if 'violation_penalty' in interfaces[i_n] \
                and interfaces[i_n]['violation_penalty'] is not None:
            interface_soft_violations += 1
            log_func = logger.info
        else:
            interface_hard_violations += 1
            log_func = logger.warning
        log_func(prepend_str+_generate_flow_viol_warning(mb.pfi, 'interface', i_n, PFV_I[i], limit, baseMVA, time))

    for i in max_viol_int_in_mb:
        i_n = PTDF.interface_keys[i]
        limit = PTDF.interface_max_limits[i]
        if 'violation_penalty' in interfaces[i_n] \
                and interfaces[i_n]['violation_penalty'] is not None:
            interface_soft_violations += 1
            log_func = logger.info
        else:
            interface_hard_violations += 1
            log_func = logger.warning
        log_func(prepend_str+_generate_flow_viol_warning(mb.pfi, 'interface', i_n, PFV_I[i], limit, baseMVA, time))

    return branch_hard_violations+interface_hard_violations, branch_soft_violations+interface_soft_violations

def _generate_flow_viol_warning(expr, e_type, bn, flow, limit, baseMVA, time):
    ret_str = "WARNING: {0} {1} is in the  monitored set".format(e_type,bn)
    if time is not None:
        ret_str += " at time {}".format(time)
    ret_str += ", but flow exceeds limit!!\n\t flow={0}, limit={1}".format(flow*baseMVA, limit*baseMVA)
    ret_str += ", model_flow={}".format(pyo.value(expr[bn])*baseMVA)
    return ret_str

def _generate_flow_monitor_message(e_type, bn, flow=None, lower_limit=None, upper_limit=None, baseMVA=None, time=None):
    ret_str = "adding {0} {1} to monitored set".format(e_type, bn)
    if time is not None:
        ret_str += " at time {}".format(time)
    if flow is not None:
        ret_str += ", flow={0}, lower limit={1}, upper limit={2}".format(flow*baseMVA, lower_limit*baseMVA, upper_limit*baseMVA)
    return ret_str

## helper for generating pf
def _iter_over_viol_set(viol_set, mb, PTDF, abs_ptdf_tol, rel_ptdf_tol):
    for i in viol_set:
        bn = PTDF.branches_keys_masked[i]
        if mb.pf[bn].expr is None:
            expr = libbranch.get_power_flow_expr_ptdf_approx(mb, bn, PTDF, abs_ptdf_tol=abs_ptdf_tol, rel_ptdf_tol=rel_ptdf_tol)
            mb.pf[bn] = expr
        yield i, bn

## helper for generating pfi
def _iter_over_int_viol_set(int_viol_set, mb, PTDF, abs_ptdf_tol, rel_ptdf_tol):
    for i in int_viol_set:
        i_n = PTDF.interface_keys[i]
        if mb.pfi[i_n].expr is None:
            expr = libbranch.get_power_flow_interface_expr_ptdf(mb, i_n, PTDF, abs_ptdf_tol=abs_ptdf_tol, rel_ptdf_tol=rel_ptdf_tol)
            mb.pfi[i_n] = expr
        yield i, i_n

def _generate_branch_thermal_bounds(mb, bn, thermal_limit):
    if bn in mb.pf_slack_pos.index_set():
        if bn not in mb.pf_slack_pos:
            neg_slack = mb.pf_slack_neg[bn]
            pos_slack = mb.pf_slack_pos[bn]
            assert len(mb.pf_slack_pos) == len(mb.pf_slack_neg)
            new_var = True
        else: # the constraint could have been added and removed
            neg_slack = mb.pf_slack_neg[bn]
            pos_slack = mb.pf_slack_pos[bn]
            new_var = False
    else:
        neg_slack = None
        pos_slack = None
        new_var = False

    return libbranch.generate_thermal_bounds(mb.pf[bn], -thermal_limit, thermal_limit, neg_slack, pos_slack), new_var

def _generate_interface_bounds(mb, i_n, minimum_limit, maximum_limit):
    if i_n in mb.pfi_slack_pos.index_set():
        if i_n not in mb.pfi_slack_pos:
            neg_slack = mb.pfi_slack_neg[i_n]
            pos_slack = mb.pfi_slack_pos[i_n]
            assert len(mb.pfi_slack_pos) == len(mb.pfi_slack_neg)
            new_var = True
        else: # the constraint could have been added and removed
            neg_slack = mb.pfi_slack_neg[i_n]
            pos_slack = mb.pfi_slack_pos[i_n]
            new_var = False
    else:
        neg_slack = None
        pos_slack = None
        new_var = False

    return libbranch.generate_thermal_bounds(mb.pfi[i_n], minimum_limit, maximum_limit, neg_slack, pos_slack), new_var

## violation adder
def add_violations(lazy_violations, flows, mb, md, solver, ptdf_options,
                    PTDF, time=None, prepend_str=""):

    if time is None:
        model = mb
    else:
        model = mb.parent_block()

    baseMVA = md.data['system']['baseMVA']

    persistent_solver = isinstance(solver, PersistentSolver)

    ## static information between runs
    rel_ptdf_tol = ptdf_options['rel_ptdf_tol']
    abs_ptdf_tol = ptdf_options['abs_ptdf_tol']

    constr = mb.ineq_pf_branch_thermal_bounds
    viol_in_mb = mb._idx_monitored
    for i, bn in _iter_over_viol_set(lazy_violations.branch_lazy_violations, mb, PTDF, abs_ptdf_tol, rel_ptdf_tol):
        thermal_limit = PTDF.branch_limits_array_masked[i]
        if flows.PFV is None:
            logger.debug(prepend_str+_generate_flow_monitor_message('branch', bn, time=time))
        else:
            logger.debug(prepend_str+_generate_flow_monitor_message('branch', bn, flows.PFV[i], -thermal_limit, thermal_limit, baseMVA, time))

        constr[bn], new_slacks = _generate_branch_thermal_bounds(mb, bn, thermal_limit)
        viol_in_mb.append(i)
        if new_slacks:
            m = model
            obj_coef = m.TimePeriodLengthHours*m.BranchLimitPenalty[bn]

            if persistent_solver:
                if m is not m.model():
                    raise RuntimeError("Cannot add lazy var for branch slacks if part of a larger model")
                ## update the objective through the add_column method
                solver.add_column(m, mb.pf_slack_pos[bn], obj_coef, [], [])
                solver.add_column(m, mb.pf_slack_neg[bn], obj_coef, [], [])
            else:
                m.BranchViolationCost[time].expr += ( obj_coef*mb.pf_slack_pos[bn] + \
                                                      obj_coef*mb.pf_slack_neg[bn] )
        if persistent_solver:
            solver.add_constraint(constr[bn])

    ## in case there's no interfaces
    if not hasattr(mb, 'ineq_pf_interface_bounds'):
        return
    constr = mb.ineq_pf_interface_bounds
    int_viol_in_mb = mb._interfaces_monitored
    for i, i_n in _iter_over_int_viol_set(lazy_violations.interface_lazy_violations, mb, PTDF, abs_ptdf_tol, rel_ptdf_tol):
        minimum_limit = PTDF.interface_min_limits[i]
        maximum_limit = PTDF.interface_max_limits[i]
        if flows.PFV_I is None:
            logger.debug(prepend_str+_generate_flow_monitor_message('interface', i_n, time=time))
        else:
            logger.debug(prepend_str+_generate_flow_monitor_message('interface', i_n, flows.PFV_I[i], minimum_limit, maximum_limit, baseMVA, time))
        constr[i_n], new_slacks = _generate_interface_bounds(mb, i_n, minimum_limit, maximum_limit)
        int_viol_in_mb.append(i)
        if new_slacks:
            m = model
            obj_coef = m.TimePeriodLengthHours*m.InterfaceLimitPenalty[i_n]

            if persistent_solver:
                if m is not m.model():
                    raise RuntimeError("Cannot add lazy var for branch slacks if part of a larger model")
                ## update the objective through the add_column method
                solver.add_column(m, mb.pfi_slack_pos[i_n], obj_coef, [], [])
                solver.add_column(m, mb.pfi_slack_neg[i_n], obj_coef, [], [])
            else:
                m.InterfaceViolationCost[time].expr += (obj_coef*mb.pfi_slack_pos[i_n] + \
                                                        obj_coef*mb.pfi_slack_neg[i_n] )
        if persistent_solver:
            solver.add_constraint(constr[i_n])


## helper for generating pf
def _iter_over_initial_set(branches, branches_in_service, PTDF):
    for bn in branches_in_service:
        branch = branches[bn]
        if 'lazy' in branch and not branch['lazy']:
            if bn in PTDF.branchname_to_index_masked_map:
                i = PTDF.branchname_to_index_masked_map[bn]
                yield i, bn
            else:
                logger.warning("Branch {0} has flag 'lazy' set to False but is excluded from monitored set based on kV limits".format(bn))

### initial monitored set adder
def add_initial_monitored_branches(mb, branches, branches_in_service, ptdf_options, PTDF):
    ## static information between runs
    rel_ptdf_tol = ptdf_options['rel_ptdf_tol']
    abs_ptdf_tol = ptdf_options['abs_ptdf_tol']

    constr = mb.ineq_pf_branch_thermal_bounds
    viol_in_mb = mb._idx_monitored
    for i, bn in _iter_over_initial_set(branches, branches_in_service, PTDF):
        thermal_limit = PTDF.branch_limits_array_masked[i]
        mb.pf[bn] = libbranch.get_power_flow_expr_ptdf_approx(mb, bn, PTDF, abs_ptdf_tol=abs_ptdf_tol, rel_ptdf_tol=rel_ptdf_tol)
        constr[bn], _ = _generate_branch_thermal_bounds(mb, bn, thermal_limit)
        viol_in_mb.append(i)

def _iter_over_initial_set_interfaces(interfaces, PTDF):
    for i_n, interface in interfaces.items():
        if 'lazy' in interface and not interface['lazy']:
            i = PTDF.interfacename_to_index_map[i_n]
            yield i, i_n

### initial monitored set adder
def add_initial_monitored_interfaces(mb, interfaces, ptdf_options, PTDF):
    ## static information between runs
    rel_ptdf_tol = ptdf_options['rel_ptdf_tol']
    abs_ptdf_tol = ptdf_options['abs_ptdf_tol']

    constr = mb.ineq_pf_interface_bounds
    int_viol_in_mb = mb._interfaces_monitored
    for i, i_n in _iter_over_initial_set_interfaces(interfaces, PTDF):
        minimum_limit = PTDF.interface_min_limits[i]
        maximum_limit = PTDF.interface_max_limits[i]
        mb.pfi[i_n] = libbranch.get_power_flow_interface_expr_ptdf(mb, i_n, PTDF, abs_ptdf_tol=abs_ptdf_tol, rel_ptdf_tol=rel_ptdf_tol)
        constr[i_n], _ = _generate_interface_bounds(mb, i_n, minimum_limit, maximum_limit)
        int_viol_in_mb.append(i)

def copy_active_to_next_time(m, b_next, PTDF_next, slacks, slacks_I):
    active_slack_tol = m._ptdf_options['active_flow_tol']

    branchname_index_map = PTDF_next.branchname_to_index_masked_map
    interfacename_index_map = PTDF_next.interfacename_to_index_map

    viol_lazy = set()
    int_viol_lazy = set()

    idx_monitored = b_next._idx_monitored
    interfaces_monitored = b_next._interfaces_monitored

    for bn, slack in slacks.items():
        if abs(slack) <= active_slack_tol:
            ## in case the topology has changed
            if bn in branchname_index_map:
                idx = branchname_index_map[bn]
                if idx not in idx_monitored:
                    viol_lazy.add(idx)

    for i_n, slack in slacks_I.items():
        if abs(slack) <= active_slack_tol:
            ## in case the topology has changed
            if i_n in interfacename_index_map:
                idx = interfacename_index_map[i_n]
                if idx not in interfaces_monitored:
                    int_viol_lazy.add(idx)

    flows = _CalculatedFlows()
    viol_lazy = _LazyViolations(branch_lazy_violations=viol_lazy,
                                interface_lazy_violations=int_viol_lazy)

    return flows, viol_lazy

def _binary_var_generator(instance):
    regulation =  bool(instance.regulation_service)
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
        ivar.domain = pyo.UnitInterval
        if persistent_solver:
            for var in ivar.itervalues():
                solver.update_var(var)

def uc_instance_binary_enforcer(model, solver):
    persistent_solver = isinstance(solver, PersistentSolver)
    for ivar in _binary_var_generator(model):
        ivar.domain = pyo.Binary
        if persistent_solver:
            for var in ivar.itervalues():
                solver.update_var(var)

def _load_pf_slacks(solver, m, t_subset):
    ## ensure the slack variables are loaded
    vars_to_load = []
    for t in t_subset:
        b = m.TransmissionBlock[t]
        vars_to_load.extend(b.pf_slack_pos.values())
        vars_to_load.extend(b.pf_slack_neg.values())
        vars_to_load.extend(b.pfi_slack_pos.values())
        vars_to_load.extend(b.pfi_slack_neg.values())
    # XpressPersistent raises an exception if
    # this list is empty
    if vars_to_load:
        solver.load_vars(vars_to_load)
