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
from egret.common.log import logger, logging
import collections.abc as abc
import egret.model_library.transmission.branch as libbranch
import pyomo.environ as pyo
import numpy as np
import copy as cp
import math

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
    def __init__(self, max_viol_add, md, prepend_str, time=None):
        self.max_viol_add = max_viol_add
        self.baseMVA = md.data['system']['baseMVA']
        self.time = time
        self.prepend_str = prepend_str
        self.violations_store = {}
        self.total_violations = 0
        self.monitored_violations = 0

    def get_violations_named(self, name):
        for key in self.violations_store:
            if key[0] == name:
                yield key[1]

    def min_flow_violation(self):
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

            if min_key == key:
                raise RuntimeError(f"Circular condition: added {key} to violations_store "
                        f"with value {self.violations_store[min_key]} only to delete it. "
                        f"violations_store: {self.violations_store}")

            del self.violations_store[min_key]
    
    def _add_violations( self, name, other_name, viol_array, viol_indices):
        while viol_indices:
            idx = np.argmax(viol_array[viol_indices])
            val = viol_array[viol_indices[idx]]
            if val < self.min_flow_violation() and len(self.violations_store) >= self.max_viol_add:
                break
            # If this violation is close in value to
            # one already in the set, it is likely
            # to be a parallel constraint.
            # If we haven't added any constraints yet
            # any(()) is False, so this won't fire
            close_to_existing = any( math.isclose( val, existing ) for existing in self.violations_store.values() )
            if close_to_existing:
                viol_indices.pop(idx)
                continue
            self._add_violation( name, other_name, viol_indices[idx], val )
            viol_indices.pop(idx)

    def check_and_add_violations(self, name, flow_array, flow_variable,
                upper_lazy_limits, upper_enforced_limits,
                lower_lazy_limits, lower_enforced_limits,
                monitored_indices, index_names, outer_name=None, PFV=None):

        if outer_name:
            # contingencies are named by cn, branch_idx, reduce to
            # branch_idx for this function
            monitored_indices = set(idx[1] for idx in monitored_indices if idx[0] == outer_name)

        ## check upper bound
        upper_viol_lazy_array = flow_array - upper_lazy_limits

        ## get the indices of the violation
        ## here filter by least violation in violations_store
        ## in the limit, this will become 0 eventually --
        upper_viol_lazy_idx = np.nonzero(upper_viol_lazy_array > self.min_flow_violation())[0]

        upper_viol_array = flow_array[upper_viol_lazy_idx] - upper_enforced_limits[upper_viol_lazy_idx]
        self._calculate_total_and_monitored_violations(upper_viol_array, upper_viol_lazy_idx, monitored_indices,
                                                        flow_variable, flow_array, index_names, upper_enforced_limits,
                                                        name, outer_name, PFV)

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
        lower_viol_lazy_idx = np.nonzero(lower_viol_lazy_array > self.min_flow_violation())[0]

        lower_viol_array =  lower_enforced_limits[lower_viol_lazy_idx] - flow_array[lower_viol_lazy_idx]
        self._calculate_total_and_monitored_violations(lower_viol_array, lower_viol_lazy_idx, monitored_indices,
                                                        flow_variable, flow_array, index_names, lower_enforced_limits,
                                                        name, outer_name, PFV)

        ## viol_lazy_idx will hold the lines we're adding
        ## this iteration -- don't want to add lines
        ## that are already in the monitored set

        # eliminate lines in the monitored set
        lower_viol_lazy_idx = list(set(lower_viol_lazy_idx).difference(monitored_indices))

        self._add_violations( name, outer_name, lower_viol_lazy_array, lower_viol_lazy_idx )


    def _calculate_total_and_monitored_violations(self, viol_array, viol_lazy_idx, monitored_indices,
                                                        flow_variable, flow_array, index_names, limits,
                                                        name, outer_name, other_flows ):
        ## viol_idx_idx will be indexed by viol_lazy_idx
        viol_idx_idx = np.nonzero(viol_array > 0)[0]
        viol_idx = frozenset(viol_lazy_idx[viol_idx_idx])

        self.total_violations += len(viol_idx)

        viol_in_mb = viol_idx.intersection(monitored_indices)
        self.monitored_violations += len(viol_in_mb)

        for i in viol_in_mb:
            element_name = index_names[i]
            thermal_limit = limits[i]
            flow = flow_array[i]
            if outer_name:
                element_name = (outer_name, element_name)
                thermal_limit += other_flows[i]
                flow += other_flows[i]
            logger.info(self.prepend_str+_generate_flow_viol_warning(flow_variable, name, element_name, flow, thermal_limit, self.baseMVA, self.time))

        ## useful debugging code
        if logger.level <= logging.DEBUG:
            for i in monitored_indices:
                element_name = index_names[i]
                thermal_limit = limits[i]
                flow = flow_array[i]
                if outer_name:
                    element_name = (outer_name, element_name)
                    thermal_limit += other_flows[i]
                    flow += other_flows[i]
                    print(f'contingency: {element_name[0]}, branch: {element_name[1]}')
                    print(f'delta: {flow_array[i]}')
                    print(f'base : {other_flows[i]}')
                    print(f'flow : {flow_array[i]+other_flows[i]}')
                    print(f'model: {pyo.value(flow_variable[element_name])}')
                    if not math.isclose(pyo.value(flow_variable[element_name]), flow_array[i]+other_flows[i]):
                        print(f'contingency: {element_name[0]}, branch_idx: {i}')
                        diff = pyo.value(flow_variable[element_name]) - (flow_array[i]+other_flows[i])
                        print(f'ABSOLUTE DIFFERENCE: { abs(diff) }')
                        flow_variable[element_name].pprint()
                        raise Exception()
                    print('')
                else:
                    print(f'{name}: {element_name}')
                    print(f'flow : {flow_array[i]}')
                    print(f'model: {pyo.value(flow_variable[element_name])}')
                    print('')


## to hold the indicies of the violations
## in the model or block
def add_monitored_flow_tracker(mb):
    mb._idx_monitored = list()
    mb._interfaces_monitored = list()
    mb._contingencies_monitored = list()

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

    if time is None: # DCOPF
        active_slack_tol = mb._ptdf_options['active_flow_tol']
    else: # Unit Commitment
        active_slack_tol = mb.parent_block()._ptdf_options['active_flow_tol']

    ## PFV -- power flow vector
    ## PFV_I -- interface power flow vector
    ## VA -- bus voltage angle vector
    PFV, PFV_I, VA = PTDF.calculate_masked_PFV(mb)

    violations_store = _MaximalViolationsStore(max_viol_add=max_viol_add, md=md, time=time, prepend_str=prepend_str)

    if len(PTDF.branches_keys_masked) > 0:
        violations_store.check_and_add_violations('branch', PFV, mb.pf,
                                            PTDF.lazy_branch_limits, PTDF.enforced_branch_limits,
                                           -PTDF.lazy_branch_limits, -PTDF.enforced_branch_limits,
                                            mb._idx_monitored, PTDF.branches_keys_masked)
    
    if len(PTDF.interface_keys) > 0:
        violations_store.check_and_add_violations('interface', PFV_I, mb.pfi,
                                            PTDF.lazy_interface_max_limits, PTDF.enforced_interface_max_limits,
                                            PTDF.lazy_interface_min_limits, PTDF.enforced_interface_min_limits,
                                            mb._interfaces_monitored, PTDF.interface_keys)

    if PTDF.contingencies and \
           violations_store.total_violations == 0:
        ## NOTE: checking contingency constraints in general could be very expensive
        ##       we probably want to delay doing so until we have a nearly transmission feasible
        ##       solution

        ## For each contingency, we'll only calculate the difference in flow,
        ## and check this against the difference in bounds, i.e.,

        ## power_flow_contingency == PFV + PFV_delta_c
        ## -rate_c <= power_flow_contingency <= +rate_c
        ## <===>
        ## -rate_c - PFV <= PFV_delta_c <= +rate_c - PFV
        ## <===>
        ## contingency_limits_lower <= PFV_delta_c <= contingency_limits_upper
        ## and
        ## contingency_limits_lower == -rate_c - PFV; contingency_limits_upper == rate_c - PFV

        ## In this way, we avoid (number of contingenies) adds PFV+PFV_delta_c

        logger.debug("Checking contingency flows...")
        lazy_contingency_limits_upper = PTDF.lazy_contingency_limits - PFV
        lazy_contingency_limits_lower = -PTDF.lazy_contingency_limits - PFV
        enforced_contingency_limits_upper = PTDF.enforced_contingency_limits - PFV
        enforced_contingency_limits_lower = -PTDF.enforced_contingency_limits - PFV
        for cn in PTDF.contingency_compensators:
            PFV_delta = PTDF.calculate_masked_PFV_delta(cn, PFV, VA)
            violations_store.check_and_add_violations('contingency', PFV_delta, mb.pfc,
                                                      lazy_contingency_limits_upper, enforced_contingency_limits_upper,
                                                      lazy_contingency_limits_lower, enforced_contingency_limits_lower,
                                                      mb._contingencies_monitored, PTDF.branches_keys_masked,
                                                      outer_name = cn, PFV = PFV)

    logger.debug(f"branches_monitored: {mb._idx_monitored}\n"
                 f"interfaces_monitored: {mb._interfaces_monitored}\n"
                 f"contingencies_monitored: {mb._contingencies_monitored}\n"
                 f"Violations being added: {violations_store.violations_store}\n"
                 f"Violations in model: {violations_store.monitored_violations}\n")

    viol_lazy = _LazyViolations(branch_lazy_violations=set(violations_store.get_violations_named('branch')),
                                interface_lazy_violations=set(violations_store.get_violations_named('interface')),
                                contingency_lazy_violations=set(violations_store.get_violations_named('contingency')))
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

    ## get the lines we're monitoring
    idx_monitored = mb._idx_monitored
    interfaces_monitored = mb._interfaces_monitored
    contingencies_monitored = mb._contingencies_monitored

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

    for name, constr in mb.ineq_pf_contingency_branch_thermal_bounds.items():
        slack = constr.slack()
        if slack_tol <= abs(slack):
            logger.debug(prepend_str+_generate_flow_monitor_remove_message('contingeny', name, abs(slack), baseMVA, time))
            constr_to_remove.append(constr)
            ## remove the index from the lines we're monitoring
            contingencies_monitored.remove((name[0], branchname_index_map[name[1]])) ## TODO: name?

    msg = prepend_str+"removing {} inactive transmission constraint(s)".format(len(constr_to_remove))
    if time is not None:
        msg += " at time {}".format(time)
    logger.debug(msg)

    for constr in constr_to_remove:
        if persistent_solver:
            solver.remove_constraint(constr)
        del constr
    return len(constr_to_remove)

def _generate_flow_viol_warning(expr, e_type, bn, flow, limit, baseMVA, time):
    ret_str = "WARNING: {0} {1} is in the  monitored set".format(e_type,bn)
    if time is not None:
        ret_str += " at time {}".format(time)
    ret_str += ", but flow exceeds limit!!\n\t flow={:.2f}, limit={:.2f}".format(flow*baseMVA, limit*baseMVA)
    ret_str += ", model_flow={:.2f}".format(pyo.value(expr[bn])*baseMVA)
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

## helper for generating pfc
def _iter_over_cont_viol_set(cont_viol_set, mb, PTDF, abs_ptdf_tol, rel_ptdf_tol):
    for (cn, i_b) in cont_viol_set:
        bn = PTDF.branches_keys_masked[i_b]
        if (cn, bn) not in mb._contingency_set:
            mb._contingency_set.add((cn,bn))
        if mb.pfc[cn, bn].expr is None:
            expr = libbranch.get_contingency_power_flow_expr_ptdf_approx(mb, cn, bn, PTDF,
                                                        abs_ptdf_tol=abs_ptdf_tol, rel_ptdf_tol=rel_ptdf_tol)
            mb.pfc[cn, bn] = expr
        yield cn, bn, i_b

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
        # initialize to 0.
        neg_slack.value = 0.
        pos_slack.value = 0.
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
        # initialize to 0.
        neg_slack.value = 0.
        pos_slack.value = 0.
    else:
        neg_slack = None
        pos_slack = None
        new_var = False

    return libbranch.generate_thermal_bounds(mb.pfi[i_n], minimum_limit, maximum_limit, neg_slack, pos_slack), new_var

def _generate_contingency_bounds(mb, cn, minimum_limit, maximum_limit):
    if cn in mb.pfc_slack_pos.index_set():
        if cn not in mb.pfc_slack_pos:
            neg_slack = mb.pfc_slack_neg[cn]
            pos_slack = mb.pfc_slack_pos[cn]
            assert len(mb.pfc_slack_pos) == len(mb.pfc_slack_neg)
            new_var = True
        else: # the constraint could have been added and removed
            neg_slack = mb.pfc_slack_neg[cn]
            pos_slack = mb.pfc_slack_pos[cn]
            new_var = False
        # initialize to 0.
        neg_slack.value = 0.
        pos_slack.value = 0.
    else:
        neg_slack = None
        pos_slack = None
        new_var = False

    return libbranch.generate_thermal_bounds(mb.pfc[cn], minimum_limit, maximum_limit, neg_slack, pos_slack), new_var

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

    _add_interface_violations(lazy_violations, flows, mb, md, solver, ptdf_options,
                                PTDF, model, baseMVA, persistent_solver, rel_ptdf_tol, abs_ptdf_tol,
                                time, prepend_str)
    _add_contingency_violations(lazy_violations, flows, mb, md, solver, ptdf_options,
                                PTDF, model, baseMVA, persistent_solver, rel_ptdf_tol, abs_ptdf_tol,
                                time, prepend_str)

def _add_interface_violations(lazy_violations, flows, mb, md, solver, ptdf_options,
                              PTDF, model, baseMVA, persistent_solver, rel_ptdf_tol, abs_ptdf_tol,
                              time, prepend_str):
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
        print(f"adding constraint for interface {i_n}")

def _add_contingency_violations(lazy_violations, flows, mb, md, solver, ptdf_options,
                                PTDF, model, baseMVA, persistent_solver, rel_ptdf_tol, abs_ptdf_tol,
                                time, prepend_str):
    ## in case there's no contingencies
    if not hasattr(mb, 'ineq_pf_contingency_branch_thermal_bounds'):
        return
    constr = mb.ineq_pf_contingency_branch_thermal_bounds
    contingencies_monitored = mb._contingencies_monitored
    for cn, bn, i_b in _iter_over_cont_viol_set(lazy_violations.contingency_lazy_violations, mb, PTDF, abs_ptdf_tol, rel_ptdf_tol):
        emergency_thermal_limit = PTDF.contingency_limits_array_masked[i_b]
        logger.debug(prepend_str+_generate_flow_monitor_message('contingency', (cn,bn), time=time))
        constr[cn,bn], new_slacks = _generate_contingency_bounds(mb, (cn,bn), -emergency_thermal_limit, emergency_thermal_limit)
        contingencies_monitored.append((cn, i_b))
        if new_slacks:
            m = model
            obj_coef = pyo.value(m.TimePeriodLengthHours*m.ContingencyLimitPenalty)

            if persistent_solver:
                if m is not m.model():
                    raise RuntimeError("Cannot add lazy var for branch slacks if part of a larger model")
                ## update the objective through the add_column method
                solver.add_column(m, mb.pfc_slack_pos[cn,bn], obj_coef, [], [])
                solver.add_column(m, mb.pfc_slack_neg[cn,bn], obj_coef, [], [])
            else:
                m.ContingencyViolationCost[time].expr += (obj_coef*mb.pfc_slack_pos[cn,bn] + \
                                                          obj_coef*mb.pfc_slack_neg[cn,bn] )
        if persistent_solver:
            solver.add_constraint(constr[cn,bn])

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

def copy_active_to_next_time(m, b_next, PTDF_next, slacks, slacks_I, slacks_C):
    active_slack_tol = m._ptdf_options['active_flow_tol']

    branchname_index_map = PTDF_next.branchname_to_index_masked_map
    interfacename_index_map = PTDF_next.interfacename_to_index_map

    viol_lazy = set()
    int_viol_lazy = set()
    cont_viol_lazy = set()

    idx_monitored = b_next._idx_monitored
    interfaces_monitored = b_next._interfaces_monitored
    contingencies_monitored = b_next._contingencies_monitored

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

    for cn, slack in slacks_C.items():
        if abs(slack) <= active_slack_tol:
            ## in case the topology has changed
            c, bn = cn
            if bn in branchname_index_map:
                bi = branchname_index_map[bn]
                if (c, bi) not in contingencies_monitored:
                    cont_viol_lazy.add((c, bi))

    flows = _CalculatedFlows()
    viol_lazy = _LazyViolations(branch_lazy_violations=viol_lazy,
                                interface_lazy_violations=int_viol_lazy,
                                contingency_lazy_violations=cont_viol_lazy)

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
        vars_to_load.extend(b.pfc_slack_pos.values())
        vars_to_load.extend(b.pfc_slack_neg.values())
    # XpressPersistent raises an exception if
    # this list is empty
    if vars_to_load:
        solver.load_vars(vars_to_load)
