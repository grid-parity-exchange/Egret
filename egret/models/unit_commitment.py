#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
This module provides functions that create the models for some typical
unit commitment formulations

# TODO: documentation
'''

from egret.model_library.unit_commitment.uc_model_generator \
        import UCFormulation, generate_model 
from egret.common.log import logger
from egret.model_library.transmission.tx_utils import unscale_ModelData_to_pu
from pyomo.solvers.plugins.solvers.persistent_solver import PersistentSolver
from pyomo.solvers.plugins.solvers.gurobi_persistent import GurobiPersistent 

import egret.common.lazy_ptdf_utils as lpu
import egret.data.ptdf_utils as ptdf_utils
import pyomo.environ as pe
import numpy as np

def _get_uc_model(model_data, formulation_list, relax_binaries, **kwargs):
    formulation = UCFormulation(*formulation_list)
    return generate_model(model_data, formulation, relax_binaries, **kwargs)

def create_tight_unit_commitment_model(model_data,
                                       network_constraints='ptdf_power_flow',
                                       relaxed=False,
                                       **kwargs):
    '''
    Create a new unit commitment model based on the "Tight" formulation from
    B. Knueven, J. Ostrowski, and J.-P. Watson. "On Mixed Integer Programming
    Formulations for the Unit Commitment Problem" (2018). Available:
    http://www.optimization-online.org/DB_FILE/2018/11/6930.pdf

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''

    formulation_list = [
                        'garver_3bin_vars',
                        'garver_power_vars',
                        'garver_power_avail_vars',
                        'pan_guan_gentile_KOW_generation_limits',
                        'damcikurt_ramping',
                        'KOW_production_costs_tightened',
                        'rajan_takriti_UT_DT',
                        'KOW_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_compact_unit_commitment_model(model_data,
                                         network_constraints='ptdf_power_flow',
                                         relaxed=False,
                                         **kwargs):
    '''
    Create a new unit commitment model based on the "Compact" formulation from
    B. Knueven, J. Ostrowski, and J.-P. Watson. "On Mixed Integer Programming
    Formulations for the Unit Commitment Problem" (2018). Available:
    http://www.optimization-online.org/DB_FILE/2018/11/6930.pdf

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''

    formulation_list = [
                         'garver_3bin_vars',
                         'garver_power_vars',
                         'garver_power_avail_vars',
                         'MLR_generation_limits',
                         'damcikurt_ramping',
                         'HB_production_costs',
                         'rajan_takriti_UT_DT',
                         'MLR_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_KOW_unit_commitment_model(model_data,
                                     network_constraints='ptdf_power_flow',
                                     relaxed=False,
                                     **kwargs):
    '''
    Create a new unit commitment model based on the formulation from
    B. Knueven, J. Ostrowski, and J.-P. Watson. "A Novel Matching 
    Formulation for Startup Costs in Unit Commitment" (2018).
    Available: http://www.optimization-online.org/DB_FILE/2017/03/5897.pdf

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''

    formulation_list = [
                         'garver_3bin_vars',
                         'garver_power_vars', 
                         'MLR_reserve_vars',
                         'MLR_generation_limits',
                         'MLR_ramping', 
                         'KOW_production_costs',
                         'rajan_takriti_UT_DT', 
                         'KOW_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_ALS_unit_commitment_model(model_data,
                                     network_constraints='ptdf_power_flow',
                                     relaxed=False,
                                     **kwargs):
    '''
    Create a new unit commitment model based on the formulation from
    Atakan, Semih, Guglielmo Lulli, and Suvrajeet Sen. "A state transition 
    MIP formulation for the unit commitment problem." IEEE Transactions on 
    Power Systems 33.1 (2018): 736-748.

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''

    formulation_list = [
                         'ALS_state_transition_vars',
                         'garver_power_vars',
                         'CA_power_avail_vars',
                         'OAV_generation_limits',
                         'ALS_damcikurt_ramping',
                         'basic_production_costs_envelope',
                         'rajan_takriti_UT_DT_2bin',
                         'ALS_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_MLR_unit_commitment_model(model_data,
                                     network_constraints='ptdf_power_flow',
                                     relaxed=False,
                                     **kwargs):
    '''
    Create a new unit commitment model based on the formulation from
    Morales-España, Germán, Jesus M. Latorre, and Andres Ramos. "Tight and
    compact MILP formulation for the thermal unit commitment problem." IEEE
    Transactions on Power Systems 28.4 (2013): 4897-4908.

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''
    formulation_list = [
                         'garver_3bin_vars',
                         'garver_power_vars',
                         'MLR_reserve_vars',
                         'MLR_generation_limits',
                         'MLR_ramping',
                         'HB_production_costs',
                         'rajan_takriti_UT_DT',
                         'MLR_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_random1_unit_commitment_model(model_data,
                                         network_constraints='ptdf_power_flow',
                                         relaxed=False,
                                         **kwargs):
    '''
    Create a new unit commitment model based on the "Random1" formulation from
    B. Knueven, J. Ostrowski, and J.-P. Watson. "On Mixed Integer Programming
    Formulations for the Unit Commitment Problem" (2018). Available:
    http://www.optimization-online.org/DB_FILE/2018/11/6930.pdf

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''

    formulation_list = [
                        'ALS_state_transition_vars',
                        'garver_power_vars',
                        'garver_power_avail_vars',
                        'MLR_generation_limits',
                        'ALS_damcikurt_ramping',
                        'HB_production_costs',
                        'rajan_takriti_UT_DT_2bin',
                        'MLR_startup_costs2',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_random2_unit_commitment_model(model_data,
                                         network_constraints='ptdf_power_flow',
                                         relaxed=False,
                                         **kwargs):
    '''
    Create a new unit commitment model based on the "Random2" formulation from
    B. Knueven, J. Ostrowski, and J.-P. Watson. "On Mixed Integer Programming
    Formulations for the Unit Commitment Problem" (2018). Available:
    http://www.optimization-online.org/DB_FILE/2018/11/6930.pdf

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''

    formulation_list = [
                        'garver_3bin_vars',
                        'garver_power_vars',
                        'garver_power_avail_vars',
                        'pan_guan_gentile_KOW_generation_limits',
                        'OAV_ramping_enhanced_2period',
                        'HB_production_costs',
                        'rajan_takriti_UT_DT',
                        'MLR_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_OAV_unit_commitment_model(model_data,
                                     network_constraints='ptdf_power_flow',
                                     relaxed=False,
                                     **kwargs):
    '''
    Create a new unit commitment model based on the formulation from
    Ostrowski, James, Miguel F. Anjos, and Anthony Vannelli. "Tight mixed 
    integer linear programming formulations for the unit commitment problem." 
    IEEE Transactions on Power Systems 27.1 (2012): 39-46.

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''

    formulation_list = [
                        'garver_3bin_vars',
                        'basic_power_vars',
                        'CA_power_avail_vars',
                        'OAV_generation_limits',
                        'OAV_ramping_enhanced',
                        'HB_production_costs',
                        'rajan_takriti_UT_DT',
                        'CA_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_OAV_tighter_unit_commitment_model(model_data,
                                             network_constraints='ptdf_power_flow',
                                            relaxed=False,
                                            **kwargs):
    '''
    Create a new unit commitment model based on the formulation from
    Ostrowski, James, Miguel F. Anjos, and Anthony Vannelli. "Tight mixed 
    integer linear programming formulations for the unit commitment problem." 
    IEEE Transactions on Power Systems 27.1 (2012): 39-46.

    with the two-period ramping can multi-period variable-upper-bound
    constraints therein.

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''

    formulation_list = [
                         'garver_3bin_vars',
                         'basic_power_vars',
                         'CA_power_avail_vars',
                         'OAV_generation_limits_enhanced',
                         'OAV_ramping_enhanced_2period',
                         'HB_production_costs',
                         'rajan_takriti_UT_DT',
                         'CA_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_OAV_original_unit_commitment_model(model_data,
                                              network_constraints='ptdf_power_flow',
                                              relaxed=False,
                                              **kwargs):
    '''
    Create a new unit commitment model based on the "original" formulation from
    Ostrowski, James, Miguel F. Anjos, and Anthony Vannelli. "Tight mixed 
    integer linear programming formulations for the unit commitment problem." 
    IEEE Transactions on Power Systems 27.1 (2012): 39-46.

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''

    formulation_list = [
                         'garver_3bin_vars',
                         'basic_power_vars',
                         'CA_power_avail_vars',
                         'OAV_generation_limits',
                         'arroyo_conejo_ramping',
                         'HB_production_costs',
                         'DEKT_UT_DT',
                         'CA_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_OAV_up_downtime_unit_commitment_model(model_data,
                                                 network_constraints='ptdf_power_flow',
                                                 relaxed=False,
                                                 **kwargs):
    '''
    Create a new unit commitment model based on the "up/downtime" formulation from
    Ostrowski, James, Miguel F. Anjos, and Anthony Vannelli. "Tight mixed 
    integer linear programming formulations for the unit commitment problem." 
    IEEE Transactions on Power Systems 27.1 (2012): 39-46.

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''

    formulation_list = [
                         'garver_3bin_vars',
                         'basic_power_vars',
                         'CA_power_avail_vars',
                         'OAV_generation_limits',
                         'arroyo_conejo_ramping',
                         'HB_production_costs',
                         'rajan_takriti_UT_DT',
                         'CA_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_CA_unit_commitment_model(model_data,
                                    network_constraints='ptdf_power_flow',
                                    relaxed=False,
                                    **kwargs):
    '''
    Create a new unit commitment model based on the formulation from
    Carrión, Miguel, and José M. Arroyo. "A computationally efficient
    mixed-integer linear formulation for the thermal unit commitment
    problem." IEEE Transactions on power systems 21.3 (2006): 1371-1378.

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''

    formulation_list = [
                         'CA_1bin_vars',
                         'basic_power_vars',
                         'CA_power_avail_vars',
                         'CA_generation_limits',
                         'CA_ramping_limits',
                         'CA_production_costs',
                         'CA_UT_DT',
                         'CA_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed, **kwargs)

def create_CHP_unit_commitment_model(model_data,
                                     network_constraints='ptdf_power_flow',
                                     relaxed=False,
                                     **kwargs):
    '''
    Create a new unit commitment model based on the "extensive form" convex hull
    pricing formulation from Knueven, Ostrowski, Castillo, and Watson (2019) 
    "A computationally efficient algorithm for computing convex hull prices".
    pre-print available: http://www.optimization-online.org/DB_FILE/2019/09/7370.pdf

    NOTE: This model asserts that certain products (storage, dual-fuel units, 
          ancillary services) are not part of the problem, so as to accurately
          return convex hull prices.

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''
    from egret.model_library.unit_commitment.thermal_convex_hull import \
            add_convex_hull_for_all_units

    formulation_list = [
                         'garver_3bin_vars',
                         'garver_power_vars',
                         'MLR_reserve_vars',
                         'gentile_generation_limits',
                         'arroyo_conejo_ramping',
                         'KOW_production_costs_tightened',
                         'rajan_takriti_UT_DT_2bin',
                         'pochet_wolsey_startup_costs',
                         network_constraints,
                       ]

    model = _get_uc_model(model_data, formulation_list, relaxed)

    add_convex_hull_for_all_units(model)

    return model

def create_super_tight_unit_commitment_model(model_data,
                                             network_constraints='ptdf_power_flow',
                                             relaxed=False,
                                             **kwargs):
    '''
    Create a new unit commitment formulation based using the tightest formulation available
    for each component, subject to avoiding quadratic (or worse) blow up in problem size

    Used as the "master problem" pricing formulation in Knueven, Ostrowski, Castillo, 
    and Watson (2019) "A computationally efficient algorithm for computing convex hull prices".
    pre-print available: http://www.optimization-online.org/DB_FILE/2019/09/7370.pdf

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    network_constraints : str (optional)
        Set of network constraints to use. 
        The default option uses lazy ptdfs, in this case no network 
        constraints are added to the model initially.
    relaxed : bool (optional)
        If True, creates a model with the binary variables relaxed to [0,1].
        Default is False.
    kwargs : dictionary (optional):
        Additional arguments for egret.model_library.unit_commitment.uc_model_generator.generate_model

    Returns
    -------
        pyomo.environ.ConcreteModel unit commitment model

    '''
    formulation_list = [
                         'garver_3bin_vars',
                         'garver_power_vars',
                         'MLR_reserve_vars',
                         'pan_guan_gentile_KOW_generation_limits',
                         'damcikurt_ramping',
                         'KOW_production_costs_super_tight',
                         'rajan_takriti_UT_DT_2bin',
                         'KOW_startup_costs',
                         network_constraints,
                       ]
    return _get_uc_model(model_data, formulation_list, relaxed)


def _lazy_ptdf_warmstart_copy_violations(m, md, t_subset, solver, ptdf_options, prepend_str):
    if len(t_subset) > 1:
        raise Exception("_lazy_ptdf_warmstart_copy_violations only handles t_subset with one element -- somebody should fix that!")
    slacks = None
    flows = dict()
    viol_lazy = dict()
    total_flow_constr_added = 0
    for t_o in m.TimePeriods:
        if t_o in t_subset:
            continue

        if slacks is None:
            slacks = dict()
            slacks_I = dict()
            slacks_C = dict()
            for t in t_subset:
                b_ = m.TransmissionBlock[t]

                ## only load the slacks once
                for bn, constr in b_.ineq_pf_branch_thermal_bounds.items():
                    slacks[bn] = constr.slack()
                ## if the interface slack variable is active, we will have 0
                ## slack in the constraint
                for i_n, constr in b_.ineq_pf_interface_bounds.items():
                    slacks_I[i_n] = constr.slack()
                for cn, constr in b_.ineq_pf_contingency_branch_thermal_bounds.items():
                    slacks_C[cn] = constr.slack()

        b_other = m.TransmissionBlock[t_o]
        PTDF_other = b_other._PTDF

        flows[t_o], viol_lazy[t_o] = lpu.copy_active_to_next_time(m,  b_other, PTDF_other, slacks, slacks_I, slacks_C)

        logger.debug(prepend_str+"adding {0} flow constraints at time {1}".format(len(viol_lazy[t_o]),t_o))
        lpu.add_violations(viol_lazy[t_o], flows[t_o], b_other, md, solver, ptdf_options, PTDF_other, time=t_o, prepend_str=prepend_str)
        total_flow_constr_added += len(viol_lazy[t_o])
    logger.info(prepend_str+"added {0} flow constraint(s)".format(total_flow_constr_added))

def _lazy_ptdf_solve(m, solver, persistent_solver, symbolic_solver_labels, solver_tee, vars_to_load, solve_method_options=None):
    if solve_method_options is None:
        solve_method_options = dict()
    if persistent_solver:
        results = solver.solve(m, tee=solver_tee, load_solutions=False, save_results=False, **solve_method_options)
        solver.load_vars(vars_to_load)
    else:
        results = solver.solve(m, tee=solver_tee, symbolic_solver_labels=symbolic_solver_labels, load_solutions=False, **solve_method_options)
        m.solutions.load_from(results)

def _lazy_ptdf_normal_terminatation(all_viol_in_model, results, i, prepend_str):
    if all_viol_in_model:
        logger.warning(prepend_str+'WARNING: Terminating with monitored violations! Result is not transmission feasible.')
        return lpu.LazyPTDFTerminationCondition.FLOW_VIOLATION, results, i
    return lpu.LazyPTDFTerminationCondition.NORMAL, results, i

def _lazy_ptdf_uc_solve_loop(m, md, solver, timelimit, solver_tee=True, symbolic_solver_labels=False, iteration_limit=100000, vars_to_load=None, 
        solve_method_options=None, add_all_lazy_violations=False, warmstart_loop=False, t_subset=None, vars_to_load_t_subset=None, prepend_str=""):
    '''
    The lazy PTDF unit commitment solver loop. This function iteratively
    adds violated transmission constraints until either the result is
    transmission feasible or we're tracking every violated constraint
    in the model

    Parameters
    ----------
    m : pyomo.environ.ConcreteModel
        An egret unit commitment model
    md : egret.data.ModelData
        An egret ModelData object
    solver : pyomo.opt.solver
        A pyomo solver object
    solver_tee : bool (optional)
        For displaying the solver log (default is True)
    symbolic_solver_labels : bool (optional)
        Use symbolic solver labels when writing to the solver (default is False)
    iteration_limit : int (optional)
        Number of iterations before a hard termination (default is 100000)
    vars_to_load : None, list (optional)
        Applies only to persistent solvers. If None, every primal variable is loaded.
        If a list, then should be a list of pyomo Vars to be loaded into the pyomo model
        at termination. Default is None.
    solve_method_options : None, dict (optional)
        Additional options passed into the solve method. Default is None
    add_all_lazy_violations : bool (optional)
        If True, on the termination iteration, additional violations or near violations
        will be added to the model before returning (and before re-solving). Used
        for the last lazy iteration on an LP relaxation of unit commitment
    warmstart_loop : bool (optional)
        If True, excutes a warmstart loop based on the time periods in t_subset
    t_subset : None, list (optional)
        Subset of TimePeriods to use for the warmstart loop
    vars_to_load_t_subset : None, list (optional)
        Subset of vars_to_load which are just from t_subset
    prepend_str : str (optional)
        String to prepend to all messages from the solve loop

    Returns
    -------
    egret.common.lazy_ptdf_utils.LazyPTDFTerminationCondition : the termination status
    pyomo.opt.results.SolverResults : The results object from the pyomo solver
    int : The number of iterations before termination

    '''

    persistent_solver = isinstance(solver, PersistentSolver)
    duals = hasattr(m, 'dual')

    results = None 

    ptdf_options = m._ptdf_options

    flows = dict()
    viol_num = dict()
    mon_viol_num = dict()
    viol_lazy = dict()

    if warmstart_loop:
        if t_subset is None:
            t_subset = [max(m.TotalDemand, key=m.TotalDemand.__getitem__)]
        time_periods = t_subset
        if vars_to_load_t_subset is None:
            vars_to_load_t_subset = vars_to_load
        vars_to_load_time_periods = vars_to_load_t_subset
        ## it doesn't make sense for both of these to be true
        assert add_all_lazy_violations is False
    else:
        time_periods = m.TimePeriods
        vars_to_load_time_periods = vars_to_load


    for i in range(iteration_limit):
        for t in time_periods:
            b = m.TransmissionBlock[t]

            PTDF = b._PTDF

            flows[t], viol_num[t], mon_viol_num[t], viol_lazy[t] = \
                    lpu.check_violations(b, md, PTDF, ptdf_options['max_violations_per_iteration'], time=t, prepend_str=prepend_str)

        total_viol_num = sum(viol_num.values())
        total_mon_viol_num = sum(mon_viol_num.values())

        ## this flag is for if we found violations **and** every violation is in the model
        all_viol_in_model = (total_viol_num > 0) and (total_viol_num == total_mon_viol_num)

        ## this flag is for if we're going to terminate this iteration,
        ## either because there are no violations in this solution
        ## **or** because every violation is already in the model
        terminate_this_iter = (total_viol_num == 0) or all_viol_in_model

        iter_status_str = prepend_str+"iteration {0}, found {1} violation(s)".format(i,total_viol_num)
        if total_mon_viol_num:
            iter_status_str += ", {} of which are already monitored".format(total_mon_viol_num)

        logger.info(iter_status_str)

        if terminate_this_iter and not add_all_lazy_violations:
            if warmstart_loop:
                if persistent_solver:
                    lpu._load_pf_slacks(solver, m, t_subset)
                _lazy_ptdf_warmstart_copy_violations(m, md, time_periods, solver, ptdf_options, prepend_str)
                results = _lazy_ptdf_solve(m, solver, persistent_solver, symbolic_solver_labels, solver_tee, vars_to_load, solve_method_options)
            if persistent_solver and duals and (results is not None) and (vars_to_load is None):
                solver.load_duals()
            return _lazy_ptdf_normal_terminatation(all_viol_in_model, results, i, prepend_str)

        total_flow_constr_added = 0
        for t in time_periods:
            b = m.TransmissionBlock[t]

            PTDF = b._PTDF

            lpu.add_violations(viol_lazy[t], flows[t], b, md, solver, ptdf_options, PTDF, time=t, prepend_str=prepend_str)
            total_flow_constr_added += len(viol_lazy[t])

        logger.info(prepend_str+"iteration {0}, added {1} flow constraint(s)".format(i,total_flow_constr_added))

        ## NOTE: Here we should not load additional variables as
        ##       we've added more constraints to the model
        if terminate_this_iter and add_all_lazy_violations:
            return _lazy_ptdf_normal_terminatation(all_viol_in_model, results, i, prepend_str)

        results = _lazy_ptdf_solve(m, solver, persistent_solver, symbolic_solver_labels, solver_tee, vars_to_load_time_periods, solve_method_options)

    else:
        logger.warning(prepend_str+'WARNING: Exiting on maximum iterations for lazy PTDF model. Result is not transmission feasible.')
        if warmstart_loop:
            if persistent_solver:
                lpu._load_pf_slacks(solver, m, t_subset)
            _lazy_ptdf_warmstart_copy_violations(m, md, time_periods, solver, ptdf_options)
            results = _lazy_ptdf_solve(m, solver, persistent_solver, symbolic_solver_labels, solver_tee, vars_to_load, solve_method_options)
        if persistent_solver and duals and (results is not None) and (vars_to_load is None):
            solver.load_duals()
        return lpu.LazyPTDFTerminationCondition.ITERATION_LIMIT, results, i

def _outer_lazy_ptdf_solve_loop(m, solver, mipgap, timelimit, solver_tee, symbolic_solver_labels, solver_options, solve_method_options, relaxed, set_instance=True):

    from egret.common.solver_interface import _solve_model
    import time

    egret_metasolver_status = dict()

    start_time = time.time()

    ## cache here the variables that need to be 
    ## loaded to check transimission feasbility
    ## for a persistent solver
    max_demand_time = max(m.TotalDemand, key=m.TotalDemand.__getitem__)
    t_subset = [max_demand_time, ]
    if isinstance(solver, PersistentSolver) or (isinstance(solver,str) and 'persistent' in solver):
        vars_to_load = list()
        vars_to_load_t_subset = list()
        for t in m.TimePeriods:
            b = m.TransmissionBlock[t]
            if isinstance(b.p_nw, pe.Var):
                vars_to_load.extend(b.p_nw.values())
                if t in t_subset:
                    vars_to_load_t_subset.extend(b.p_nw.values())
            else:
                vars_to_load = None
                vars_to_load_t_subset = None
                break
    else:
        vars_to_load = None
        vars_to_load_t_subset = None

    lp_iter_limit = m._ptdf_options['lp_iteration_limit']
    lp_warmstart_iter_limit = m._ptdf_options['pre_lp_iteration_limit']
    model_data = m.model_data

    ## if this is a MIP, iterate though a few times with just the LP relaxation
    if not relaxed and lp_iter_limit > 0:

        lpu.uc_instance_binary_relaxer(m, None)
        m, results_init, solver = _solve_model(m,solver,mipgap,timelimit,solver_tee,symbolic_solver_labels,solver_options,solve_method_options, return_solver=True, vars_to_load = vars_to_load, set_instance=set_instance)
        if lp_warmstart_iter_limit > 0:
            lp_warmstart_termination_cond, results, lp_warmstart_iterations = \
                    _lazy_ptdf_uc_solve_loop(m, model_data, solver, timelimit, solver_tee=solver_tee,iteration_limit=lp_warmstart_iter_limit, vars_to_load_t_subset = vars_to_load_t_subset, vars_to_load=vars_to_load, t_subset=t_subset, warmstart_loop=True, prepend_str="[LP warmstart phase] ")
            egret_metasolver_status['lp_warmstart_termination_cond'] = lp_warmstart_termination_cond
            egret_metasolver_status['lp_warmstart_iterations'] = lp_warmstart_iterations

        lp_termination_cond, results, lp_iterations = \
                _lazy_ptdf_uc_solve_loop(m, model_data, solver, timelimit, solver_tee=solver_tee,iteration_limit=lp_iter_limit, vars_to_load = vars_to_load, add_all_lazy_violations=True, prepend_str="[LP phase] ")
        ## if the initial solve was transmission feasible, then
        ## we never re-solved the problem
        if results is None:
            results = results_init

        egret_metasolver_status['lp_termination_cond'] = lp_termination_cond
        egret_metasolver_status['lp_iterations'] = lp_iterations

        if m._ptdf_options['lp_cleanup_phase']:
            tot_removed = 0
            for t,b in m.TransmissionBlock.items():
                tot_removed += lpu.remove_inactive(b, solver, t, prepend_str="[LP cleanup phase] ")
            logger.info("[LP cleanup phase] removed {0} inactive flow constraint(s)".format(tot_removed))

        lpu.uc_instance_binary_enforcer(m, solver)

        ## solve the MIP after enforcing binaries
        results_init = solver.solve(m, tee=solver_tee, load_solutions=False)
        if isinstance(solver, PersistentSolver):
            solver.load_vars(vars_to_load)
        else:
            m.solutions.load_from(results_init)

    ## else if relaxed or lp_iter_limit == 0, do an initial solve
    else:
        m, results_init, solver = _solve_model(m,solver,mipgap,timelimit,solver_tee,symbolic_solver_labels,solver_options, solve_method_options, return_solver=True, vars_to_load=vars_to_load, set_instance=set_instance)

    iter_limit = m._ptdf_options['iteration_limit']
    
    if relaxed and lp_warmstart_iter_limit > 0:
        lp_termination_cond, results, lp_iterations = \
                _lazy_ptdf_uc_solve_loop(m, model_data, solver, timelimit, solver_tee=solver_tee,iteration_limit=lp_warmstart_iter_limit, vars_to_load_t_subset = vars_to_load_t_subset, vars_to_load=vars_to_load, t_subset=t_subset, warmstart_loop=True, prepend_str="[LP warmstart phase] ")
    termination_cond, results, iterations = _lazy_ptdf_uc_solve_loop(m, model_data, solver, timelimit, solver_tee=solver_tee, iteration_limit=iter_limit, vars_to_load=vars_to_load, prepend_str=("[LP phase] " if relaxed else "[MIP phase] "))
    ## if the initial solve was transmission feasible, then
    ## we never re-solved the problem
    if results is None:
        results = results_init

    egret_metasolver_status['termination_cond'] = termination_cond
    egret_metasolver_status['iterations'] = iterations

    if isinstance(solver, PersistentSolver) and vars_to_load is not None:
        solver.load_vars()
        if hasattr(m, "dual"):
            solver.load_duals()
        if hasattr(m, "slack"):
            solver.load_slacks()

    egret_metasolver_status['time'] = time.time() - start_time
    results.egret_metasolver = egret_metasolver_status

    return m, results, solver

def _time_series_dict(values):
    return {'data_type':'time_series', 'values':values}

def _preallocated_list(other_iter):
    return [ None for _ in other_iter ]

def _save_uc_results(m, relaxed):
    from pyomo.environ import value

    # dual suffix is on top-level model
    if relaxed:
        dual = m.model().dual

    md = m.model_data

    # save results data to ModelData object
    thermal_gens = dict(md.elements(element_type='generator', generator_type='thermal'))
    renewable_gens = dict(md.elements(element_type='generator', generator_type='renewable'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    interfaces = dict(md.elements(element_type='interface'))
    contingencies = dict(md.elements(element_type='contingency'))
    storage = dict(md.elements(element_type='storage'))
    zones = dict(md.elements(element_type='zone'))
    areas = dict(md.elements(element_type='area'))
    pg_security_constraints = dict(md.elements(element_type='security_constraint', security_constraint_type='pg'))
    dc_branches = dict(md.elements(element_type='dc_branch'))

    data_time_periods = md.data['system']['time_keys']
    reserve_requirement = ('reserve_requirement' in md.data['system'])

    regulation = bool(m.regulation_service)
    spin = bool(m.spinning_reserve)
    nspin = bool(m.non_spinning_reserve)
    supp = bool(m.supplemental_reserve)
    flex = bool(m.flexible_ramping)

    fs = bool(m.fuel_supply)
    fc = bool(m.fuel_consumption)

    ## All prices are in $/(MW*time_period) by construction
    ## time_period_length_hours == hours/time_period, so
    ##    $/(MW*time_period)/time_period_length_hours
    ## == $/(MW*time_period) * (time_period/hours)
    ## == $/(MW*hours) == $/MWh.
    ## All dual values are divided by this quantity
    ## so as to report out $/MWh.
    time_period_length_hours = value(m.TimePeriodLengthHours)

    ## all of the potential constraints that could limit maximum output
    ## Not all unit commitment models have these constraints, so first
    ## we need check if they're on the model object
    ramp_up_avail_potential_constrs = [
                                      'EnforceMaxAvailableRampUpRates',
                                      'AncillaryServiceRampUpLimit',
                                      'power_limit_from_start',
                                      'power_limit_from_stop',
                                      'power_limit_from_start_stop',
                                      'power_limit_from_start_stops',
                                      'max_power_limit_from_starts',
                                      'EnforceMaxAvailableRampDownRates',
                                      'EnforceMaxCapacity',
                                      'OAVUpperBound',
                                      'EnforceGenerationLimits',
                                     ]
    ramp_up_avail_constrs = []
    for constr in ramp_up_avail_potential_constrs:
        if hasattr(m, constr):
            ramp_up_avail_constrs.append(getattr(m, constr))

    for g,g_dict in thermal_gens.items():
        pg_dict = _preallocated_list(data_time_periods)
        if reserve_requirement:
            rg_dict = _preallocated_list(data_time_periods)
        commitment_dict = _preallocated_list(data_time_periods)
        commitment_cost_dict = _preallocated_list(data_time_periods)
        production_cost_dict = _preallocated_list(data_time_periods)
        ramp_up_avail_dict = _preallocated_list(data_time_periods)


        if regulation:
            reg_prov = _preallocated_list(data_time_periods)
            reg_up_supp = _preallocated_list(data_time_periods)
            reg_dn_supp = _preallocated_list(data_time_periods)
        if spin:
            spin_supp = _preallocated_list(data_time_periods)
        if nspin:
            nspin_supp = _preallocated_list(data_time_periods)
        if supp:
            supp_supp = _preallocated_list(data_time_periods)
        if flex:
            flex_up_supp = _preallocated_list(data_time_periods)
            flex_dn_supp = _preallocated_list(data_time_periods)
        gfs = (fs and (g in m.FuelSupplyGenerators))
        if gfs:
            fuel_consumed = _preallocated_list(data_time_periods)
        gdf = (fc and (g in m.DualFuelGenerators))
        if gdf:
            aux_fuel_consumed = _preallocated_list(data_time_periods)
        gdsf = (gdf and (g in m.SingleFireDualFuelGenerators))
        if gdsf:
            aux_fuel_indicator = _preallocated_list(data_time_periods)


        for dt, mt in enumerate(m.TimePeriods):
            pg_dict[dt] = value(m.PowerGenerated[g,mt])
            if reserve_requirement:
                rg_dict[dt] = value(m.ReserveProvided[g,mt])
            if relaxed:
                commitment_dict[dt] = value(m.UnitOn[g,mt])
            else:
                commitment_dict[dt] = int(round(value(m.UnitOn[g,mt])))
            commitment_cost_dict[dt] = value(m.ShutdownCost[g,mt])
            if g in m.DualFuelGenerators:
                commitment_cost_dict[dt] += value(m.DualFuelCommitmentCost[g,mt])
                production_cost_dict[dt] = value(m.DualFuelProductionCost[g,mt])
            else:
                commitment_cost_dict[dt] += value(m.NoLoadCost[g,mt]+m.StartupCost[g,mt])
                production_cost_dict[dt] = value(m.ProductionCost[g,mt])

            if regulation:
                if g in m.AGC_Generators:
                    if relaxed:
                        reg_prov[dt] = value(m.RegulationOn[g,mt])
                    else:
                        reg_prov[dt] = int(round(value(m.RegulationOn[g,mt])))
                    reg_up_supp[dt] = value(m.RegulationReserveUp[g,mt])
                    reg_dn_supp[dt] = value(m.RegulationReserveDn[g,mt])
                    commitment_cost_dict[dt] += value(m.RegulationCostCommitment[g,mt])
                    production_cost_dict[dt] += value(m.RegulationCostGeneration[g,mt])
                else:
                    reg_prov[dt] = 0.
                    reg_up_supp[dt] = 0.
                    reg_dn_supp[dt] = 0.

            if spin:
                spin_supp[dt] = value(m.SpinningReserveDispatched[g,mt])
                production_cost_dict[dt] += value(m.SpinningReserveCostGeneration[g,mt])

            if nspin:
                if g in m.NonSpinGenerators:
                    nspin_supp[dt] = value(m.NonSpinningReserveDispatched[g,mt])
                    production_cost_dict[dt] += value(m.NonSpinningReserveCostGeneration[g,mt])
                else:
                    nspin_supp[dt] = 0.
            if supp:
                supp_supp[dt] = value(m.SupplementalReserveDispatched[g,mt])
                production_cost_dict[dt] += value(m.SupplementalReserveCostGeneration[g,mt])
            if flex:
                flex_up_supp[dt] = value(m.FlexUpProvided[g,mt])
                flex_dn_supp[dt] = value(m.FlexDnProvided[g,mt])
            if gfs:
                fuel_consumed[dt] = value(m.PrimaryFuelConsumed[g,mt])
            if gdsf:
                aux_fuel_indicator[dt] = value(m.UnitOnAuxFuel[g,mt])
            if gdf:
                aux_fuel_consumed[dt] = value(m.AuxiliaryFuelConsumed[g,mt])

            ## pyomo doesn't add constraints that are skiped to the index set, so we also
            ## need check here if the index exists.
            slack_list = []
            for constr in ramp_up_avail_constrs:
                if (g,mt) in constr:
                    slack_list.append(constr[g,mt].slack())

            ramp_up_avail_dict[dt] = min( slack_list )


        g_dict['pg'] = _time_series_dict(pg_dict)
        if reserve_requirement:
            g_dict['rg'] = _time_series_dict(rg_dict)
        g_dict['commitment'] = _time_series_dict(commitment_dict)
        g_dict['commitment_cost'] = _time_series_dict(commitment_cost_dict)
        g_dict['production_cost'] = _time_series_dict(production_cost_dict)
        if regulation:
            g_dict['reg_provider'] = _time_series_dict(reg_prov)
            g_dict['reg_up_supplied'] = _time_series_dict(reg_up_supp)
            g_dict['reg_down_supplied'] = _time_series_dict(reg_dn_supp)
        if spin:
            g_dict['spinning_supplied'] = _time_series_dict(spin_supp)
        if nspin:
            g_dict['non_spinning_supplied'] = _time_series_dict(nspin_supp)
        if supp:
            g_dict['supplemental_supplied'] = _time_series_dict(supp_supp)
        if flex:
            g_dict['flex_up_supplied'] = _time_series_dict(flex_up_supp)
            g_dict['flex_down_supplied'] = _time_series_dict(flex_dn_supp)
        if gfs:
            g_dict['fuel_consumed'] = _time_series_dict(fuel_consumed)
        if gdsf:
            g_dict['aux_fuel_status'] = _time_series_dict(aux_fuel_indicator)
        if gdf:
            g_dict['aux_fuel_consumed'] = _time_series_dict(aux_fuel_consumed)
        g_dict['headroom'] = _time_series_dict(ramp_up_avail_dict)

    for g,g_dict in renewable_gens.items():
        pg_dict = _preallocated_list(data_time_periods)
        for dt, mt in enumerate(m.TimePeriods):
            pg_dict[dt] = value(m.NondispatchablePowerUsed[g,mt])
        g_dict['pg'] = _time_series_dict(pg_dict)

    for s,s_dict in storage.items():
        state_of_charge_dict = _preallocated_list(data_time_periods)
        p_discharge_dict = _preallocated_list(data_time_periods)
        p_charge_dict = _preallocated_list(data_time_periods)
        operational_cost_dict = _preallocated_list(data_time_periods)
        for dt, mt in enumerate(m.TimePeriods):
            p_discharge_dict[dt] = value(m.PowerOutputStorage[s,mt])
            p_charge_dict[dt] = value(m.PowerInputStorage[s,mt])
            operational_cost_dict[dt] = value(m.StorageCost[s,mt])
            state_of_charge_dict[dt] = value(m.SocStorage[s,mt])

        s_dict['p_discharge'] = _time_series_dict(p_discharge_dict)
        s_dict['p_charge'] = _time_series_dict(p_charge_dict)
        s_dict['operational_cost'] = _time_series_dict(operational_cost_dict)
        s_dict['state_of_charge'] = _time_series_dict(state_of_charge_dict)

    for sc, sc_dict in pg_security_constraints.items():
        sc_violation = None
        sc_flow = _preallocated_list(data_time_periods)
        for dt, mt in enumerate(m.TimePeriods):
            b = m.TransmissionBlock[mt]
            sc_flow[dt] = value(b.pgSecurityExpression[sc])
            if sc in b.pgRelaxedSecuritySet:
                if sc_violation is None:
                    sc_violation = _preallocated_list(data_time_periods)
                sc_violation[dt] = value(b.pgSecuritySlackPos[sc] - b.pgSecuritySlackNeg[sc])
        sc_dict['pf'] = _time_series_dict(sc_flow)
        if sc_violation is not None:
            sc_dict['pf_violation'] = _time_series_dict(sc_violation)

    ## NOTE: UC model currently has no notion of separate loads

    if m.power_balance == 'btheta_power_flow':
        for l,l_dict in branches.items():
            pf_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                pf_dict[dt] = value(m.TransmissionBlock[mt].pf[l])
            l_dict['pf'] = _time_series_dict(pf_dict)
            if l in m.BranchesWithSlack:
                pf_violation_dict = _preallocated_list(data_time_periods)
                for dt, (mt, b) in enumerate(m.TransmissionBlock.items()):
                    if l in b.pf_slack_pos:
                        pf_violation_dict[dt] = value(b.pf_slack_pos[l] - b.pf_slack_neg[l])
                    else:
                        pf_violation_dict[dt] = 0.
                l_dict['pf_violation'] = _time_series_dict(pf_violation_dict)

        for b,b_dict in buses.items():
            va_dict = _preallocated_list(data_time_periods)
            p_balance_violation_dict = _preallocated_list(data_time_periods)
            pl_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                va_dict[dt] = value(m.TransmissionBlock[mt].va[b])
                p_balance_violation_dict[dt] = value(m.LoadGenerateMismatch[b,mt])
                pl_dict[dt] = value(m.TransmissionBlock[mt].pl[b])
            b_dict['va'] = _time_series_dict(va_dict)
            b_dict['p_balance_violation'] = _time_series_dict(p_balance_violation_dict)
            b_dict['pl'] = _time_series_dict(pl_dict)
            if relaxed:
                lmp_dict = _preallocated_list(data_time_periods)
                for dt, mt in enumerate(m.TimePeriods):
                    lmp_dict[dt] = value(dual[m.TransmissionBlock[mt].eq_p_balance[b]])/time_period_length_hours
                b_dict['lmp'] = _time_series_dict(lmp_dict)

        for i,i_dict in interfaces.items():
            pf_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                pf_dict[dt] = value(m.TransmissionBlock[mt].pfi[i])
            i_dict['pf'] = _time_series_dict(pf_dict)
            if i in m.InterfacesWithSlack:
                pf_violation_dict = _preallocated_list(data_time_periods)
                for dt, (mt, b) in enumerate(m.TransmissionBlock.items()):
                    if i in b.pfi_slack_pos:
                        pf_violation_dict[dt] = value(b.pfi_slack_pos[i] - b.pfi_slack_neg[i])
                    else:
                        pf_violation_dict[dt] = 0.
                i_dict['pf_violation'] = _time_series_dict(pf_violation_dict)

        for k,k_dict in dc_branches.items():
            pf_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                pf_dict[dt] = value(m.HVDCLinePower[k,mt])
            k_dict['pf'] = _time_series_dict(pf_dict)

    elif m.power_balance == 'ptdf_power_flow':
        flows_dict = dict()
        interface_flows_dict = dict()
        voltage_angle_dict = dict()
        if relaxed:
            lmps_dict = dict()
        for mt in m.TimePeriods:
            b = m.TransmissionBlock[mt]
            PTDF = b._PTDF

            branches_idx = PTDF.branches_keys
            PFV, PFV_I, VA = PTDF.calculate_PFV(b)

            flows_dict[mt] = dict()
            for i,bn in enumerate(branches_idx):
                flows_dict[mt][bn] = PFV[i]

            interface_idx = PTDF.interface_keys
            interface_flows_dict[mt] = dict()
            for i, i_n in enumerate(interface_idx):
                interface_flows_dict[mt][i_n] = PFV_I[i]

            buses_idx = PTDF.buses_keys
            voltage_angle_dict[mt] = dict()
            for i,bn in enumerate(buses_idx):
                voltage_angle_dict[mt][bn] = VA[i]

            if relaxed:
                LMP = PTDF.calculate_LMP(b, dual, b.eq_p_balance)
                lmps_dict[mt] = dict()
                for i,bn in enumerate(buses_idx):
                    lmps_dict[mt][bn] = LMP[i]

        for i,i_dict in interfaces.items():
            pf_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                pf_dict[dt] = interface_flows_dict[mt][i]
            i_dict['pf'] = _time_series_dict(pf_dict)
            if i in m.InterfacesWithSlack:
                pf_violation_dict = _preallocated_list(data_time_periods)
                for dt, (mt, b) in enumerate(m.TransmissionBlock.items()):
                    if i in b.pfi_slack_pos:
                        pf_violation_dict[dt] = value(b.pfi_slack_pos[i] - b.pfi_slack_neg[i])
                    else:
                        pf_violation_dict[dt] = 0.
                i_dict['pf_violation'] = _time_series_dict(pf_violation_dict)

        if contingencies:
            for dt, (mt, b) in enumerate(m.TransmissionBlock.items()):
                contingency_flows = PTDF.calculate_monitored_contingency_flows(b)
                for (cn,bn), flow in contingency_flows.items():
                    c_dict = contingencies[cn]
                    if 'monitored_branches' not in c_dict:
                        c_dict['monitored_branches'] = _time_series_dict([{} for _ in data_time_periods])
                    monitored_branches = c_dict['monitored_branches']['values'][dt]
                    monitored_branches[bn] = {'pf' : flow}
                    if (cn,bn) in b.pfc_slack_pos:
                        monitored_branches[bn]['pf_violation'] = value(b.pfc_slack_pos[cn,bn]-b.pfc_slack_neg[cn,bn])

        for l,l_dict in branches.items():
            pf_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                ## if the key doesn't exist, it is because that line was out
                pf_dict[dt] = flows_dict[mt].get(l, 0.)
            l_dict['pf'] = _time_series_dict(pf_dict)
            if l in m.BranchesWithSlack:
                pf_violation_dict = _preallocated_list(data_time_periods)
                for dt, (mt, b) in enumerate(m.TransmissionBlock.items()):
                    if l in b.pf_slack_pos:
                        pf_violation_dict[dt] = value(b.pf_slack_pos[l] - b.pf_slack_neg[l])
                    else:
                        pf_violation_dict[dt] = 0.
                l_dict['pf_violation'] = _time_series_dict(pf_violation_dict)

        for b,b_dict in buses.items():
            va_dict = _preallocated_list(data_time_periods)
            p_balance_violation_dict = _preallocated_list(data_time_periods)
            pl_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                p_balance_violation_dict[dt] = value(m.LoadGenerateMismatch[b,mt])
                pl_dict[dt] = value(m.TransmissionBlock[mt].pl[b])
                va_dict[dt] = voltage_angle_dict[mt][b]
            b_dict['p_balance_violation'] = _time_series_dict(p_balance_violation_dict)
            b_dict['pl'] = _time_series_dict(pl_dict)
            b_dict['va'] = _time_series_dict(va_dict)
            if relaxed:
                lmp_dict = _preallocated_list(data_time_periods)
                for dt, mt in enumerate(m.TimePeriods):
                    lmp_dict[dt] = lmps_dict[mt][b]/time_period_length_hours
                b_dict['lmp'] = _time_series_dict(lmp_dict)

        for k,k_dict in dc_branches.items():
            pf_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                pf_dict[dt] = value(m.HVDCLinePower[k,mt])
            k_dict['pf'] = _time_series_dict(pf_dict)

    elif m.power_balance == 'power_balance_constraints':
        for l,l_dict in branches.items():
            pf_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                pf_dict[dt] = value(m.LinePower[l,mt])
            l_dict['pf'] = _time_series_dict(pf_dict)
            if l in m.BranchesWithSlack:
                pf_violation_dict = _preallocated_list(data_time_periods)
                for dt, mt in enumerate(m.TimePeriods):
                    pf_violation_dict[dt] = value(m.BranchSlackPos[l,mt] - m.BranchSlackNeg[l,mt])
                l_dict['pf_violation'] = _time_series_dict(pf_violation_dict)

        for b,b_dict in buses.items():
            va_dict = _preallocated_list(data_time_periods)
            p_balance_violation_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                va_dict[dt] = value(m.Angle[b,mt])
                p_balance_violation_dict[dt] = value(m.LoadGenerateMismatch[b,mt])
            b_dict['va'] = _time_series_dict(va_dict)
            b_dict['p_balance_violation'] = _time_series_dict(p_balance_violation_dict)
            if relaxed:
                lmp_dict = _preallocated_list(data_time_periods)
                for dt, mt in enumerate(m.TimePeriods):
                    lmp_dict[dt] = value(dual[m.PowerBalance[b,mt]])/time_period_length_hours
                b_dict['lmp'] = _time_series_dict(lmp_dict)

        for i,i_dict in interfaces.items():
            pf_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                pf_dict[dt] = value(m.InterfaceFlow[i,mt])
            i_dict['pf'] = _time_series_dict(pf_dict)
            if i in m.InterfacesWithSlack:
                pf_violation_dict = _preallocated_list(data_time_periods)
                for dt, mt in enumerate(m.TimePeriods):
                    pf_violation_dict[dt] = value(m.InterfaceSlackPos[i,mt] - m.InterfaceSlackNeg[i,mt])
                i_dict['pf_violation'] = _time_series_dict(pf_violation_dict)

        for k,k_dict in dc_branches.items():
            pf_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                pf_dict[dt] = value(m.HVDCLinePower[k,mt])
            k_dict['pf'] = _time_series_dict(pf_dict)

    elif m.power_balance in ['copperplate_power_flow', 'copperplate_relaxed_power_flow']:
        sys_dict = md.data['system']
        p_viol_dict = _preallocated_list(data_time_periods)
        for dt, mt in enumerate(m.TimePeriods):
            p_viol_dict[dt] = sum(value(m.LoadGenerateMismatch[b,mt]) for b in m.Buses)
        sys_dict['p_balance_violation'] = _time_series_dict(p_viol_dict)
        if relaxed:
            p_price_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                p_price_dict[dt] = value(dual[m.TransmissionBlock[mt].eq_p_balance])/time_period_length_hours
            sys_dict['p_price'] = _time_series_dict(p_price_dict)
    else:
        raise Exception("Unrecongized network type "+m.power_balance)


    if reserve_requirement:
        ## populate the system attributes
        sys_dict = md.data['system']
        sr_s_dict = _preallocated_list(data_time_periods)
        for dt, mt in enumerate(m.TimePeriods):
            sr_s_dict[dt] = value(m.ReserveShortfall[mt])
        sys_dict['reserve_shortfall'] = _time_series_dict(sr_s_dict)
        if relaxed:
            sr_p_dict = _preallocated_list(data_time_periods)
            for dt, mt in enumerate(m.TimePeriods):
                ## TODO: if the 'relaxed' flag is set, we should automatically
                ##       pick a formulation which uses the MLR reserve constraints
                sr_p_dict[dt] = value(dual[m.EnforceReserveRequirements[mt]])/time_period_length_hours
            sys_dict['reserve_price'] = _time_series_dict(sr_p_dict)


    ## TODO: Can the code above this be re-factored in a similar way?
    ## as we add more zonal reserve products, they can be added here
    _zonal_reserve_map = dict()
    _system_reserve_map = dict()
    if spin:
        _zonal_reserve_map['spinning_reserve_requirement'] = {'shortfall' : 'spinning_reserve_shortfall',
                                                              'price'     : 'spinning_reserve_price',
                                                              'shortfall_m' : m.ZonalSpinningReserveShortfall,
                                                              'balance_m' : m.EnforceZonalSpinningReserveRequirement,
                                                             }
        _system_reserve_map['spinning_reserve_requirement'] = {'shortfall' : 'spinning_reserve_shortfall',
                                                               'price'     : 'spinning_reserve_price',
                                                               'shortfall_m' : m.SystemSpinningReserveShortfall,
                                                               'balance_m' : m.EnforceSystemSpinningReserveRequirement,
                                                              }
    if nspin:
        _zonal_reserve_map['non_spinning_reserve_requirement'] = {'shortfall' : 'non_spinning_reserve_shortfall',
                                                                  'price'     : 'non_spinning_reserve_price',
                                                                  'shortfall_m' : m.ZonalNonSpinningReserveShortfall,
                                                                  'balance_m' : m.EnforceNonSpinningZonalReserveRequirement,
                                                                 }
        _system_reserve_map['non_spinning_reserve_requirement'] = {'shortfall' : 'non_spinning_reserve_shortfall',
                                                                   'price'     : 'non_spinning_reserve_price',
                                                                   'shortfall_m' : m.SystemNonSpinningReserveShortfall,
                                                                   'balance_m' : m.EnforceSystemNonSpinningReserveRequirement,
                                                                   }
    if regulation:
        _zonal_reserve_map['regulation_up_requirement'] = {'shortfall' : 'regulation_up_shortfall',
                                                           'price'    : 'regulation_up_price',
                                                           'shortfall_m' : m.ZonalRegulationUpShortfall,
                                                           'balance_m' : m.EnforceZonalRegulationUpRequirements,
                                                          }
        _system_reserve_map['regulation_up_requirement'] = {'shortfall' : 'regulation_up_shortfall',
                                                            'price'    : 'regulation_up_price',
                                                            'shortfall_m' : m.SystemRegulationUpShortfall,
                                                            'balance_m' : m.EnforceSystemRegulationUpRequirement,
                                                           }
        _zonal_reserve_map['regulation_down_requirement'] = {'shortfall' : 'regulation_down_shortfall',
                                                             'price'     : 'regulation_down_price',
                                                             'shortfall_m' : m.ZonalRegulationDnShortfall,
                                                             'balance_m' : m.EnforceZonalRegulationDnRequirements,
                                                            }
        _system_reserve_map['regulation_down_requirement'] = {'shortfall' : 'regulation_down_shortfall',
                                                              'price'     : 'regulation_down_price',
                                                              'shortfall_m' : m.SystemRegulationDnShortfall,
                                                              'balance_m' : m.EnforceSystemRegulationDnRequirement,
                                                             } 
    if flex:
        _zonal_reserve_map['flexible_ramp_up_requirement'] = { 'shortfall' : 'flexible_ramp_up_shortfall',
                                                              'price' : 'flexible_ramp_up_price',
                                                              'shortfall_m' : m.ZonalFlexUpShortfall,
                                                              'balance_m' : m.ZonalFlexUpRequirementConstr,
                                                             }
        _system_reserve_map['flexible_ramp_up_requirement'] = {'shortfall' : 'flexible_ramp_up_shortfall',
                                                               'price' : 'flexible_ramp_up_price',
                                                               'shortfall_m' : m.SystemFlexUpShortfall,
                                                               'balance_m' : m.SystemFlexUpRequirementConstr,
                                                              }
        _zonal_reserve_map['flexible_ramp_down_requirement'] = {'shortfall' : 'flexible_ramp_down_shortfall',
                                                                'price'    : 'flexible_ramp_down_price',
                                                                'shortfall_m' : m.ZonalFlexDnShortfall,
                                                                'balance_m' : m.ZonalFlexDnRequirementConstr,
                                                               }
        _system_reserve_map['flexible_ramp_down_requirement'] = {'shortfall' : 'flexible_ramp_down_shortfall',
                                                                 'price'    : 'flexible_ramp_down_price',
                                                                 'shortfall_m' : m.SystemFlexDnShortfall,
                                                                 'balance_m' : m.SystemFlexDnRequirementConstr,
                                                                }
    if supp:
        _zonal_reserve_map['supplemental_reserve_requirement'] = {'shortfall' : 'supplemental_shortfall',
                                                                 'price' : 'supplemental_price',
                                                                 'shortfall_m' : m.ZonalSupplementalReserveShortfall,
                                                                 'balance_m' : m.EnforceZonalSupplementalReserveRequirement,
                                                                 }

        _system_reserve_map['supplemental_reserve_requirement'] = {'shortfall' : 'supplemental_shortfall',
                                                                   'price' : 'supplemental_price',
                                                                   'shortfall_m' : m.SystemSupplementalReserveShortfall,
                                                                   'balance_m' : m.EnforceSystemSupplementalReserveRequirement,
                                                                   }

    def _populate_zonal_reserves(elements_dict, string_handle):
        for e,e_dict in elements_dict.items():
            me = string_handle+e
            for req, req_dict in _zonal_reserve_map.items():
                if req in e_dict:
                    req_shortfall_dict = _preallocated_list(data_time_periods)
                    for dt, mt in enumerate(m.TimePeriods):
                        req_shortfall_dict[dt] = value(req_dict['shortfall_m'][me,mt])
                    e_dict[req_dict['shortfall']] = _time_series_dict(req_shortfall_dict)
                    if relaxed:
                        req_price_dict = _preallocated_list(data_time_periods)
                        for dt, mt in enumerate(m.TimePeriods):
                            req_price_dict[dt] = value(dual[req_dict['balance_m'][me,mt]])/time_period_length_hours
                        e_dict[req_dict['price']] = _time_series_dict(req_price_dict)

    def _populate_system_reserves(sys_dict):
        for req, req_dict in _system_reserve_map.items():
            if req in sys_dict:
                req_shortfall_dict = _preallocated_list(data_time_periods)
                for dt, mt in enumerate(m.TimePeriods):
                    req_shortfall_dict[dt] = value(req_dict['shortfall_m'][mt])
                sys_dict[req_dict['shortfall']] = _time_series_dict(req_shortfall_dict)
                if relaxed:
                    req_price_dict = _preallocated_list(data_time_periods)
                    for dt, mt in enumerate(m.TimePeriods):
                        req_price_dict[dt] = value(dual[req_dict['balance_m'][mt]])/time_period_length_hours
                    sys_dict[req_dict['price']] = _time_series_dict(req_price_dict)
    
    _populate_zonal_reserves(areas, 'area_')
    _populate_zonal_reserves(zones, 'zone_')

    _populate_system_reserves(md.data['system'])

    if fs:
        fuel_supplies = dict(md.elements(element_type='fuel_supply'))
        for f, f_dict in fuel_supplies.items():
            fuel_consumed = _preallocated_list(data_time_periods)
            fuel_supply_type = f_dict['fuel_supply_type']
            if fuel_supply_type == 'instantaneous':
                for dt, mt in enumerate(m.TimePeriods):
                    fuel_consumed[dt] = value(m.TotalFuelConsumedAtInstFuelSupply[f,mt])
            else:
                logger.warning('WARNING: unrecongized fuel_supply_type {} for fuel_supply {}'.format(fuel_supply_type, f))
            f_dict['fuel_consumed'] = _time_series_dict(fuel_consumed)

    md.data['system']['total_cost'] = value(m.TotalCostObjective)

    unscale_ModelData_to_pu(md, inplace=True)

    return md

def _solve_unit_commitment(m, solver, mipgap, timelimit, solver_tee, symbolic_solver_labels, solver_options, solve_method_options,relaxed, set_instance=True):
    from egret.common.solver_interface import _solve_model
    model_data = m.model_data
    network = ('branch' in model_data.data['elements']) and bool(len(model_data.data['elements']['branch']))
    if m.power_balance == 'ptdf_power_flow' and m._ptdf_options['lazy'] and network:
        return _outer_lazy_ptdf_solve_loop(m, solver, mipgap, timelimit, solver_tee, symbolic_solver_labels, solver_options, solve_method_options,relaxed, set_instance=set_instance )
    else:
        return _solve_model(m,solver,mipgap,timelimit,solver_tee,symbolic_solver_labels,solver_options,solve_method_options, return_solver=True, set_instance=set_instance)

def solve_unit_commitment(model_data,
                          solver,
                          mipgap = 0.001,
                          timelimit = None,
                          solver_tee = True,
                          symbolic_solver_labels = False,
                          solver_options = None,
                          solve_method_options = None,
                          uc_model_generator = create_tight_unit_commitment_model,
                          relaxed = False,
                          return_model = False,
                          return_results = False,
                          **kwargs):
    '''
    Create and solve a new unit commitment model

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
        # TODO: describe the required and optional attributes
    solver : str or pyomo.opt.base.solvers.OptSolver
        Either a string specifying a pyomo solver name, or an instanciated pyomo solver
    mipgap : float (optional)
        Mipgap to use for unit commitment solve; default is 0.001
    timelimit : float (optional)
        Time limit for unit commitment run. Default of None results in no time
        limit being set -- runs until mipgap is satisfied
    solver_tee : bool (optional)
        Display solver log. Default is True.
    symbolic_solver_labels : bool (optional)
        Use symbolic solver labels. Useful for debugging; default is False.
    solver_options : dict (optional)
        Other options to pass into the solver. Default is dict().
    solve_method_options : dict (optional)
        Other options to pass into the pyomo solve method. Default is dict().
    uc_model_generator : function (optional)
        Function for generating the unit commitment model. Default is 
        egret.models.unit_commitment.create_tight_unit_commitment_model
    relaxed : bool (optional)
        If True, creates a relaxed unit commitment model
    return_model : bool (optional)
        If True, returns the pyomo model object
    return_results : bool (optional)
        If True, returns the pyomo results object
    kwargs : dictionary (optional)
        Additional arguments for building model
    '''

    m = uc_model_generator(model_data, relaxed=relaxed, **kwargs)

    if relaxed:
        m.dual = pe.Suffix(direction=pe.Suffix.IMPORT)

    m, results, solver = _solve_unit_commitment(m, solver, mipgap, timelimit, solver_tee, symbolic_solver_labels, solver_options, solve_method_options,relaxed )

    md = _save_uc_results(m, relaxed)
    
    if return_model and return_results:
        return md, m, results
    elif return_model:
        return md, m
    elif return_results:
        return md, results
    return md

if __name__ == '__main__':
    from egret.data.model_data import ModelData

    filen = "tests/uc_test_instances/tiny_uc_tc_2.json"
    md = ModelData.read(filen)
    md_results = solve_unit_commitment(md, "cbc")
