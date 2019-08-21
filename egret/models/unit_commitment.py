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
from egret.model_library.transmission.tx_utils import \
        scale_ModelData_to_pu, unscale_ModelData_to_pu
import egret.common.lazy_ptdf_utils as lpu
import egret.data.data_utils as data_utils
import numpy as np

def _get_uc_model(model_data, formulation_list, relax_binaries, **kwargs):
    formulation = UCFormulation(*formulation_list)
    md = model_data.clone_in_service()
    scale_ModelData_to_pu(md, inplace=True)
    return generate_model(md, formulation, relax_binaries, **kwargs)

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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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
        # TODO: describe the required and optional attributes
    network_constraints : str (optional)
        Set of network constraints to use. The default option uses a B-\\theta
        "DC" network.
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

def _lazy_ptdf_uc_solve_loop(m, md, solver, timelimit, solver_tee=True, symbolic_solver_labels=False, iteration_limit=100000, warn_on_max_iter=True, add_all_lazy_violations=False):
    from pyomo.solvers.plugins.solvers.persistent_solver import PersistentSolver

    persistent_solver = isinstance(solver, PersistentSolver)
    duals = hasattr(m, 'dual')

    ptdf_options = m._ptdf_options

    PVF = dict()
    viol_num = dict()
    viols_tup = dict()
    all_viol_in_model = dict()

    for i in range(iteration_limit):
        for t in m.TimePeriods:
            b = m.TransmissionBlock[t]

            PTDF = b._PTDF

            ## these should only need to be gathered once
            if PTDF.lazy_branch_limits is None: 
                branch_limits = PTDF.branch_limits_array
                rel_flow_tol = ptdf_options['rel_flow_tol']
                abs_flow_tol = ptdf_options['abs_flow_tol']
                lazy_rel_flow_tol = ptdf_options['lazy_rel_flow_tol']

                ## if lazy_rel_flow_tol < 0, this narrows the branch limits
                PTDF.lazy_branch_limits = branch_limits*(1+lazy_rel_flow_tol)
                ## only enforce the relative and absolute, within tollerance
                PTDF.enforced_branch_limits = np.maximum(branch_limits*(1+rel_flow_tol), branch_limits+abs_flow_tol)

            PVF[t], viol_num[t], viols_tup[t] = lpu.check_violations(b, PTDF)

        total_viol_num = sum(val for val in viol_num.values())
        print("iteration {0}, found {1} violation(s)".format(i,total_viol_num))
        if total_viol_num <= 0 and not add_all_lazy_violations:
            if persistent_solver and duals:
                solver.load_duals()
            return lpu.LazyPTDFTerminationCondition.NORMAL

        all_times_all_viol_in_model = True
        for t in m.TimePeriods:
            b = m.TransmissionBlock[t]

            PTDF = b._PTDF

            all_viol_in_model[t] = lpu.add_violations(viols_tup[t], PVF[t], b, md, solver, ptdf_options, PTDF, time=t)
            if all_times_all_viol_in_model and (not all_viol_in_model[t]):
                all_times_all_viol_in_model = False

        if all_times_all_viol_in_model or (total_viol_num <= 0.):
            if all_times_all_viol_in_model:
                print('WARNING: Terminating with monitored violations!')
                print('         Result is not transmission feasible.')
            if persistent_solver and duals:
                solver.load_duals()
            return lpu.LazyPTDFTerminationCondition.FLOW_VIOLATION

        if persistent_solver:
            solver.solve(m, tee=solver_tee, load_solutions=False, save_results=False)
            solver.load_vars()
        else:
            m.preprocess()
            solver.solve(m, tee=solver_tee, symbolic_solver_labels=symbolic_solver_labels)

    else:
        if warn_on_max_iter:
            print('WARNING: Exiting on maximum iterations for lazy PTDF model.')
            print('         Result is not transmission feasible.')
        if persistent_solver and duals:
            solver.load_duals()
        return lpu.LazyPTDFTerminationCondition.ITERATION_LIMIT


def _time_series_dict(values):
    return {'data_type':'time_series', 'values':values}

def solve_unit_commitment(model_data,
                          solver,
                          mipgap = 0.001,
                          timelimit = None,
                          solver_tee = True,
                          symbolic_solver_labels = False,
                          options = None,
                          uc_model_generator = create_tight_unit_commitment_model,
                          relaxed = False,
                          return_model = False,
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
    options : dict (optional)
        Other options to pass into the solver. Default is dict().
    uc_model_generator : function (optional)
        Function for generating the unit commitment model. Default is 
        egret.models.unit_commitment.create_tight_unit_commitment_model
    relaxed : bool (optional)
        If True, creates a relaxed unit commitment model
    return_model : bool (optional)
        If True, returns the pyomo model object
    kwargs : dictionary (optional)
        Additional arguments for building model
    '''

    import pyomo.environ as pe
    from pyomo.environ import value
    from egret.common.solver_interface import _solve_model

    m = uc_model_generator(model_data, relaxed=relaxed, **kwargs)

    network = ('branch' in model_data.data['elements']) and bool(len(model_data.data['elements']['branch']))

    if relaxed:
        m.dual = pe.Suffix(direction=pe.Suffix.IMPORT)
        m.rc = pe.Suffix(direction=pe.Suffix.IMPORT)
    elif m.power_balance == 'ptdf_power_flow' and m._ptdf_options['lazy'] and network:
        ## if we were asked to, serialize the PTDF matrices we calculated
        data_utils.write_ptdf_potentially_to_file(m._ptdf_options, m._PTDFs)
        ## BK -- be a bit more aggressive bring in PTDF constraints
        ##       from the relaxation
        relax_add_flow_tol = 0.05
        m._ptdf_options['lazy_rel_flow_tol'] -= relax_add_flow_tol

        iter_limit = m._ptdf_options['lp_iteration_limit']

        lpu.uc_instance_binary_relaxer(m, None)
        m, results, solver = _solve_model(m,solver,mipgap,timelimit,solver_tee,symbolic_solver_labels,options, return_solver=True)
        lp_termination_cond = _lazy_ptdf_uc_solve_loop(m, model_data, solver, timelimit, solver_tee=solver_tee,iteration_limit=iter_limit, warn_on_max_iter=False, add_all_lazy_violations=True)

        lpu.uc_instance_binary_enforcer(m, solver)
        m._ptdf_options['lazy_rel_flow_tol'] += relax_add_flow_tol

        iter_limit = m._ptdf_options['iteration_limit']
        solver.solve(m, tee=solver_tee)
        termination_cond = _lazy_ptdf_uc_solve_loop(m, model_data, solver, timelimit, solver_tee=solver_tee, iteration_limit=iter_limit, warn_on_max_iter=True)

    else:
        m, results, solver = _solve_model(m,solver,mipgap,timelimit,solver_tee,symbolic_solver_labels,options, return_solver=True)

    md = m.model_data

    # save results data to ModelData object
    thermal_gens = dict(md.elements(element_type='generator', generator_type='thermal'))
    renewable_gens = dict(md.elements(element_type='generator', generator_type='renewable'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    storage = dict(md.elements(element_type='storage'))
    zones = dict(md.elements(element_type='zone'))
    areas = dict(md.elements(element_type='area'))

    data_time_periods = md.data['system']['time_indices']
    reserve_requirement = ('reserve_requirement' in md.data['system'])

    regulation = False
    spin = False
    nspin = False
    supp = False
    flex = False
    if hasattr(m, 'regulation_service'):
        regulation = True
    if hasattr(m, 'spinning_reserve'):
        spin = True
    if hasattr(m, 'non_spinning_reserve'):
        nspin = True
    if hasattr(m, 'supplemental_reserve'):
        supp = True
    if hasattr(m, 'flexible_ramping'):
        flex = True

    fs = False
    if hasattr(m, 'fuel_supply'):
        fs = True
    fc = False
    if hasattr(m, 'fuel_consumption'):
        fc = True

    for g,g_dict in thermal_gens.items():
        pg_dict = {}
        if reserve_requirement:
            rg_dict = {}
        commitment_dict = {}
        commitment_cost_dict = {}
        production_cost_dict = {}
        ramp_up_avail_dict = {}

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
                                          'EnforceMaxAvailableRampDownRates',
                                          'EnforceMaxCapacity',
                                         ]
        ramp_up_avail_constrs = []
        for constr in ramp_up_avail_potential_constrs:
            if hasattr(m, constr):
                ramp_up_avail_constrs.append(getattr(m, constr))

        if regulation:
            reg_prov = {}
            reg_up_supp = {}
            reg_dn_supp = {}
        if spin:
            spin_supp = {}
        if nspin:
            nspin_supp = {}
        if supp:
            supp_supp = {}
        if flex:
            flex_up_supp = {}
            flex_dn_supp = {}
        gfs = (fs and (g in m.FuelSupplyGenerators))
        if gfs:
            fuel_consumed = {}
        gdf = (fc and (g in m.DualFuelGenerators))
        if gdf:
            aux_fuel_consumed = {}
        gdsf = (gdf and (g in m.SingleFireDualFuelGenerators))
        if gdsf:
            aux_fuel_indicator = {}


        for dt, mt in zip(data_time_periods,m.TimePeriods):
            pg_dict[dt] = value(m.PowerGenerated[g,mt])
            if reserve_requirement:
                rg_dict[dt] = value(m.ReserveProvided[g,mt])
            commitment_dict[dt] = value(m.UnitOn[g,mt])
            commitment_cost_dict[dt] = value(m.ShutdownCost[g,mt])
            if g in m.DualFuelGenerators:
                commitment_cost_dict[dt] += value(m.DualFuelCommitmentCost[g,mt])
                production_cost_dict[dt] = value(m.DualFuelProductionCost[g,mt])
            else:
                commitment_cost_dict[dt] += value(m.NoLoadCost[g,mt]+m.StartupCost[g,mt])
                production_cost_dict[dt] = value(m.ProductionCost[g,mt])

            if regulation:
                if g in m.AGC_Generators:
                    reg_prov[dt] = value(m.RegulationOn[g,mt])
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
        pg_dict = {}
        for dt, mt in zip(data_time_periods,m.TimePeriods):
            pg_dict[dt] = value(m.NondispatchablePowerUsed[g,mt])
        g_dict['pg'] = _time_series_dict(pg_dict)

    for s,s_dict in storage.items():
        state_of_charge_dict = {}
        p_discharge_dict = {}
        p_charge_dict = {}
        operational_cost_dict = {}
        for dt, mt in zip(data_time_periods,m.TimePeriods):
            p_discharge_dict[dt] = value(m.PowerOutputStorage[s,mt])
            p_charge_dict[dt] = value(m.PowerInputStorage[s,mt])
            operational_cost_dict[dt] = value(m.StorageCost[s,mt])
            state_of_charge_dict[dt] = value(m.SocStorage[s.mt])

        s_dict['p_discharge'] = _time_series_dict(p_discharge_dict)
        s_dict['p_charge'] = _time_series_dict(p_charge_dict)
        s_dict['operational_cost'] = _time_series_dict(operational_cost_dict)
        s_dict['state_of_charge'] = _time_series_dict(state_of_charge_dict)

    ## NOTE: UC model currently has no notion of separate loads

    if m.power_balance == 'btheta_power_flow':
        for l,l_dict in branches.items():
            pf_dict = {}
            for dt, mt in zip(data_time_periods,m.TimePeriods):
                pf_dict[dt] = value(m.TransmissionBlock[mt].pf[l])
            l_dict['pf'] = _time_series_dict(pf_dict)

        for b,b_dict in buses.items():
            va_dict = {}
            p_balance_violation_dict = {}
            pl_dict = {}
            for dt, mt in zip(data_time_periods,m.TimePeriods):
                va_dict[dt] = value(m.TransmissionBlock[mt].va[b])
                p_balance_violation_dict[dt] = value(m.LoadGenerateMismatch[b,mt])
                pl_dict[dt] = value(m.TransmissionBlock[mt].pl[b])
            b_dict['va'] = _time_series_dict(va_dict)
            b_dict['p_balance_violation'] = _time_series_dict(p_balance_violation_dict)
            b_dict['pl'] = _time_series_dict(pl_dict)
            if relaxed:
                lmp_dict = {}
                for dt, mt in zip(data_time_periods,m.TimePeriods):
                    lmp_dict[dt] = value(m.dual[m.TransmissionBlock[mt].eq_p_balance[b]])
                b_dict['lmp'] = _time_series_dict(lmp_dict)

    elif m.power_balance == 'ptdf_power_flow':
        flows_dict = dict()
        if relaxed:
            lmps_dict = dict()
        for mt in m.TimePeriods:
            b = m.TransmissionBlock[mt]
            flows_dict[mt] = dict()
            PTDF = b._PTDF
            PTDFM = PTDF.PTDFM
            branches_idx = PTDF.branches_keys
            NWV = np.array([pe.value(b.p_nw[bus]) for bus in PTDF.bus_iterator()])
            PFV = np.dot(PTDFM, NWV)
            for i,bn in enumerate(branches_idx):
                flows_dict[mt][bn] = PFV[i]
            if relaxed:
                PFD = np.zeros(len(branches_idx))
                for i,bn in enumerate(branches_idx):
                    if bn in b.ineq_pf_branch_thermal_lb:
                        PFD[i] += value(m.dual[b.ineq_pf_branch_thermal_lb[bn]])
                    if bn in b.ineq_pf_branch_thermal_ub:
                        PFD[i] += value(m.dual[b.ineq_pf_branch_thermal_ub[bn]])
                ## TODO: PFD is likely to be sparse, implying we just need a few
                ##       rows of the PTDF matrix (or columns in its transpose).
                LMPC = np.dot(-PTDFM.T, PFD)
                LMPE = value(m.dual[b.eq_p_balance])
                buses_idx = PTDF.buses_keys
                lmps_dict[mt] = dict()
                for i,bn in enumerate(buses_idx):
                    lmps_dict[mt][bn] = LMPE + LMPC[i]

        for l,l_dict in branches.items():
            pf_dict = {}
            for dt, mt in zip(data_time_periods,m.TimePeriods):
                ## if the key doesn't exist, it is because that line was out
                pf_dict[dt] = flows_dict[mt].get(l, 0.)
            l_dict['pf'] = _time_series_dict(pf_dict)

        for b,b_dict in buses.items():
            va_dict = {}
            p_balance_violation_dict = {}
            pl_dict = {}
            for dt, mt in zip(data_time_periods,m.TimePeriods):
                p_balance_violation_dict[dt] = value(m.LoadGenerateMismatch[b,mt])
                pl_dict[dt] = value(m.TransmissionBlock[mt].pl[b])
            b_dict['p_balance_violation'] = _time_series_dict(p_balance_violation_dict)
            b_dict['pl'] = _time_series_dict(pl_dict)
            if relaxed:
                lmp_dict = {}
                for dt, mt in zip(data_time_periods,m.TimePeriods):
                    lmp_dict[dt] = lmps_dict[mt][b]
                b_dict['lmp'] = _time_series_dict(lmp_dict)

    elif m.power_balance == 'power_balance_constraints':
        for l,l_dict in branches.items():
            pf_dict = {}
            for dt, mt in zip(data_time_periods,m.TimePeriods):
                pf_dict[dt] = value(m.LinePower[l,mt])
            l_dict['pf'] = _time_series_dict(pf_dict)

        for b,b_dict in buses.items():
            va_dict = {}
            p_balance_violation_dict = {}
            for dt, mt in zip(data_time_periods,m.TimePeriods):
                va_dict[dt] = value(m.Angle[b,mt])
                p_balance_violation_dict[dt] = value(m.LoadGenerateMismatch[b,mt])
            b_dict['va'] = _time_series_dict(va_dict)
            b_dict['p_balance_violation'] = _time_series_dict(p_balance_violation_dict)
            if relaxed:
                lmp_dict = {}
                for dt, mt in zip(data_time_periods,m.TimePeriods):
                    lmp_dict[dt] = value(m.dual[m.PowerBalance[b,mt]])
                b_dict['lmp'] = _time_series_dict(lmp_dict)
    elif m.power_balance in ['copperplate_power_flow', 'copperplate_relaxed_power_flow']:
        sys_dict = md.data['system']
        p_viol_dict = {}
        for dt, mt in zip(data_time_periods,m.TimePeriods):
            p_viol_dict[dt] = sum(value(m.LoadGenerateMismatch[b,mt]) for b in m.Buses)
        sys_dict['p_balance_violation'] = _time_series_dict(p_viol_dict)
        if relaxed:
            p_price_dict = {}
            for dt, mt in zip(data_time_periods,m.TimePeriods):
                p_price_dict[dt] = value(m.dual[m.TransmissionBlock[mt].eq_p_balance])
            sys_dict['p_price'] = _time_series_dict(p_price_dict)
    else:
        raise Exception("Unrecongized network type "+m.power_balance)


    if reserve_requirement:
        ## populate the system attributes
        sys_dict = md.data['system']
        sr_s_dict = {}
        for dt, mt in zip(data_time_periods,m.TimePeriods):
            sr_s_dict[dt] = value(m.ReserveShortfall[mt])
        sys_dict['reserve_shortfall'] = _time_series_dict(sr_s_dict)
        if relaxed:
            sr_p_dict = {}
            for dt, mt in zip(data_time_periods,m.TimePeriods):
                ## TODO: if the 'relaxed' flag is set, we should automatically
                ##       pick a formulation which uses the MLR reserve constraints
                sr_p_dict[dt] = value(m.dual[m.EnforceReserveRequirements[mt]])
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
                    req_shortfall_dict = {}
                    for dt, mt in zip(data_time_periods, m.TimePeriods):
                        req_shortfall_dict[dt] = value(req_dict['shortfall_m'][me,mt])
                    e_dict[req_dict['shortfall']] = _time_series_dict(req_shortfall_dict)
                    if relaxed:
                        req_price_dict = {}
                        for dt, mt in zip(data_time_periods, m.TimePeriods):
                            req_price_dict[dt] = value(m.dual[req_dict['balance_m'][me,mt]])
                        e_dict[req_dict['price']] = _time_series_dict(req_price_dict)

    def _populate_system_reserves(sys_dict):
        for req, req_dict in _system_reserve_map.items():
            if req in sys_dict:
                req_shortfall_dict = {}
                for dt, mt in zip(data_time_periods, m.TimePeriods):
                    req_shortfall_dict[dt] = value(req_dict['shortfall_m'][mt])
                sys_dict[req_dict['shortfall']] = _time_series_dict(req_shortfall_dict)
                if relaxed:
                    req_price_dict = {}
                    for dt, mt in zip(data_time_periods, m.TimePeriods):
                        req_price_dict[dt] = value(m.dual[req_dict['balance_m'][mt]])
                    sys_dict[req_dict['price']] = _time_series_dict(req_price_dict)
    
    _populate_zonal_reserves(areas, 'area_')
    _populate_zonal_reserves(zones, 'zone_')

    _populate_system_reserves(md.data['system'])

    if fs:
        fuel_supplies = dict(md.elements(element_type='fuel_supply'))
        for f, f_dict in fuel_supplies.items():
            fuel_consumed = {}
            fuel_supply_type = f_dict['fuel_supply_type']
            if fuel_supply_type == 'instantaneous':
                for dt, mt in zip(data_time_periods, m.TimePeriods):
                    fuel_consumed[dt] = value(m.TotalFuelConsumedAtInstFuelSupply[f,mt])
            else:
                print('WARNING: unrecongized fuel_supply_type {} for fuel_supply {}'.format(fuel_supply_type, f))
            f_dict['fuel_consumed'] = _time_series_dict(fuel_consumed)

    md.data['system']['total_cost'] = value(m.TotalCostObjective)

    unscale_ModelData_to_pu(md, inplace=True)
    
    if return_model:
        return md, m
    return md

# if __name__ == '__main__':
#     from egret.data.model_data import ModelData
#
#     file = "tests/uc_test_instances/test_case_1.json"
#     md = ModelData()
#     md.read_from_json(file)
#     solve_unit_commitment(md, "gurobi")
