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

def _get_uc_model(model_data, formulation_list, relax_binaries):
    formulation = UCFormulation(*formulation_list)
    md = scale_ModelData_to_pu(model_data)
    return generate_model(md, formulation, relax_binaries)

def create_tight_unit_commitment_model(model_data,
                                       network_constraints='power_balance_constraints',
                                       relaxed=False):
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
    return _get_uc_model(model_data, formulation_list, relaxed)

def create_compact_unit_commitment_model(model_data,
                                         network_constraints='power_balance_constraints',
                                         relaxed=False):
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
    return _get_uc_model(model_data, formulation_list, relaxed)

def create_KOW_unit_commitment_model(model_data,
                                     network_constraints='power_balance_constraints',
                                     relaxed=False):
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
    return _get_uc_model(model_data, formulation_list, relaxed)

def create_ALS_unit_commitment_model(model_data,
                                     network_constraints='power_balance_constraints',
                                     relaxed=False):
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
    return _get_uc_model(model_data, formulation_list, relaxed)

def create_MLR_unit_commitment_model(model_data,
                                     network_constraints='power_balance_constraints',
                                     relaxed=False):

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
    return _get_uc_model(model_data, formulation_list, relaxed)

def create_random1_unit_commitment_model(model_data,
                                         network_constraints='power_balance_constraints',
                                         relaxed=False):
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
    return _get_uc_model(model_data, formulation_list, relaxed)

def create_random2_unit_commitment_model(model_data,
                                         network_constraints='power_balance_constraints',
                                         relaxed=False):
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
    return _get_uc_model(model_data, formulation_list, relaxed)

def create_OAV_unit_commitment_model(model_data,
                                     network_constraints='power_balance_constraints',
                                     relaxed=False):
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
    return _get_uc_model(model_data, formulation_list, relaxed)

def create_OAV_tighter_unit_commitment_model(model_data,
                                             network_constraints='power_balance_constraints',
                                             relaxed=False):
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
    return _get_uc_model(model_data, formulation_list, relaxed)

def create_OAV_original_unit_commitment_model(model_data,
                                              network_constraints='power_balance_constraints',
                                              relaxed=False):
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
    return _get_uc_model(model_data, formulation_list, relaxed)

def create_OAV_up_downtime_unit_commitment_model(model_data,
                                                 network_constraints='power_balance_constraints',
                                                 relaxed=False):
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
    return _get_uc_model(model_data, formulation_list, relaxed)

def create_CA_unit_commitment_model(model_data,
                                    network_constraints='power_balance_constraints',
                                    relaxed=False):
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
    return _get_uc_model(model_data, formulation_list, relaxed)

def _set_options(solver, mipgap, timelimit, other_options):
    solver_name = solver.name

    if 'gurobi' in solver_name:
        solver.options.MIPGap = mipgap
        if timelimit is not None:
            solver.options.TimeLimit = timelimit
    elif 'cplex' in solver_name:
        solver.options.mip_tolerances_mipgap = mipgap
    elif 'glpk' in solver_name:
        solver.options.mipgap = mipgap
    elif 'cbc' in solver_name:
        solver.options.ratioGap = mipgap
    else:
        raise Exception('Solver {0} not recognized'.format(solver_name))

    for key, opt in other_options:
        solver.options[key] = opt

def _time_series_dict(values):
    return {'data_type':'time_series', 'values':values}

def solve_unit_commitment(model_data,
                          solver,
                          mipgap = 0.001,
                          timelimit = None,
                          solver_tee = True,
                          symbolic_solver_labels = False,
                          options = dict(),
                          uc_model_generator=create_tight_unit_commitment_model,
                          relaxed=False):
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

    '''
    import pyomo.environ as pe
    from pyomo.opt import SolverFactory, TerminationCondition
    from pyomo.solvers.plugins.solvers.persistent_solver import PersistentSolver

    ## termination conditions which are acceptable
    safe_termination_conditions = [ 
                                   TerminationCondition.maxTimeLimit,
                                   TerminationCondition.maxIterations,
                                   TerminationCondition.minFunctionValue,
                                   TerminationCondition.minStepLength,
                                   TerminationCondition.globallyOptimal,
                                   TerminationCondition.locallyOptimal,
                                   TerminationCondition.feasible,
                                   TerminationCondition.optimal,
                                   TerminationCondition.maxEvaluations,
                                   TerminationCondition.other,
                                  ]

    m = uc_model_generator(model_data, relaxed=relaxed)

    if relaxed:
        m.dual = pe.Suffix(direction=pe.Suffix.IMPORT)

    if isinstance(solver, str):
        solver = SolverFactory(solver)
    elif isinstance(solver, pyomo.opt.base.OptSolver):
        pass
    else:
        raise Exception('solver must be string or an instanciated pyomo solver')

    _set_options(solver, mipgap, timelimit, options) 

    if isinstance(solver, PersistentSolver):
        solver.set_instance(m, symbolic_solver_labels=symbolic_solver_labels)
        results = solver.solve(m, timelimit=timelimit, tee=solver_tee)
    else:
        results = solver.solve(m, timelimit=timelimit, tee=solver_tee, \
                              symbolic_solver_labels=symbolic_solver_labels)

    if results.solver.termination_condition not in safe_termination_conditions:
        raise Exception('Problem encountered during solve, termination_condition {}'.format(results.solver.terminataion_condition))

    md = m.model_data

    thermal_gens = dict(md.elements(element_type='generator', generator_type='thermal'))
    renewable_gens = dict(md.elements(element_type='generator', generator_type='renewable'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    storage = dict(md.elements(element_type='storage'))
    zones = dict(md.elements(element_type='zone'))
    areas = dict(md.elements(element_type='area'))

    data_time_periods = md.data['system']['time_indices']
    reserve_requirement = ('reserve_requirement' in md.data['system'])
    value = pe.value

    for g,g_dict in thermal_gens.items():
        pg_dict = {}
        if reserve_requirement:
            rg_dict = {}
        commitment_dict = {}
        commitment_cost_dict = {}
        production_cost_dict = {}
        for dt, mt in zip(data_time_periods,m.TimePeriods):
            if value(m.ThermalGeneratorForcedOutage[g,mt]):
                pg_dict[dt] = 0.
                if reserve_requirement:
                    rg_dict[dt] = 0.
                spin_g_dict[dt] = 0.
                commitment_dict[dt] = 0
                commitment_cost_dict[dt] = 0.
                production_cost_dict[dt] = 0.
            else:
                pg_dict[dt] = value(m.PowerGenerated[g,mt])
                if reserve_requirement:
                    rg_dict[dt] = value(m.ReserveProvided[g,mt])
                commitment_dict[dt] = value(m.UnitOn[g,mt])
                commitment_cost_dict[dt] = value(m.StartupCost[g,mt]+m.ShutdownCost[g,mt]+\
                                        m.MinimumProductionCost[g]*m.UnitOn[g,mt]*m.TimePeriodLengthHours)
                production_cost_dict[dt] = value(m.ProductionCost[g,mt])

                ## NOTE: we may want different hooks for these as we better integrate the
                ##       ancillary service model
                ## TODO: need add hooks for different kinds of ancillary services
                if m.ancillary_services:
                    pass


        g_dict['pg'] = _time_series_dict(pg_dict)
        if reserve_requirement:
            g_dict['rg'] = _time_series_dict(rg_dict)
        g_dict['commitment'] = _time_series_dict(commitment_dict)
        g_dict['commitment_cost'] = _time_series_dict(commitment_cost_dict)
        g_dict['production_cost'] = _time_series_dict(production_cost_dict)

    for g,g_dict in renewable_gens.items():
        pg_dict = {}
        for dt, mt in zip(data_time_periods,m.TimePeriods):
            if value(m.NondispatchableGeneratorForcedOutage[g,mt]):
                pg_dict[dt] = 0.
            else:
                pg_dict[dt] = value(m.NondispatchablePowerUsed[g,mt])
        g_dict['pg'] = _time_series_dict(pg_dict)

    for s,s_dict in storage.items():
        state_of_charge_dict = {}
        p_discharge_dict = {}
        p_charge_dict = {}
        operational_cost_dict = {}
        for dt, mt in zip(data_time_periods,m.TimePeriods):
            if value(m.StorageForceOutage[g,mt]):
                p_discharge_dict[dt] = 0.
                p_charge_dict[dt] = 0.
                operational_cost_dict[dt] = 0.
                # NOTE: if it goes offline, we'll assume it reverts to its minimum SOC 
                # TODO: Is this reasonable? 
                state_of_charge_dict[dt] = s_dict['minimum_state_of_charge']
            else:
                p_discharge_dict[dt] = value(m.PowerOutputStorage[s,mt])
                p_charge_dict[dt] = value(m.PowerInputStorage[s,mt])
                operational_cost_dict[dt] = value(m.StorageCost[s,mt])
                state_of_charge_dict[dt] = value(m.SocStorage[s.mt])

        s_dict['p_discharge'] = _time_series_dict(p_discharge_dict)
        s_dict['p_charge'] = _time_series_dict(p_charge_dict)
        s_dict['operational_cost'] = _time_series_dict(operational_cost_dict)
        s_dict['state_of_charge'] = _time_series_dict(state_of_charge_dict)

    ## NOTE: UC model currently has no notion of separate loads

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


    if m.ancillary_services:
        ## TODO: Can the code above this be re-factored in a similar way?
        ## as we add more zonal reserve products, they can be added here
        _zonal_reserve_map = { 'spinning_reserve_requirement' : { 'shortfall' : 'spinning_reserve_shortfall',
                                                                  'price'     : 'spinning_reserve_price',
                                                                  'shortfall_m' : m.ZonalSpinningReserveShortfall,
                                                                  'balance_m' : m.EnforceZonalSpinningReserveRequirement,
                                                                 },
                               'regulation_up_requirement' : { 'shorfall' : 'regulation_up_shortfall',
                                                               'price'    : 'regulation_up_price',
                                                               'shortfall_m' : m.ZonalRegulationUpShorfall,
                                                               'balance_m' : m.EnforceZonalRegulationUpRequirements,
                                                              },
                               'regulation_down_requirement' : { 'shortfall' : 'regulation_down_shortfall',
                                                                 'price'     : 'regulation_down_price',
                                                                 'shortfall_m' : m.ZonalRegulationDnShortfall,
                                                                 'balance_m' : m.EnforceZonalRegulationDnRequirements,
                                                                },
                               'flexible_ramp_up_requirement' : { 'shorfall' : 'flexible_ramp_up_shortfall',
                                                                  'price' : 'flexible_ramp_up_price',
                                                                  'shortfall_m' : m.ZonalFlexUpShortfall,
                                                                  'balance_m' : m.ZonalFlexUpRequirementConstr,
                                                                },
                               'flexible_ramp_down_requirement' : { 'shorfall' : 'flexible_ramp_down_shortfall',
                                                                    'price'    : 'flexible_ramp_down_price',
                                                                    'shortfall_m' : m.ZonalFlexDnShortfall,
                                                                    'balance_m' : m.ZonalFlexDnRequirementConstr,
                                                                   },
                               }

        ## as we add more system reserve products, they can be added here
        _system_reserve_map = { 'spinning_reserve_requirement' : { 'shortfall' : 'spinning_reserve_shortfall',
                                                                  'price'     : 'spinning_reserve_price',
                                                                  'shortfall_m' : m.SystemSpinningReserveShortfall,
                                                                  'balance_m' : m.EnforceSystemSpinningReserveRequirement,
                                                                 },
                               'regulation_up_requirement' : { 'shorfall' : 'regulation_up_shortfall',
                                                               'price'    : 'regulation_up_price',
                                                               'shortfall_m' : m.SystemRegulationUpShortfall,
                                                               'balance_m' : m.EnforceSystemRegulationUpRequirement,
                                                              },
                               'regulation_down_requirement' : { 'shortfall' : 'regulation_down_shortfall',
                                                                 'price'     : 'regulation_down_price',
                                                                 'shortfall_m' : m.SystemRegulationDnShortfall,
                                                                 'balance_m' : m.EnforceSystemRegulationDnRequirement,
                                                                },
                               'flexible_ramp_up_requirement' : { 'shorfall' : 'flexible_ramp_up_shortfall',
                                                                  'price' : 'flexible_ramp_up_price',
                                                                  'shortfall_m' : m.SystemFlexUpShortfall,
                                                                  'balance_m' : m.SystemFlexUpRequirementConstr,
                                                                },
                               'flexible_ramp_down_requirement' : { 'shorfall' : 'flexible_ramp_down_shortfall',
                                                                    'price'    : 'flexible_ramp_down_price',
                                                                    'shortfall_m' : m.SystemFlexDnShortfall,
                                                                    'balance_m' : m.SystemFlexDnRequirementConstr,
                                                                   },
                               }


        def _populate_zonal_reserves(elements_dict, string_handle):
            for e,e_dict in elements_dict.items():
                me = string_handle+e
                for req, req_dict in _zonal_reserve_map.items():
                    if req in e_dict:
                        req_shortfall_dict = {}
                        for dt, mt in zip(data_time_periods, m.TimePeriods):
                            req_shortfall_dict[dt] = value(req_dict['shorfall_m'][me,mt])
                        e_dict[req_dict['shorfall']] = _time_series_dict(req_shortfall_dict)
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
                        req_shortfall_dict[dt] = value(req_dict['shorfall_m'][mt])
                    sys_dict[req_dict['shortfall']] = _time_series_dict(req_shortfall_dict)
                    if relaxed:
                        req_price_dict = {}
                        for dt, mt in zip(data_time_periods, m.TimePeriods):
                            req_price_dict[dt] = value(m.dual[req_dict['balance_m'][mt]])
                        sys_dict[req_dict['price']] = _time_series_dict(req_price_dict)
        
        _populate_zonal_reserves(areas, 'area_')
        _populate_zonal_reserves(zones, 'zone_')

        _populate_system_reserves(md.data['system'])
    ## end if

    unscale_ModelData_to_pu(md, inplace=True)
    
    return md

