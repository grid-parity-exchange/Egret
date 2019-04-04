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
        scale_ModelData_to_pu

def _get_uc_model(model_data, formulation_list):
    formulation = UCFormulation(*formulation_list)
    md = scale_ModelData_to_pu(model_data)
    return generate_model(md, formulation)

def create_tight_unit_commitment_model(model_data, network_constraints='power_balance_constraints'):
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
    return _get_uc_model(model_data, formulation_list)

def create_compact_unit_commitment_model(model_data, network_constraints='power_balance_constraints'):
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
    return _get_uc_model(model_data, formulation_list)

def create_KOW_unit_commitment_model(model_data, network_constraints='power_balance_constraints'):
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
    return _get_uc_model(model_data, formulation_list)

def create_ALS_unit_commitment_model(model_data, network_constraints='power_balance_constraints'):
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
    return _get_uc_model(model_data, formulation_list)

def create_MLR_unit_commitment_model(model_data, network_constraints='power_balance_constraints'):
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
    return _get_uc_model(model_data, formulation_list)

def create_random1_unit_commitment_model(model_data, network_constraints='power_balance_constraints'):
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
    return _get_uc_model(model_data, formulation_list)

def create_random2_unit_commitment_model(model_data, network_constraints='power_balance_constraints'):
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
    return _get_uc_model(model_data, formulation_list)
