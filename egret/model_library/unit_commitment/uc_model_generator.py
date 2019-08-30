#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
Helpers for constructing unit commitment models
'''

from egret.model_library.unit_commitment import \
        params, status_vars, power_vars, reserve_vars, \
        non_dispatchable_vars, generation_limits, \
        ramping_limits, production_costs, \
        uptime_downtime, startup_costs, \
        services, power_balance, reserve_requirement, \
        objective, fuel_supply, fuel_consumption

from collections import namedtuple
import pyomo.environ as pe

## tools for generating models

UCFormulation = namedtuple('UCFormulation',
                            ['status_vars',
                             'power_vars',
                             'reserve_vars',
                             'generation_limits',
                             'ramping_limits',
                             'production_costs',
                             'uptime_downtime',
                             'startup_costs',
                             'network_constraints',
                             ]
                            )

def generate_model( model_data, uc_formulation, relax_binaries=False, ptdf_options=None ):
    """
    returns a UC uc_formulation as an abstract model with the 
    components specified in a UCFormulation, with the option
    to relax the binary variables.

    Parameters
    ----------
    model_data : egret.data.ModelData
    uc_formulation : egret.model_components.model_generator.UCFormulation
        The named tuple with the specified formulation
    relax_binaries : bool, optional
        Relaxes all binary variables in the constructed model, resulting in a continuous problem.
        Default is False.
    ptdf_options : dict, optional
        Dictionary of options for ptdf transmission model

    Returns
    -------
        pyomo.environ.ConcreteModel : The unit commitment formulation specified with the data
                                      from model_data
    """

    return _generate_model( model_data, *_get_formulation_from_UCFormulation( uc_formulation ), relax_binaries , ptdf_options )

def _generate_model( model_data,
                    _status_vars,
                    _power_vars,
                    _reserve_vars,
                    _non_dispatchable_vars,
                    _generation_limits,
                    _ramping_limits,
                    _production_costs,
                    _uptime_downtime,
                    _startup_costs,
                    _power_balance, 
                    _reserve_requirement, 
                    _objective, 
                    _relax_binaries = False,
                    _ptdf_options = None,
                    ):

    
    model = pe.ConcreteModel()

    # hold on to this for later
    model.model_data = model_data
    
    model.name = "UnitCommitment"
    
    ## munge ptdf_options, if necessary
    if _power_balance in ['ptdf_power_flow']:
        import egret.common.lazy_ptdf_utils as lpu
        if _ptdf_options is None:
            _ptdf_options = dict()
        lpu.populate_default_ptdf_options(_ptdf_options)

        baseMVA = model_data.data['system']['baseMVA']
        lpu.check_and_scale_ptdf_options(_ptdf_options, baseMVA)

        model._ptdf_options = _ptdf_options

    ## enforece time 1 ramp rates
    model.enforce_t1_ramp_rates = True
    
    ## to relax binaries
    model.relax_binaries = _relax_binaries

    params.load_params(model, model_data)
    getattr(status_vars, _status_vars)(model)
    getattr(power_vars, _power_vars)(model)
    getattr(reserve_vars, _reserve_vars)(model)
    getattr(non_dispatchable_vars, _non_dispatchable_vars)(model)
    getattr(generation_limits, _generation_limits)(model)
    getattr(ramping_limits, _ramping_limits)(model)
    getattr(production_costs, _production_costs)(model)
    getattr(uptime_downtime, _uptime_downtime)(model)
    getattr(startup_costs, _startup_costs)(model)
    services.storage_services(model)
    services.ancillary_services(model)
    getattr(power_balance, _power_balance)(model)
    getattr(reserve_requirement, _reserve_requirement)(model)

    if 'fuel_supply' in model_data.data['elements'] and bool(model_data.data['elements']['fuel_supply']):
        fuel_consumption.fuel_consumption_model(model)
        fuel_supply.fuel_supply_model(model)

    getattr(objective, _objective)(model)

    return model

def _get_formulation_from_UCFormulation( uc_formulation ):
    return [  uc_formulation.status_vars,
              uc_formulation.power_vars,
              uc_formulation.reserve_vars,
              'file_non_dispatchable_vars',
              uc_formulation.generation_limits,
              uc_formulation.ramping_limits,
              uc_formulation.production_costs,
              uc_formulation.uptime_downtime,
              uc_formulation.startup_costs,
              uc_formulation.network_constraints,
              'MLR_reserve_constraints' if uc_formulation.reserve_vars in ['MLR_reserve_vars'] else 'CA_reserve_constraints',
              'basic_objective',
            ]
