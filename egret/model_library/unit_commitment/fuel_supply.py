#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

from pyomo.environ import *
import math

from .uc_utils import add_model_attr, build_uc_time_mapping
from .status_vars import _is_relaxed

## For now, this function will depend on having one of the production costs with
## PiecewiseProduction variables, and start-up costs which use auxillary variables
@add_model_attr('fuel_supply', requires = {'data_loader': None,
                                           'status_vars': None,
                                           'power_vars': None,
                                           'production_costs': ['CA_production_costs',
                                                                'wu_production_costs',
                                                                'KOW_production_costs',
                                                                'KOW_production_costs_tightened',
                                                                'KOW_production_costs_super_tight',
                                                                ],
                                           'startup_costs' : ['KOW_startup_costs',
                                                              'MLR_startup_costs',
                                                              'MLR_startup_costs2',
                                                             ],
                                           })
def fuel_supply_model(model):
    '''
    Defines a fuel supply model for generators. For now, just considers
    'instantaneous' supply, i.e., like NG
    '''
    md = model.model_data

    system = md.data['system']
    time_keys = system['time_indices']
    TimeMapper = build_uc_time_mapping(time_keys)


    ## generator fuel consumption model
    thermal_gen_attrs = md.attributes(element_type='generator', generator_type='thermal')

    model.FuelConsumed = Var(model.FuelSupplyGenerators, model.TimePeriods, within=NonNegativeReals)
    model.AuxFuelConsumed = Var(model.AuxFuelSupplyGenerators, model.TimePeriods, within=NonNegativeReals)

    def _fuel_consumed_function(m, g, i):
        return thermal_gen_attrs['p_fuel'][g]['values'][i][1]*m.TimePeriodLengthHours
    def _aux_fuel_consumed_function(m, g, i):
        return thermal_gen_attrs['aux_p_fuel'][g]['values'][i][1]*m.TimePeriodLengthHours

    def production_fuel_consumed_rule(m, g, t):
        if g in m.DualFuelGenerators:
            p_cost_vals = thermal_gen_attrs['p_cost'][g]['values']
            return sum((_fuel_consumed_function(m,g,i+1) - _fuel_consumed_function(m,g,i))/(p_cost_vals[i+1][0] - p_cost_vals[i][0])*m.PriPiecewiseProduction[g,t,i] \
                     for i in range(m.NumProductionSegments[g,'p_cost']))
        else:    
            if (g,t) in m.PiecewiseGeneratorTimeIndexSet:
                return sum( (_fuel_consumed_function(m,g,i+1) - _fuel_consumed_function(m,g,i))/(m.PowerGenerationPiecewisePoints[g,t][i+1] - m.PowerGenerationPiecewisePoints[g,t][i]) * m.PiecewiseProduction[g,t,i] for i in range(len(m.PowerGenerationPiecewisePoints[g,t])-1))
            elif (g,t) in m.LinearGeneratorTimeIndexSet:
                i=0
                return (_fuel_consumed_function(m,g,i+1) - _fuel_consumed_function(m,g,i))/(m.PowerGenerationPiecewisePoints[g,t][i+1] - m.PowerGenerationPiecewisePoints[g,t][i]) * m.PowerGeneratedAboveMinimum[g,t]
            else:
                return 0.
    model.ProductionFuelConsumed = Expression(model.FuelSupplyGenerators, model.TimePeriods, rule=production_fuel_consumed_rule)

    def aux_production_fuel_consumed_rule(m, g, t):
        assert (g in m.DualFuelGenerators)
        aux_cost_vals = thermal_gen_attrs['aux_p_cost'][g]['values']
        return sum((_aux_fuel_consumed_function(m,g,i+1) - _aux_fuel_consumed_function(m,g,i))/(aux_cost_vals[i+1][0] - aux_cost_vals[i][0])*m.AuxPiecewiseProduction[g,t,i] \
                     for i in range(m.NumProductionSegments[g,'aux_p_cost']))
    model.AuxProductionFuelConsumed = Expression(model.AuxFuelSupplyGenerators, model.TimePeriods, rule=aux_production_fuel_consumed_rule)

    def _startup_fuel_consumed_function(m, g, i):
        return thermal_gen_attrs['startup_fuel'][g][i][1]
    def _aux_startup_fuel_consumed_function(m, g, i):
        return thermal_gen_attrs['aux_startup_fuel'][g][i][1]

    model_startup_costs = model.startup_costs
    def startup_fuel_consumed_rule(m,g,t):
        if g in m.DualFuelGenerators:
            if 'startup_fuel' in thermal_gen_attrs and g in thermal_gen_attrs['startup_fuel']:
                last_fuel_consumed = _startup_fuel_consumed_function(m,g,-1)
                return last_fuel_consumed*m.UnitStartPriFuel[g,t] + \
                            sum( (_startup_fuel_consumed_function(m,g,s-1) - last_fuel_consumed) * \
                            sum( m.PriStartupIndicator[g,tp,t] for tp in m.ValidShutdownTimePeriods[g] \
                              if (list(m.ScaledStartupLags[g])[s-1] <= t - tp < (list(m.ScaledStartupLags[g])[s])) ) \
                            for s in m.StartupCostIndices[g] if s < len(m.StartupCostIndices[g]))
        else:
            if 'startup_fuel' in thermal_gen_attrs and g in thermal_gen_attrs['startup_fuel']:
                if model_startup_costs == 'KOW_startup_costs':
                    last_fuel_consumed = _startup_fuel_consumed_function(m,g,-1)
                    return last_fuel_consumed*m.UnitStart[g,t] + \
                                sum( (_startup_fuel_consumed_function(m,g,s-1) - last_fuel_consumed) * \
                                sum( m.StartupIndicator[g,tp,t] for tp in m.ValidShutdownTimePeriods[g] \
                                  if (list(m.ScaledStartupLags[g])[s-1] <= t - tp < (list(m.ScaledStartupLags[g])[s])) ) \
                                for s in m.StartupCostIndices[g] if s < len(m.StartupCostIndices[g]))
                elif model_startup_costs == 'MLR_startup_costs':
                    return sum(_startup_fuel_consumed_function(m,g,s-1)*m.delta[g,s,t] for s in m.StartupCostIndices[g])
                elif model_startup_costs == 'MLR_startup_costs2':
                    last_fuel_consumed = _startup_fuel_consumed_function(m,g,-1)
                    return last_fuel_consumed*m.UnitStart[g,t] + \
                               sum((_startup_fuel_consumed_function(m,g,s-1) - last_fuel_consumed)*m.delta[g,s,t] 
                                    for s in m.StartupCostIndices[g] if s < len(m.StartupCostIndices[g]))
                else:
                    raise Exception('fuel_supply model requires one of KOW_startup_costs, MLR_startup_costs, or MLR_startup_costs2')
            else:
                return 0.
    model.StartupFuelConsumed = Expression(model.FuelSupplyGenerators, model.TimePeriods, rule=startup_fuel_consumed_rule)

    def aux_startup_fuel_consumed_rule(m,g,t):
        assert (g in m.DualFuelGenerators)
        if 'aux_startup_fuel' in thermal_gen_attrs and g in thermal_gen_attrs['aux_startup_fuel']:
            last_fuel_consumed = _aux_startup_fuel_consumed_function(m,g,-1)
            return last_fuel_consumed*m.UnitStartAuxFuel[g,t] + \
                        sum( (_aux_startup_fuel_consumed_function(m,g,s-1) - last_fuel_consumed) * \
                        sum( m.AuxStartupIndicator[g,tp,t] for tp in m.ValidShutdownTimePeriods[g] \
                          if (list(m.AuxScaledStartupLags[g])[s-1] <= t - tp < (list(m.AuxScaledStartupLags[g])[s])) ) \
                        for s in m.AuxStartupCostIndices[g] if s < len(m.AuxStartupCostIndices[g]))
        else:
            return 0.
    model.AuxStartupFuelConsumed = Expression(model.AuxFuelSupplyGenerators, model.TimePeriods, rule=startup_fuel_consumed_rule)

    def fuel_consumed_constr_rule(m, g, t):
        ######  total fuel consumed   ==   fuel consumed above minimum + fuel consumed for being on and operating at p_min + fuel consumed for startup
        if g in m.DualFuelGenerators:
            return m.FuelConsumed[g,t] == m.ProductionFuelConsumed[g,t] + _fuel_consumed_function(m,g,0)*m.UnitOnPriFuel[g,t] + m.StartupFuelConsumed[g,t]
        else:
            return m.FuelConsumed[g,t] == m.ProductionFuelConsumed[g,t] + _fuel_consumed_function(m,g,0)*m.UnitOn[g,t] + m.StartupFuelConsumed[g,t]
    model.FuelConsumedConstr = Constraint(model.FuelSupplyGenerators, model.TimePeriods, rule=fuel_consumed_constr_rule)

    def aux_fuel_consumed_constr_rule(m, g, t):
        ######  total fuel consumed   ==   fuel consumed above minimum + fuel consumed for being on and operating at p_min + fuel consumed for startup
        assert (g in m.DualFuelGenerators)
        return m.AuxFuelConsumed[g,t] == m.AuxProductionFuelConsumed[g,t] + _aux_fuel_consumed_function(m,g,0)*m.UnitOnAuxFuel[g,t] + m.AuxStartupFuelConsumed[g,t]
    model.AuxFuelConsumedConstr = Constraint(model.FuelSupplyGenerators, model.TimePeriods, rule=fuel_consumed_constr_rule)
    ## end generator fuel consumption model

    ## instantaneous fuel supply model
    inst_fuel_supply_attrs = md.attributes(element_type='fuel_supply', fuel_supply_type='instantaneous')

    model.InstantaneousFuelSupplies = Set(initialize=inst_fuel_supply_attrs['names'])

    def gens_using_inst_fuel_supply_init(m, f):
        for g in m.FuelSupplyGenerators:
            if thermal_gen_attrs['fuel_supply'][g] == f:
                yield g
    model.ThermalGeneratorsUsingInstFuelSupply = Set(model.InstantaneousFuelSupplies, initialize=gens_using_inst_fuel_supply_init,)

    def gens_using_inst_fuel_supply_as_aux_init(m, f):
        for g in m.AuxFuelSupplyGenerators:
            if thermal_gen_attrs['aux_fuel_supply'][g] == f:
                yield g
    model.ThermalGeneratorsUsingInstFuelSupplyAsAux = Set(model.InstantaneousFuelSupplies, initialize=gens_using_inst_fuel_supply_as_aux_init)

    model.InstFuelSupply = Param(model.InstantaneousFuelSupplies, model.TimePeriods, initialize=TimeMapper(inst_fuel_supply_attrs['fuel_available']))

    def total_fuel_consumed_expr(m, f, t):
        return sum(m.FuelConsumed[g,t] for g in m.ThermalGeneratorsUsingInstFuelSupply[f]) \
                + sum(m.AuxFuelConsumed[g,t] for g in m.ThermalGeneratorsUsingInstFuelSupplyAsAux[f])
    model.TotalFuelConsumedAtInstFuelSupply = Expression(model.InstantaneousFuelSupplies, model.TimePeriods, rule=total_fuel_consumed_expr)

    def total_fuel_consumed_rule(m, f, t):
        if m.ThermalGeneratorsUsingInstFuelSupply[f]:
            return m.TotalFuelConsumedAtInstFuelSupply[f,t] <= m.InstFuelSupply[f,t]
        else:
            if t == m.TimePeriods.first():
                print('WARNING: no generators attached to fuel_supply {}'.format(f))
            return Constraint.Feasible
    model.FuelLimitConstr = Constraint(model.InstantaneousFuelSupplies, model.TimePeriods, rule=total_fuel_consumed_rule)
    ## end instantaneous fuel supply model
