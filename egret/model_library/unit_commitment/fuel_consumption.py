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

from .uc_utils import add_model_attr
from .status_vars import _is_relaxed

@add_model_attr('fuel_consumption', requires = {'data_loader': None,
                                                'status_vars': ['garver_3bin_vars','garver_2bin_vars','garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
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
def fuel_consumption_model(model):
    '''
    Defines extra variables and constraints for tracking fuel consumption
    for generators
    '''

    md = model.model_data

    relaxed = _is_relaxed(model)
    ## generator fuel consumption model
    thermal_gen_attrs = md.attributes(element_type='generator', generator_type='thermal')


    model.FuelConsumedCommitment = Var(model.FuelSupplyGenerators, model.TimePeriods, within=NonNegativeReals)
    model.FuelConsumedProduction = Var(model.FuelSupplyGenerators, model.TimePeriods, within=NonNegativeReals)

    def _fuel_consumed_function(m, g, t, i):
        return m.PowerGenerationPiecewiseFuelValues[g,t][i] * m.TimePeriodLengthHours

    def production_fuel_consumed_rule(m, g, t):
        if (g,t) in m.PiecewiseGeneratorTimeIndexSet:
            return sum( (_fuel_consumed_function(m,g,t,i+1) - _fuel_consumed_function(m,g,t,i))/(m.PowerGenerationPiecewisePoints[g,t][i+1] - m.PowerGenerationPiecewisePoints[g,t][i]) * m.PiecewiseProduction[g,t,i] for i in range(len(m.PowerGenerationPiecewisePoints[g,t])-1))
        elif (g,t) in m.LinearGeneratorTimeIndexSet:
            i=0
            return (_fuel_consumed_function(m,g,t,i+1) - _fuel_consumed_function(m,g,t,i))/(m.PowerGenerationPiecewisePoints[g,t][i+1] - m.PowerGenerationPiecewisePoints[g,t][i]) * m.PowerGeneratedAboveMinimum[g,t]
        else:
            return 0.

    model.ProductionFuelConsumed = Expression(model.FuelSupplyGenerators, model.TimePeriods, rule=production_fuel_consumed_rule)

    def production_fuel_consumed_constr(m,g,t):
        return m.FuelConsumedProduction[g,t] == m.ProductionFuelConsumed[g,t]
    model.ProductionFuelConsumedConstr = Constraint(model.FuelSupplyGenerators, model.TimePeriods, rule=production_fuel_consumed_constr)

    def _startup_fuel_consumed_function(m, g, i):
        return thermal_gen_attrs['startup_fuel'][g][i][1]

    model_startup_costs = model.startup_costs
    def startup_fuel_consumed_rule(m,g,t):
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

    def fuel_commitment_consumed_rule(m,g,t):
        return m.MinimumFuelConsumption[g,t]*m.TimePeriodLengthHours*m.UnitOn[g,t] 
    model.CommitmentFuelConsumed = Expression(model.FuelSupplyGenerators, model.TimePeriods, rule=fuel_commitment_consumed_rule)

    def commitment_fuel_consumed_constr(m,g,t):
        return m.FuelConsumedCommitment[g,t] == m.StartupFuelConsumed[g,t] + m.CommitmentFuelConsumed[g,t]
    model.FuelConsumedCommitmentConstr = Constraint(model.FuelSupplyGenerators, model.TimePeriods, rule=commitment_fuel_consumed_constr)

    ## end generator fuel consumption model

    ## dual fuel and fuel supply shared expression for fuel supply models

    model.PrimaryFuelConsumedCommitment = Var(model.DualFuelGenerators, model.TimePeriods, within=NonNegativeReals)
    model.PrimaryFuelConsumedProduction = Var(model.DualFuelGenerators, model.TimePeriods, within=NonNegativeReals)

    def primary_fuel_consumed(m, g, t):
        if g in m.DualFuelGenerators:
            return m.PrimaryFuelConsumedCommitment[g,t] + m.PrimaryFuelConsumedProduction[g,t]
        else:
            return m.FuelConsumedProduction[g,t] + m.FuelConsumedCommitment[g,t]
    model.PrimaryFuelConsumed = Expression(model.FuelSupplyGenerators, model.TimePeriods, rule=primary_fuel_consumed)

    ## if we don't have any dual fuel generators,
    ## we don't need to do anything else
    if not model.DualFuelGenerators:
        return

    ## DUAL FUEL GENERATORS

    ## load and verify some parameters
    dual_fuel_attrs = md.attributes(element_type='generator', generator_type='thermal', aux_fuel_capable=True)

    model.UnitSwitchOperating = Param(model.DualFuelGenerators, within=Boolean, default=False, initialize=dual_fuel_attrs.get('aux_fuel_online_switching', dict()))

    model.UnitFuelBlending = Param(model.DualFuelGenerators, within=Boolean, default=False, initialize=dual_fuel_attrs.get('aux_fuel_blending', dict()))

    def verify_dual_fuel_consistency(m, g):
        if value(m.UnitFuelBlending[g]) and not value(m.UnitSwitchOperating[g]):
            print("DATA ERROR: Dual fuel generators which are fuel blending capable must also be online fuel switching capable")
            print("Generator {} had aux_fuel_blending = True and aux_fuel_online_switching=False".format(g))
            assert(False)
    model.VerifyDualFuelBlendingSwitchingConsistency = BuildAction(model.DualFuelGenerators, rule=verify_dual_fuel_consistency)
    
    def verify_initial_fuel_defined(m, g):
        if not value(m.UnitSwitchOperating[g]) and value(m.UnitOnT0[g]):
            if 'aux_fuel_supply_initial' not in dual_fuel_attrs or g not in dual_fuel_attrs['aux_fuel_supply_initial']:
                print("DATA ERROR: Couldn't find initial fuel for dual fuel generator "+g)
                assert(False)
    model.VerifyUnitInitialFuelDefined = BuildAction(model.DualFuelGenerators, rule=verify_initial_fuel_defined)

    model.AuxiliaryFuelConsumedCommitment = Var(model.DualFuelGenerators, model.TimePeriods, within=NonNegativeReals)

    model.AuxiliaryFuelConsumedProduction = Var(model.DualFuelGenerators, model.TimePeriods, within=NonNegativeReals)

    def pri_aux_fuel_consumed_commitment_rule(m, g, t):
        return m.FuelConsumedCommitment[g,t] == m.PrimaryFuelConsumedCommitment[g,t] + m.AuxiliaryFuelConsumedCommitment[g,t]
    model.PriAuxFuelConsumedCommitmentConstr = Constraint(model.DualFuelGenerators, model.TimePeriods, rule=pri_aux_fuel_consumed_commitment_rule)

    def pri_aux_fuel_consumed_production_rule(m, g, t):
        return m.FuelConsumedProduction[g,t] == m.PrimaryFuelConsumedProduction[g,t] + m.AuxiliaryFuelConsumedProduction[g,t]
    model.PriAuxFuelConsumedProductionConstr = Constraint(model.DualFuelGenerators, model.TimePeriods, rule=pri_aux_fuel_consumed_production_rule)


    def auxiliary_fuel_consumed(m, g, t):
        return m.AuxiliaryFuelConsumedCommitment[g,t] + m.AuxiliaryFuelConsumedProduction[g,t]
    model.AuxiliaryFuelConsumed = Expression(model.DualFuelGenerators, model.TimePeriods, rule=auxiliary_fuel_consumed)

    model.PrimaryFuelCost = Param(model.DualFuelGenerators, initialize=dual_fuel_attrs['fuel_cost'], within=NonNegativeReals)

    model.AuxiliaryFuelCost = Param(model.DualFuelGenerators, initialize=dual_fuel_attrs['aux_fuel_cost'], within=NonNegativeReals)

    model.NonFuelNoLoadCost = Param(model.DualFuelGenerators, initialize=dual_fuel_attrs.get('non_fuel_no_load_cost', dict()), default=0., within=Reals)

    model.NonFuelStartupCost = Param(model.DualFuelGenerators, initialize=dual_fuel_attrs.get('non_fuel_startup_cost', dict()), default=0., within=Reals)

    def dual_fuel_startup_running_cost(m,g,t):
        return m.PrimaryFuelCost[g]*m.PrimaryFuelConsumedCommitment[g,t] \
                + m.AuxiliaryFuelCost[g]*m.AuxiliaryFuelConsumedCommitment[g,t] \
                + m.NonFuelNoLoadCost[g]*m.UnitOn[g,t] \
                + m.NonFuelStartupCost[g]*m.UnitStart[g,t]
    model.DualFuelCommitmentCost = Expression(model.DualFuelGenerators, model.TimePeriods, rule=dual_fuel_startup_running_cost)

    def dual_fuel_production_cost(m,g,t):
        return m.PrimaryFuelCost[g]*m.PrimaryFuelConsumedProduction[g,t] \
                + m.AuxiliaryFuelCost[g]*m.AuxiliaryFuelConsumedProduction[g,t]
    model.DualFuelProductionCost = Expression(model.DualFuelGenerators, model.TimePeriods, rule=dual_fuel_production_cost)

    ## SINGLE-FIRE DUAL-FUEL UNITS 
    def init_single_fire(m):
        for g in m.DualFuelGenerators:
            if value(m.UnitFuelBlending[g]):
                continue
            else:
                yield g
    model.SingleFireDualFuelGenerators = Set(within=model.DualFuelGenerators, initialize=init_single_fire)

    model.UnitOnPriFuel = Var(model.SingleFireDualFuelGenerators, model.TimePeriods, within=Binary)
    model.UnitOnAuxFuel = Var(model.SingleFireDualFuelGenerators, model.TimePeriods, within=Binary)

    def single_fire_rule(m, g, t):
        return m.UnitOn[g,t] == m.UnitOnPriFuel[g,t] + m.UnitOnAuxFuel[g,t]
    model.UnitOnLink = Constraint(model.SingleFireDualFuelGenerators, model.TimePeriods, rule=single_fire_rule)

    def init_fuel_ub(m,g):
        return (thermal_gen_attrs['p_fuel'][g]['values'][-1][1])*m.TimePeriodLengthHours + thermal_gen_attrs['startup_fuel'][g][-1][1]
    model.FuelConsumedUB = Param(model.SingleFireDualFuelGenerators, initialize=init_fuel_ub)

    def enforce_single_fire_primary(m, g, t):
        return m.PrimaryFuelConsumed[g,t] <= m.FuelConsumedUB[g]*m.UnitOnPriFuel[g,t]
    model.PrimaryFuelConstr = Constraint(model.SingleFireDualFuelGenerators, model.TimePeriods, rule=enforce_single_fire_primary)

    def enforce_single_fire_auxiliary(m, g, t):
        return m.AuxiliaryFuelConsumed[g,t] <= m.FuelConsumedUB[g]*m.UnitOnAuxFuel[g,t]
    model.AuxiliaryFuelConstr = Constraint(model.SingleFireDualFuelGenerators, model.TimePeriods, rule=enforce_single_fire_auxiliary)


    ### DAUL-FUEL UNITS WHICH CAN ONLY FUEL-SWITCH OFFLINE

    def init_offline_switching(m):
        for g in m.SingleFireDualFuelGenerators:
            if value(m.UnitSwitchOperating[g]):
                continue
            else:
                yield g
    model.OfflineSwitchingDualFuelGenerators = Set(within=model.SingleFireDualFuelGenerators, initialize=init_offline_switching)

    model.UnitStartPriFuel = Var(model.OfflineSwitchingDualFuelGenerators, model.TimePeriods, within=Binary)
    model.UnitStartAuxFuel = Var(model.OfflineSwitchingDualFuelGenerators, model.TimePeriods, within=Binary)

    model.UnitStopPriFuel = Var(model.OfflineSwitchingDualFuelGenerators, model.TimePeriods, within=Binary)
    model.UnitStopAuxFuel = Var(model.OfflineSwitchingDualFuelGenerators, model.TimePeriods, within=Binary)

    def startup_fire_rule(m,g,t):
        return m.UnitStart[g,t] == m.UnitStartPriFuel[g,t] + m.UnitStartAuxFuel[g,t]
    model.UnitStartLink = Constraint(model.OfflineSwitchingDualFuelGenerators, model.TimePeriods, rule=startup_fire_rule)

    def shutdown_fire_rule(m,g,t):
        return m.UnitStop[g,t] == m.UnitStopPriFuel[g,t] + m.UnitStopAuxFuel[g,t]
    model.UnitStopLink = Constraint(model.OfflineSwitchingDualFuelGenerators, model.TimePeriods, rule=shutdown_fire_rule)

    ## These enforce that if the generator starts with a certain fuel type it
    ## must maintain it until shutdown.
    def pri_logical_constr(m, g, t):
        if t==value(m.InitialTime):
            if value(m.UnitOnT0[g]):
                pri_fuel_init = 1 - dual_fuel_attrs['aux_fuel_supply_initial'][g]
                return m.UnitOnPriFuel[g,t] - pri_fuel_init == m.UnitStartPriFuel[g,t] - m.UnitStopPriFuel[g,t]
            else:
                return Constraint.Skip
        else:
            return m.UnitOnPriFuel[g,t] - m.UnitOnPriFuel[g,t-1] == m.UnitStartPriFuel[g,t] - m.UnitStopPriFuel[g,t]
    model.PriLogicalConst = Constraint(model.OfflineSwitchingDualFuelGenerators, model.TimePeriods, rule=pri_logical_constr)

    def aux_logical_constr(m, g, t):
        if t==value(m.InitialTime):
            if value(m.UnitOnT0[g]):
                aux_fuel_init = dual_fuel_attrs['aux_fuel_supply_initial'][g]
                return m.UnitOnAuxFuel[g,t] - aux_fuel_init == m.UnitStartAuxFuel[g,t] - m.UnitStopAuxFuel[g,t]
            else:
                return Constraint.Skip
        else:
            return m.UnitOnAuxFuel[g,t] - m.UnitOnAuxFuel[g,t-1] == m.UnitStartAuxFuel[g,t] - m.UnitStopAuxFuel[g,t]
    model.AuxLogicalConst = Constraint(model.OfflineSwitchingDualFuelGenerators, model.TimePeriods, rule=aux_logical_constr)

