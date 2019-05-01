#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## file for all the ancillary services
from pyomo.environ import *
import math

from .uc_utils import add_model_attr 

@add_model_attr('storage_service', requires = {'data_loader': None,
                                            })
def storage_services(model):
    '''
    Defines a storage component
    '''

    ##############################
    # Storage decision variables #
    ##############################

    # binary variables for storage (input/output are semicontinuous)
    # NOTE: two variables for this because presumably the battery can be "turned off",
    #       neither inputing or outputing. If we only have one variable, then 
    #       if the minimum charge/discharge rates are not zero, we cannot model
    #       this condition.
    model.InputStorage = Var(model.Storage, model.TimePeriods, within=Binary)
    model.OutputStorage = Var(model.Storage, model.TimePeriods, within=Binary)

    def input_output_complementarity_rule(m,s,t):
        return m.InputStorage[s,t] + m.OutputStorage[s,t] <= 1
    model.InputOutputComplementarity = Constraint(model.Storage, model.TimePeriods, rule=input_output_complementarity_rule)

    # amount of output power of each storage unit, at each time period, on the grid side
    def power_output_storage_bounds_rule(m, s, t):
        return (0, m.MaximumPowerOutputStorage[s])
    model.PowerOutputStorage = Var(model.Storage, model.TimePeriods, within=NonNegativeReals, bounds=power_output_storage_bounds_rule)

    # amount of input power of each storage unit, at each time period, on the grid side
    def power_input_storage_bounds_rule(m, s, t):
        return (0, m.MaximumPowerInputStorage[s])
    model.PowerInputStorage = Var(model.Storage, model.TimePeriods, within=NonNegativeReals, bounds=power_input_storage_bounds_rule)

    # state of charge of each storage unit, at each time period
    model.SocStorage = Var(model.Storage, model.TimePeriods, within=PercentFraction)

    def min_soc_rule(model, m, t):
        return model.SocStorage[m,t] >= model.MinimumSocStorage[m]
    model.SocMinimum = Constraint(model.Storage, model.TimePeriods, rule=min_soc_rule)

    #######################################
    # energy storage bounding constraints #
    #######################################
    # NOTE: The expressions below are what we really want - however, due to a pyomo design feature, we have to split it into two constraints:
    # m.MinimumPowerInputStorage[g] * m.InputStorage[g, t] <= m.StoragePowerInput[g,t] <= m.MaximumPowerInputStorage[g] * m.InputStorage[g, t]
    # m.MinimumPowerOutputStorage[g] * m.OutputStorage[g, t] <= m.StoragePowerOutput[g,t] <= m.MaximumPowerOutputStorage[g] * m.OutputStorage[g, t]

    def enforce_storage_input_limits_rule_part_a(m, s, t):
        return m.MinimumPowerInputStorage[s] * (m.InputStorage[s, t]) <= m.PowerInputStorage[s,t]

    model.EnforceStorageInputLimitsPartA = Constraint(model.Storage, model.TimePeriods, rule=enforce_storage_input_limits_rule_part_a)

    def enforce_storage_input_limits_rule_part_b(m, s, t):
        return m.PowerInputStorage[s,t] <= m.MaximumPowerInputStorage[s] * (m.InputStorage[s, t])

    model.EnforceStorageInputLimitsPartB = Constraint(model.Storage, model.TimePeriods, rule=enforce_storage_input_limits_rule_part_b)

    def enforce_storage_output_limits_rule_part_a(m, s, t):
        return m.MinimumPowerOutputStorage[s] * m.OutputStorage[s, t] <= m.PowerOutputStorage[s,t]

    model.EnforceStorageOutputLimitsPartA = Constraint(model.Storage, model.TimePeriods, rule=enforce_storage_output_limits_rule_part_a)

    def enforce_storage_output_limits_rule_part_b(m, s, t):
        return m.PowerOutputStorage[s,t] <= m.MaximumPowerOutputStorage[s] * m.OutputStorage[s, t]

    model.EnforceStorageOutputLimitsPartB = Constraint(model.Storage, model.TimePeriods, rule=enforce_storage_output_limits_rule_part_b)

    #####################################
    # energy storage ramping contraints #
    #####################################

    def enforce_ramp_up_rates_power_output_storage_rule(m, s, t):
        if t == m.InitialTime:
            return m.PowerOutputStorage[s, t] <= m.StoragePowerOutputOnT0[s] + m.ScaledNominalRampUpLimitStorageOutput[s]
        else:
            return m.PowerOutputStorage[s, t] <= m.PowerOutputStorage[s, t-1] + m.ScaledNominalRampUpLimitStorageOutput[s]

    model.EnforceStorageOutputRampUpRates = Constraint(model.Storage, model.TimePeriods, rule=enforce_ramp_up_rates_power_output_storage_rule)

    def enforce_ramp_down_rates_power_output_storage_rule(m, s, t):
        if t == m.InitialTime:
            return m.PowerOutputStorage[s, t] >= m.StoragePowerOutputOnT0[s] - m.ScaledNominalRampDownLimitStorageOutput[s]
        else:
            return m.PowerOutputStorage[s, t] >= m.PowerOutputStorage[s, t-1] - m.ScaledNominalRampDownLimitStorageOutput[s]

    model.EnforceStorageOutputRampDownRates = Constraint(model.Storage, model.TimePeriods, rule=enforce_ramp_down_rates_power_output_storage_rule)

    def enforce_ramp_up_rates_power_input_storage_rule(m, s, t):
        if t == m.InitialTime:
            return m.PowerInputStorage[s, t] <= m.StoragePowerInputOnT0[s] + m.ScaledNominalRampUpLimitStorageInput[s]
        else:
            return m.PowerInputStorage[s, t] <= m.PowerInputStorage[s, t-1] + m.ScaledNominalRampUpLimitStorageInput[s]

    model.EnforceStorageInputRampUpRates = Constraint(model.Storage, model.TimePeriods, rule=enforce_ramp_up_rates_power_input_storage_rule)

    def enforce_ramp_down_rates_power_input_storage_rule(m, s, t):
        if t == m.InitialTime:
            return m.PowerInputStorage[s, t] >= m.StoragePowerInputOnT0[s] - m.ScaledNominalRampDownLimitStorageInput[s]
        else:
            return m.PowerInputStorage[s, t] >= m.PowerInputStorage[s, t-1] - m.ScaledNominalRampDownLimitStorageInput[s]

    model.EnforceStorageInputRampDownRates = Constraint(model.Storage, model.TimePeriods, rule=enforce_ramp_down_rates_power_input_storage_rule)

    ##########################################
    # storage energy conservation constraint #
    ##########################################

    def energy_conservation_rule(m, s, t):
        # storage s, time t
        if t == m.InitialTime:
            return m.SocStorage[s, t] == m.StorageSocOnT0[s]  + \
                (-m.PowerOutputStorage[s, t]/m.OutputEfficiencyEnergy[s] + m.PowerInputStorage[s,t]*m.InputEfficiencyEnergy[s])*m.TimePeriodLengthHours/m.MaximumEnergyStorage[s]
        else:
            return m.SocStorage[s, t] == m.SocStorage[s, t-1]*m.ScaledRetentionRate[s]  + \
                (-m.PowerOutputStorage[s, t]/m.OutputEfficiencyEnergy[s] + m.PowerInputStorage[s,t]*m.InputEfficiencyEnergy[s])*m.TimePeriodLengthHours/m.MaximumEnergyStorage[s]
    model.EnergyConservation = Constraint(model.Storage, model.TimePeriods, rule=energy_conservation_rule)

    ##################################
    # storage end-point constraints  #
    ##################################

    def storage_end_point_soc_rule(m, s):
        # storage s, last time period
        return m.SocStorage[s, value(m.NumTimePeriods)] == m.EndPointSocStorage[s]
    model.EnforceEndPointSocStorage = Constraint(model.Storage, rule=storage_end_point_soc_rule)

    def storage_cost_rule(m, s, t):
        return m.ChargeCost[s]*m.PowerInputStorage[s,t]*m.TimePeriodLengthHours + \
                m.DischargeCost[s]*m.PowerOutputStorage[s,t]*m.TimePeriodLengthHours
    model.StorageCost = Expression(model.Storage, model.TimePeriods, rule=storage_cost_rule)

    return
## end storage_services


## TODO: add hooks into model_data
## TODO: add hooks for ThermalGeneratorForcedOutage
## NEW: general ancillary service function. These cannot be separated
##      because when multilple services are active they have interdependent constraints
##      (mostly involving capacity and ramping).
## NOTE: when moving to a Real-Time market, ramping limits need to be also considered
##       It seems MISO did not in its DA as of 2009 [1], but definitely does in RT as of 2016 [2].
##       This model is, however, based more on CAISO, which has ramping limits for all markets,[2],[3, p. 2-36].
##
##       [1] Ma, Xingwang, Haili Song, Mingguo Hong, Jie Wan, Yonghong Chen, and Eugene Zak.
##           "The security-constrained commitment and dispatch for Midwest ISO day-ahead co-
##           optimized energy and ancillary service market." In Power & Energy Society 
##           General Meeting, 2009. PES'09. IEEE, pp. 1-8. IEEE, 2009.
##       [2] Wang, Qin, and Bri-Mathias Hodge. "Enhancing power system operational flexibility 
##           with flexible ramping products: A review." IEEE Transactions on Industrial Informatics
##           13, no. NREL/JA-5D00-67471 (2017).
##       [3] Califonia ISO. Technical Bulletin 2009-06-05: Market Optimization Details. Revised Nov 19, 2009.

@add_model_attr('ancillary_service', requires = {'data_loader': None,
                                                  'status_vars': ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars'],
                                                  'power_vars': None,
                                                 })
def ancillary_services(model):
    '''
    Defines ancillary services: regulation, spinning reserve, nonspinning reserve, operational reserve, flexible ramp
    ## NOTE: As in most markets, the value of ancillary services from high to low is regulation, spinning reserve, nonspinning reserve, and operational reserve.
    ##       We allow for a higher-quality ancillary service to be subtituted for a lower-quality one
    ##       Flexible ramp is treated differently, again as it is in most markets. There is no bid for flexible ramp, and it is priced at opportunity cost
    '''

    ## begin regulation

    #################################
    # Regulation ancillary services #
    #################################
    
    # Regulation
    model.RegulationProvider = Param(model.ThermalGenerators, within=Binary, default=0) #indicates if a unit is offering regulation
    model.RegulationMinutes = Param(within=PositiveReals, default=5.)
    
    # When units are selected for regulation, their limits are bounded by the RegulationHighLimit and RegulationLowLimit
    # I'll refer to it as the "regulation band"
    def regulation_high_limit_validator(m, v, g):
        return v <= value(m.MaximumPowerOutput[g])
    def regulation_high_limit_init(m, g):
        return value(m.MaximumPowerOutput[g])   
    model.RegulationHighLimit = Param(model.ThermalGenerators, within=NonNegativeReals, validate=regulation_high_limit_validator, initialize=regulation_high_limit_init)
    
    def calculate_max_power_minus_reg_high_limit_rule(m, g):
        return m.MaximumPowerOutput[g] - m.RegulationHighLimit[g]
    model.MaxPowerOutputMinusRegHighLimit = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=calculate_max_power_minus_reg_high_limit_rule)
    
    def regulation_low_limit_validator(m, v, g):
        return (v <= value(m.RegulationHighLimit[g]) and v >= value(m.MinimumPowerOutput[g]))
    def regulation_low_limit_init(m, g):
        return value(m.MinimumPowerOutput[g])
    model.RegulationLowLimit = Param(model.ThermalGenerators, within=NonNegativeReals, validate=regulation_low_limit_validator, initialize=regulation_low_limit_init)
    
    # Regulation capacity is calculated as the max of "regulation band" and 5*AutomaticResponseRate
    model.AutomaticResponseRate = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0)
    
    def calculate_regulation_capability_rule(m, g):
        temp1 = 5 * m.AutomaticResponseRate[g]
        temp2 = (m.RegulationHighLimit[g] - m.RegulationLowLimit[g])/2
        if temp1 > temp2:
            return temp2
        else:
            return temp1
    
    model.RegulationUpCapability = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=calculate_regulation_capability_rule, default=0.0)
    model.RegulationDnCapability = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=calculate_regulation_capability_rule, default=0.0)

    model.ZonalRegulationUpRequirement = Param(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, default=0.0)
    model.SystemRegulationUpRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0)

    model.ZonalRegulationDnRequirement = Param(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, default=0.0)
    model.SystemRegulationDnRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0)

    model.RegulationPenalty = Param(within=NonNegativeReals, mutable=True)

    def set_regulation_penalty(m):
        if m.RegulationPenalty(exception=False) is None:
            m.RegulationPenalty = value(m.LoadMismatchPenalty+m.ReserveShortfallPenalty)/2.
    model.SetRegulationPenalty = BuildAction(rule=set_regulation_penalty)

    def zonal_up_bounds(m, rz, t):
        return (0, m.ZonalRegulationUpRequirement[zt,t])
    model.ZonalRegulationUpShortfall = Var(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_up_bounds)
    def zonal_dn_bounds(m, rz, t):
        return (0, m.ZonalRegulationDnRequirement[zt,t])
    model.ZonalRegulationDnShortfall = Var(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_dn_bounds)

    def system_up_bounds(m, t):
        return (0, m.SystemRegulationUpRequirement[t])
    model.SystemRegulationUpShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_up_bounds)
    def system_dn_bounds(m, t):
        return (0, m.SystemRegulationDnRequirement[t])
    model.SystemRegulationDnShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_dn_bounds)
    
    # regulation cost for
    model.RegulationOfferFixedCost = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0)
    model.RegulationOfferMarginalCost = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0)

    # constraint whether a unit can provide regulation
    def reg_on_bounds(m,g,t):
        return (0, m.RegulationProvider[g])
    model.RegulationOn = Var(model.ThermalGenerators, model.TimePeriods, within=Binary, bounds=reg_on_bounds)

    def reg_up_bounds(m,g,t):
        return (0, m.RegulationUpCapability[g])
    def reg_dn_bounds(m,g,t):
        return (0, m.RegulationDnCapability[g])
    model.RegulationReserveUp = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=reg_up_bounds)
    model.RegulationReserveDn = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=reg_dn_bounds)

    # a generator can provide regulation only when it's on
    def provide_regulation_when_unit_on_rule(m, g, t):
        return m.RegulationOn[g, t] <= m.UnitOn[g, t]
    model.EnforceRegulationOnWhenUnitOn = Constraint(model.ThermalGenerators, model.TimePeriods, rule=provide_regulation_when_unit_on_rule)

    def reg_up_rule(m,g,t):
        reg_up_limit = min(value(m.RegulationUpCapability[g]), value(m.NominalRampUpLimit[g]/60.*m.RegulationMinutes))
        return m.RegulationReserveUp[g,t] <= reg_up_limit*m.UnitOn[g,t]
    model.EnforceRegulationUpBound = Constraint(model.ThermalGenerators, model.TimePeriods, rule=reg_up_rule)

    def reg_dn_rule(m,g,t):
        reg_dn_limit = min(value(m.RegulationDnCapability[g]), value(m.NominalRampDownLimit[g]/60.*m.RegulationMinutes))
        return m.RegulationReserveDn[g,t] <= reg_dn_limit*m.UnitOn[g,t]
    model.EnforceRegulationDnBound = Constraint(model.ThermalGenerators, model.TimePeriods, rule=reg_dn_rule)

    def zonal_reg_up_provided(m,rz,t):
        return sum(m.RegulationReserveUp[g,t] for g in m.ThermalGeneratorsInReserveZone[rz]) + m.ZonalRegulationUpShortfall[rz,t] 
    model.ZonalRegulationUpProvided = Expression(model.ReserveZones, model.TimePeriods, rule=zonal_reg_up_provided)

    def enforce_zonal_reg_up_requirement_rule(m, rz, t):
        return  m.ZonalRegulationUpProvided[rz,t] >= m.ZonalRegulationUpRequirement[rz,t]
    model.EnforceZonalRegulationUpRequirements = Constraint(model.ReserveZones, model.TimePeriods, rule=enforce_zonal_reg_up_requirement_rule)

    def enforce_zonal_reg_dn_requirement_rule(m, rz, t):
        return sum(m.RegulationReserveDn[g,t] for g in m.ThermalGeneratorsInReserveZone[rz]) + \
                m.ZonalRegulationDnShortfall[rz,t] >= m.ZonalRegulationDnRequirement[rz,t]
    model.EnforceZonalRegulationDnRequirements = Constraint(model.ReserveZones, model.TimePeriods, rule=enforce_zonal_reg_dn_requirement_rule)

    ## NOTE: making sure not to double count the shortfall
    def system_reg_up_provided(m,t):
        return sum(m.RegulationReserveUp[g,t] for g in m.ThermalGenerators) + \
                m.SystemRegulationUpShortfall[t] + sum(m.ZonalRegulationUpShortfall[rz,t] for rz in m.ReserveZones) 
    model.SystemRegulationUpProvided = Expression(model.TimePeriods, rule=system_reg_up_provided)

    def enforce_system_regulation_up_requirement_rule(m, t):
        return m.SystemRegulationUpProvided[t] >= m.SystemRegulationUpRequirement[t]
    model.EnforceSystemRegulationUpRequirement = Constraint(model.TimePeriods, rule=enforce_system_regulation_up_requirement_rule)

    def enforce_system_regulation_dn_requirement_rule(m, t):
        return sum(m.RegulationReserveDn[g,t] for g in m.ThermalGenerators) + \
                m.SystemRegulationDnShortfall[t] + sum(m.ZonalRegulationDnShortfall[rz,t] for rz in m.ReserveZones) \
                >= m.SystemRegulationDnRequirement[t]
    model.EnforceSystemRegulationDnRequirement = Constraint(model.TimePeriods, rule=enforce_system_regulation_dn_requirement_rule)

    def regulation_cost_commitment(m,g,t):
        return m.RegulationOfferFixedCost[g] * m.RegulationOn[g, t]*m.TimePeriodLengthHours
    model.RegulationCostCommitment = Expression(model.ThermalGenerators, model.TimePeriods, rule=regulation_cost_commitment)

    def regulation_cost_generation(m,g,t):
        return m.RegulationOfferMarginalCost[g]*m.TimePeriodLengthHours*(m.RegulationReserveUp[g,t] + m.RegulationReserveDn[g,t])
    model.RegulationCostGeneration = Expression(model.ThermalGenerators, model.TimePeriods, rule=regulation_cost_generation)

    def regulation_cost_slacks(m,t):
        return m.TimePeriodLengthHours*m.RegulationPenalty*(
                        m.SystemRegulationUpShortfall[t] + m.SystemRegulationDnShortfall[t] \
                      + sum(m.ZonalRegulationUpShortfall[rz,t] for rz in m.ReserveZones) \
                      + sum(m.ZonalRegulationDnShortfall[rz,t] for rz in m.ReserveZones) \
                      )
    model.RegulationCostPenalty = Expression(model.TimePeriods, rule=regulation_cost_slacks)

    ## end regulation_services


    ## begin spinning reserve

    # spinning reserve response time
    model.SpinningReserveMinutes = Param(within=PositiveReals, default=10.) # in minutes, varies among ISOs

    # limit,  cost of spinning reserves
    def get_spin_bid(m,g):
        return m.MaximumPowerOutput[g] - m.MinimumPowerOutput[g]
    model.SpinningReserveCapability = Param(model.ThermalGenerators, within=NonNegativeReals, default=get_spin_bid)
    model.SpinningReservePrice = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0)
    
    # spinning reserve requirements
    model.ZonalSpinningReserveRequirement = Param(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, default=0.0)
    model.SystemSpinningReserveRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0)

    model.SpinningReservePenalty = Param(within=NonNegativeReals, mutable=True)

    def set_spinning_reserve_penalty(m):
        if m.SpinningReservePenalty(exception=False) is None:
            m.SpinningReservePenalty = value(m.RegulationPenalty+m.ReserveShortfallPenalty)/2.
    model.SetSpinningReservePenalty = BuildAction(rule=set_spinning_reserve_penalty)

    def zonal_spin_bounds(m,rz,t):
        return (0, m.ZonalSpinningReserveRequirement[rz,t])
    model.ZonalSpinningReserveShortfall = Var(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_spin_bounds)
    def system_spin_bounds(m,t):
        return (0, m.SystemSpinningReserveRequirement[t])
    model.SystemSpinningReserveShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_spin_bounds)

    # spinning reserve
    def spin_bounds(m,g,t):
        return (0,m.SpinningReserveCapability[g])
    model.SpinningReserveDispatched = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=spin_bounds)

    def spinning_reserve_available(m, g, t):
        spin_limit = min(value(m.SpinningReserveCapability[g]), value(m.NominalRampUpLimit[g]/60.*m.SpinningReserveMinutes))
        return m.SpinningReserveDispatched[g, t] <= spin_limit*m.UnitOn[g,t]
    model.SpinningReserveAvailableConstr = Constraint(model.ThermalGenerators, model.TimePeriods, rule=spinning_reserve_available)


    def zonal_spinning_reserve_provided(m, rz, t):
        return sum(m.SpinningReserveDispatched[g, t] for g in m.ThermalGeneratorsInReserveZone[rz]) + m.ZonalSpinningReserveShortfall[rz,t]
    model.ZonalSpinningReserveProvided = Expression(model.ReserveZones, model.TimePeriods, rule=zonal_spinning_reserve_provided)

    def enforce_zonal_spinning_reserve_requirement(m, rz, t):
        return  m.ZonalSpinningReserveProvided[rz,t] + m.ZonalRegulationUpProvided[rz,t] \
                >= m.ZonalSpinningReserveRequirement[rz, t] + m.ZonalRegulationUpRequirement[rz, t]
    model.EnforceZonalSpinningReserveRequirement = Constraint(model.ReserveZones, model.TimePeriods, rule=enforce_zonal_spinning_reserve_requirement)

    def system_spinning_reserve_provided(m,t):
        return sum(m.SpinningReserveDispatched[g,t] for g in m.ThermalGenerators) \
                + sum(m.ZonalSpinningReserveShortfall[rz,t] for rz in m.ReserveZones) \
                + m.SystemSpinningReserveShortfall[t]
    model.SystemSpinningReserveProvided = Expression(model.TimePeriods, rule=system_spinning_reserve_provided)

    def enforce_system_spinning_reserve_requirement(m, t):
        return m.SystemSpinningReserveProvided[t] + m.SystemRegulationUpProvided[t] \
                >= m.SystemSpinningReserveRequirement[t] + m.SystemRegulationUpRequirement[t]
    model.EnforceSystemSpinningReserveRequirement = Constraint(model.TimePeriods, rule=enforce_system_spinning_reserve_requirement)

    def compute_spinning_reserve_cost(m, g, t):
        return m.SpinningReserveDispatched[g, t] * m.SpinningReservePrice[g] * m.TimePeriodLengthHours
    model.SpinningReserveCostGeneration = Expression(model.ThermalGenerators, model.TimePeriods, rule=compute_spinning_reserve_cost)

    def spinning_reserve_cost_slacks(m,t):
        return m.TimePeriodLengthHours*m.SpinningReservePenalty*(
                      m.SystemSpinningReserveShortfall[t] \
                    + sum(m.ZonalSpinningReserveShortfall[rz,t] for rz in m.ReserveZones)
                    )
    model.SpinningReserveCostPenalty = Expression(model.TimePeriods, rule=spinning_reserve_cost_slacks)

    ## end spinning reserves

    ## begin non-spinning reserves
    # non-spinning reserve
    
    # Non-spinning reserves are assumed to be fast -- Operating reserves are slow (30 minutes)

    def validate_nonspin_bid(m,v,g):
        return v <= value(m.MaximumPowerOutput[g])
    model.NonSpinningReserveCapability = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0, validate=validate_nonspin_bid)
    model.NonSpinningReservePrice = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0)
    
    model.ZonalNonSpinningReserveRequirement = Param(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, default=0.0)
    model.SystemNonSpinningReserveRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0)

    model.NonSpinningReservePenalty = Param(within=NonNegativeReals, mutable=True)

    def set_non_spinning_reserve_penalty(m):
        if m.NonSpinningReservePenalty(exception=False) is None:
            m.NonSpinningReservePenalty = value(m.SpinningReservePenalty+m.ReserveShortfallPenalty)/2.
    model.SetNonSpinningReservePenalty = BuildAction(rule=set_non_spinning_reserve_penalty)

    def zonal_fast_bounds(m,rz,t):
        return (0, m.ZonalNonSpinningReserveRequirement[rz,t])
    model.ZonalNonSpinningReserveShortfall = Var(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_fast_bounds)
    def system_fast_bounds(m,t):
        return (0, m.SystemNonSpinningReserveRequirement[t])
    model.SystemNonSpinningReserveShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_fast_bounds)

    def nspin_bounds(m,g,t):
        return (0,m.NonSpinningReserveCapability[g])
    model.NonSpinningReserveDispatched = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=nspin_bounds)

    # non-spinning reserve
    def calculate_non_spinning_reserve_limit_rule(m, g, t):
        return m.NonSpinningReserveDispatched[g, t] <= m.NonSpinningReserveCapability[g] * (1 - m.UnitOn[g, t])
    model.CalculateNonSpinningReserveLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=calculate_non_spinning_reserve_limit_rule)

    def nspin_zonal_reserves_provided(m,rz,t):
        return sum(m.NonSpinningReserveDispatched[g,t] for g in m.ThermalGeneratorsInReserveZone[rz]) + m.ZonalNonSpinningReserveShortfall[rz,t]
    model.NonSpinningZonalReservesProvided = Expression(model.ReserveZones, model.TimePeriods, rule=nspin_zonal_reserves_provided)

    def enforce_zonal_non_spinning_reserve_rule(m, rz, t):
        return m.NonSpinningZonalReservesProvided[rz,t] + m.ZonalSpinningReserveProvided[rz,t] + m.ZonalRegulationUpProvided[rz,t] \
               >= m.ZonalNonSpinningReserveRequirement[rz, t] + m.ZonalSpinningReserveRequirement[rz,t] + m.ZonalRegulationUpRequirement[rz,t]
    model.EnforceNonSpinningZonalReserveRequirement = Constraint(model.ReserveZones, model.TimePeriods, rule=enforce_zonal_non_spinning_reserve_rule)

    def nspin_reserves_provided(m,t):
        return sum(m.NonSpinningReserveDispatched[g,t] for g in m.ThermalGenerators) \
                + sum(m.ZonalNonSpinningReserveShortfall[rz,t] for rz in m.ReserveZones) \
                + m.SystemNonSpinningReserveShortfall[t]
    model.SystemNonSpinningReserveProvided = Expression(model.TimePeriods, rule=nspin_reserves_provided)

    def enforce_system_non_spinning_reserve_requirement(m, t):
        return m.SystemNonSpinningReserveProvided[t] + m.SystemSpinningReserveProvided[t] + m.SystemRegulationUpProvided[t] \
                >= m.SystemNonSpinningReserveRequirement[t] + m.SystemSpinningReserveRequirement[t] + m.SystemRegulationUpRequirement[t]
    model.EnforceSystemNonSpinningReserveRequirement = Constraint(model.TimePeriods, rule=enforce_system_non_spinning_reserve_requirement)

    def calculate_non_spinning_reserve_cost(m, g, t):
        return m.NonSpinningReserveDispatched[g, t] * m.NonSpinningReservePrice[g] * m.TimePeriodLengthHours
    model.NonSpinningReserveCostGeneration = Expression(model.ThermalGenerators, model.TimePeriods, rule=calculate_non_spinning_reserve_cost)

    def non_spinning_reserve_cost_penalty(m,t):
        return m.TimePeriodLengthHours*m.NonSpinningReservePenalty*(
                        m.SystemNonSpinningReserveShortfall[t] \
                      + sum(m.ZonalNonSpinningReserveShortfall[rz,t] for rz in m.ReserveZones)
                      )
    model.NonSpinningReserveCostPenalty = Expression(model.TimePeriods, rule=non_spinning_reserve_cost_penalty)

    ## end non-spinning reserve

    ## begin operating reserve

    # Thirty-minute operating reserve, for generators which can start in 30 minutes
    model.OperatingReserveCapability = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0, validate=validate_nonspin_bid)
    model.OperatingReservePrice = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0)

    # Operating reserve requirement
    model.ZonalOperatingReserveRequirement = Param(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, default=0.0)
    model.SystemOperatingReserveRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0)

    model.OperatingReservePenalty = Param(within=NonNegativeReals, mutable=True)

    def set_operating_reserve_penalty(m):
        if m.OperatingReservePenalty(exception=False) is None:
            m.OperatingReservePenalty = value(m.NonSpinningReservePenalty+m.ReserveShortfallPenalty)/2.
    model.SetOperatingReservePenalty = BuildAction(rule=set_operating_reserve_penalty)

    def zonal_op_bounds(m,rz,t):
        return (0, m.ZonalOperatingReserveRequirement[rz,t])
    model.ZonalOperatingReserveShortfall = Var(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_op_bounds)
    def system_op_bounds(m,t):
        return (0, m.SystemOperatingReserveRequirement[t])
    model.SystemOperatingReserveShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_op_bounds)
    
    def op_bounds(m,g,t):
        return (0,m.OperatingReserveCapability[g])
    model.OperatingReserveDispatched = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=op_bounds)

    # thirty-minute operating reserve, for units which are off
    def calculate_operating_reserve_limit_rule(m, g, t):
        return m.OperatingReserveDispatched[g, t] + m.NonSpinningReserveDispatched[g, t] <= m.OperatingReserveCapability[g] * (1 - m.UnitOn[g, t])
    model.CalculateOperatingReserveLimits = Constraint(model.ThermalGenerators, model.TimePeriods, rule=calculate_operating_reserve_limit_rule)

    def operational_zonal_reserves_provided(m,rz,t):
        return sum(m.OperatingReserveDispatched[g,t] for g in m.ThermalGeneratorsInReserveZone[rz]) + m.ZonalOperatingReserveShortfall[rz,t]
    model.OperatingZonalReservesProvided = Expression(model.ReserveZones, model.TimePeriods, rule=operational_zonal_reserves_provided)

    def enforce_zonal_operating_reserve_requirement_rule(m, rz, t):
        return m.OperatingZonalReservesProvided[rz,t] + m.NonSpinningZonalReservesProvided[rz,t] + m.ZonalSpinningReserveProvided[rz,t] + m.ZonalRegulationUpRequirement[rz,t] \
                >= m.ZonalOperatingReserveRequirement[rz,t] + m.ZonalNonSpinningReserveRequirement[rz,t] + m.ZonalSpinningReserveRequirement[rz,t] + m.ZonalRegulationUpRequirement[rz,t]
    model.EnforceZonalOperatingReserveRequirement = Constraint(model.ReserveZones, model.TimePeriods, rule=enforce_zonal_operating_reserve_requirement_rule)

    def operational_reserves_provided(m,t):
        return sum(m.OperatingReserveDispatched[g,t] for g in m.ThermalGenerators) \
                + sum(m.ZonalOperatingReserveShortfall[rz,t] for rz in m.ReserveZones) \
                + m.SystemOperatingReserveShortfall[t]
    model.SystemOperatingReserveProvided = Expression(model.TimePeriods, rule=nspin_reserves_provided)

    def enforce_system_operating_reserve_requirement(m, t):
        return m.SystemOperatingReserveProvided[t] + m.SystemNonSpinningReserveProvided[t] + m.SystemSpinningReserveProvided[t] + m.SystemRegulationUpRequirement[t] >= \
                m.SystemOperatingReserveRequirement[t] + m.SystemNonSpinningReserveRequirement[t] + m.SystemSpinningReserveRequirement[t] + m.SystemRegulationUpRequirement[t]

    model.EnforceSystemOperatingReserveRequirement = Constraint(model.TimePeriods, rule=enforce_system_operating_reserve_requirement)

    def calculate_operating_reserve_cost_rule(m, g, t):
        return m.OperatingReserveDispatched[g, t] * m.OperatingReservePrice[g] * m.TimePeriodLengthHours
    model.OperatingReserveCostGeneration = Expression(model.ThermalGenerators, model.TimePeriods, rule=calculate_operating_reserve_cost_rule)

    def operating_reserve_cost_penalty(m,t):
        return m.TimePeriodLengthHours*m.OperatingReservePenalty*(
                        m.SystemOperatingReserveShortfall[t] \
                      + sum(m.ZonalOperatingReserveShortfall[rz,t] for rz in m.ReserveZones)
                      )
    model.OperatingReserveCostPenalty = Expression(model.TimePeriods, rule=operating_reserve_cost_penalty)

    ## end operating reserve


    ## begin flexible_ramp
    model.FlexRampMinutes = Param(within=PositiveReals, default=20.)

    model.ZonalFlexUpRequirement = Param(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, default=0.0)
    model.ZonalFlexDnRequirement = Param(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, default=0.0)

    model.SystemFlexUpRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0)
    model.SystemFlexDnRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0)

    def zonal_flex_up_bounds(m, rz, t):
        return (0, m.ZonalFlexUpRequirement[rz,t])
    model.ZonalFlexUpShortfall = Var(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_flex_up_bounds)
    def zonal_flex_dn_bounds(m, rz, t):
        return (0, m.ZonalFlexDnRequirement[rz,t])
    model.ZonalFlexDnShortfall = Var(model.ReserveZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_flex_dn_bounds)

    def system_flex_up_bounds(m, t):
        return (0, m.SystemFlexUpRequirement[t])
    model.SystemFlexUpShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_flex_up_bounds)
    def system_flex_dn_bounds(m, t):
        return (0, m.SystemFlexDnRequirement[t])
    model.SystemFlexDnShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_flex_dn_bounds)

    model.FlexRampPenalty = Param(within=NonNegativeReals, mutable=True)

    def set_flex_ramp_penalty(m):
        if m.FlexRampPenalty(exception=False) is None:
            m.FlexRampPenalty = value(m.NonSpinningReservePenalty+m.SpinningReservePenalty)/2.
    model.SetFlexRampPenalty = BuildAction(rule=set_flex_ramp_penalty)

    model.FlexUpProvided = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)
    model.FlexDnProvided = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)

    def flex_up_limit_rule(m, g, t):
        return m.FlexUpProvided[g,t] <= m.FlexRampMinutes*(m.NominalRampUpLimit[g]/60.)*m.UnitOn[g,t]
    model.FlexUpLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=flex_up_limit_rule)

    def flex_down_limit_rule(m, g, t):
        return m.FlexDnProvided[g,t] <= m.FlexRampMinutes*(m.NominalRampDownLimit[g]/60.)*m.UnitOn[g,t]
    model.FlexDnLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=flex_down_limit_rule)

    def zonal_flex_up_requirement_rule(m, rz, t):
        return sum(m.FlexUpProvided[g,t] for g in m.ThermalGeneratorsInReserveZone[rz]) + m.ZonalFlexUpShortfall[rz,t] >= m.ZonalFlexUpRequirement[t]
    model.ZonalFlexUpRequirementConstr = Constraint(model.ReserveZones, model.TimePeriods, rule=zonal_flex_up_requirement_rule)

    def zonal_flex_dn_requirement_rule(m, rz, t):
        return sum(m.FlexDnProvided[g,t] for g in m.ThermalGeneratorsInReserveZone[rz]) + m.ZonalFlexDnShortfall[rz,t] >= m.ZonalFlexDnRequirement[t]
    model.ZonalFlexDnRequirementConstr = Constraint(model.ReserveZones, model.TimePeriods, rule=zonal_flex_dn_requirement_rule)

    def system_flex_up_requirement_rule(m, t):
        return sum(m.FlexUpProvided[g,t] for g in m.ThermalGenerators) \
                 + sum(m.ZonalFlexUpShortfall[rz,t] for rz in m.ReserveZones) \
                 + m.SystemFlexUpShortfall[t] \
                 >= m.SystemFlexUpRequirement[t]
    model.SystemFlexUpRequirementConstr = Constraint(model.ReserveZones, model.TimePeriods, rule=system_flex_up_requirement_rule)

    def system_flex_dn_requirement_rule(m, t):
        return sum(m.FlexDnProvided[g,t] for g in m.ThermalGenerators) \
                 + sum(m.ZonalFlexDnShortfall[rz,t] for rz in m.ReserveZones) \
                 + m.SystemFlexDnShortfall[t] \
                 >= m.SystemFlexDnRequirement[t]
    model.SystemFlexDnRequirementConstr = Constraint(model.ReserveZones, model.TimePeriods, rule=system_flex_dn_requirement_rule)
    
    def flex_ramp_penalty_cost(m, t):
        return m.TimePeriodLengthHours*m.FlexRampPenalty*(
                        m.SystemFlexUpShortfall[t] + m.SystemFlexDnShortfall[t] \
                      + sum(m.ZonalFlexUpShortfall[rz,t]+m.ZonalFlexDnShortfall[rz,t] for rz in m.ReserveZones)
                      )
    model.FlexibleRampingCostPenalty = Expression(model.TimePeriods, rule = flex_ramp_penalty_cost)


    ## end flexible ramp

    ## Ancillary service capacity limits (enhance for ramping, start-up/shutdown)

    def ancillary_service_capacity_limit_upper(m, g, t):
        return m.PowerGenerated[g,t] + m.FlexUpProvided[g,t] + m.RegulationReserveUp[g,t] + m.SpinningReserveDispatched[g,t] <= \
                    m.MaximumPowerOutput[g]*m.UnitOn[g,t] - (m.MaximumPowerOutput[g] - m.RegulationHighLimit[g])*m.RegulationOn[g,t]
    model.AncillaryServiceCapacityLimitUpper = Constraint(model.ThermalGenerators, model.TimePeriods, rule=ancillary_service_capacity_limit_upper)

    def ancillary_service_capacity_limit_lower(m, g, t):
        return m.PowerGeneratedAboveMinimum[g,t] - m.FlexDnProvided[g,t] - m.RegulationReserveDn[g,t] >= (m.RegulationLowLimit[g] - m.MinimumPowerOutput[g])*m.RegulationOn[g,t]
    model.AncillaryServiceCapacityLimitLower = Constraint(model.ThermalGenerators, model.TimePeriods, rule=ancillary_service_capacity_limit_lower)

    ## NOTE: ScaledNominalRampUpLimit/ScaledNominalRampDownLimit and ScaledStartupRampLimit/ScaledShutdownRampLimit
    ##       are not appropriate in the ramp sharing constraints that follow.
    ##       In particular, we need to possibly allow these to be larger than (MaximumPowerGenerated -
    ##       MinimumPowerGenerated), which these ramp limts do not allow for tightness and less error checking
    ##       in the base UC/ED constrants

    def as_ramp_up(m,g):
        return m.NominalRampUpLimit[g]*m.TimePeriodLengthHours
    model.AS_ScaledNominalRampUpLimit = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=as_ramp_up, mutable=True)

    def as_ramp_down(m,g):
        return m.NominalRampDownLimit[g]*m.TimePeriodLengthHours
    model.AS_ScaledNominalRampDownLimit = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=as_ramp_down, mutable=True)

    def as_startup_ramp(m,g):
        return (m.StartupRampLimit[g] - m.MinimumPowerOutput[g])*m.TimePeriodLengthHours
    model.AS_ScaledStartupRampLessMin = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=as_startup_ramp, mutable=True)

    def as_shutdown_ramp(m,g):
        return (m.ShutdownRampLimit[g] - m.MinimumPowerOutput[g])*m.TimePeriodLengthHours
    model.AS_ScaledShutdownRampLessMin = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=as_shutdown_ramp, mutable=True)

    def ancillary_service_ramp_up_limit(m,g,t):
        if t == m.InitialTime:
            return m.MaximumPowerAvailableAboveMinimum[g, t] - (m.PowerGeneratedT0[g]-m.MinimumPowerOutput[g]*m.UnitOnT0[g]) \
                     + (m.TimePeriodLengthMinutes/m.RegulationMinutes)*m.RegulationReserveUp[g,t] \
                     + (m.TimePeriodLengthMinutes/m.SpinningReserveMinutes)*m.SpinningReserveDispatched[g,t] \
                     + (m.TimePeriodLengthMinutes/m.FlexRampMinutes)*m.FlexUpProvided[g,t] \
                  <= \
                    m.AS_ScaledNominalRampUpLimit[g]*m.UnitOn[g,t] + \
    		    (m.AS_ScaledStartupRampLessMin[g] - m.AS_ScaledNominalRampUpLimit[g])*m.UnitStart[g,t] 
        else: ## average the regulation and spin over the two time periods, which is what is done in CAISO
            return m.MaximumPowerAvailableAboveMinimum[g, t] - m.PowerGeneratedAboveMinimum[g, t-1] \
                     + (m.TimePeriodLengthMinutes/m.RegulationMinutes)*(m.RegulationReserveUp[g,t]+m.RegulationReserveUp[g,t-1])/2. \
                     + (m.TimePeriodLengthMinutes/m.SpinningReserveMinutes)*(m.SpinningReserveDispatched[g,t]+m.SpinningReserveDispatched[g,t-1])/2. \
                     + (m.TimePeriodLengthMinutes/m.FlexRampMinutes)*m.FlexUpProvided[g,t] \
                  <= \
                    m.AS_ScaledNominalRampUpLimit[g]*m.UnitOn[g,t] + \
    		    (m.AS_ScaledStartupRampLessMin[g] - m.AS_ScaledNominalRampUpLimit[g])*m.UnitStart[g,t] 
    model.AncillaryServiceRampUpLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=ancillary_service_ramp_up_limit)

    ## NOTE: for the regulation and flexible down services, these subtract from power generated at t, so they get added here
    def ancillary_service_ramp_dn_limit(m,g,t):
        if t == m.InitialTime:
            if not m.enforce_t1_ramp_rates:
                return Constraint.Skip
            else:
                return (m.PowerGeneratedT0[g] - m.MinimumPowerOutput[g]*m.UnitOnT0[g]) - m.PowerGeneratedAboveMinimum[g, t] \
                          + (m.TimePeriodLengthMinutes/m.RegulationMinutes)*m.RegulationReserveDn[g,t] \
                          + (m.TimePeriodLengthMinutes/m.FlexRampMinutes)*m.FlexDnProvided[g,t] \
                       <= \
                         m.AS_ScaledNominalRampDownLimit[g]*m.UnitOnT0[g] + \
                         (m.AS_ScaledShutdownRampLessMin[g] - m.AS_ScaledNominalRampDownLimit[g])*m.UnitStop[g,t]
        else:
            return m.PowerGeneratedAboveMinimum[g, t-1] - m.PowerGeneratedAboveMinimum[g, t] \
                     + (m.TimePeriodLengthMinutes/m.RegulationMinutes)*(m.RegulationReserveDn[g,t]+m.RegulationReserveDn[g,t-1])/2. \
                     + (m.TimePeriodLengthMinutes/m.FlexRampMinutes)*m.FlexDnProvided[g,t] \
                  <= \
                    m.AS_ScaledNominalRampDownLimit[g]*m.UnitOn[g,t-1] + \
                    (m.AS_ScaledShutdownRampLessMin[g] - m.AS_ScaledNominalRampDownLimit[g])*m.UnitStop[g,t]
    model.AncillaryServiceRampDnLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=ancillary_service_ramp_dn_limit)
