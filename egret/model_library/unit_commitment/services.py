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

from .uc_utils import add_model_attr, uc_time_helper 
from .status_vars import _is_relaxed

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
    if _is_relaxed(model):
        model.InputStorage = Var(model.Storage, model.TimePeriods, within=UnitInterval)
        model.OutputStorage = Var(model.Storage, model.TimePeriods, within=UnitInterval)
    else:
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
        if value(m.ScaledNominalRampUpLimitStorageOutput[s]) >= \
                value(m.MaximumPowerOutputStorage[s]-m.MinimumPowerOutputStorage[s]):
            return Constraint.Skip
        if t == m.InitialTime:
            return m.PowerOutputStorage[s, t] <= m.StoragePowerOutputOnT0[s] + m.ScaledNominalRampUpLimitStorageOutput[s]
        else:
            return m.PowerOutputStorage[s, t] <= m.PowerOutputStorage[s, t-1] + m.ScaledNominalRampUpLimitStorageOutput[s]

    model.EnforceStorageOutputRampUpRates = Constraint(model.Storage, model.TimePeriods, rule=enforce_ramp_up_rates_power_output_storage_rule)

    def enforce_ramp_down_rates_power_output_storage_rule(m, s, t):
        if value(m.ScaledNominalRampDownLimitStorageOutput[s]) >= \
                value(m.MaximumPowerOutputStorage[s]-m.MinimumPowerOutputStorage[s]):
            return Constraint.Skip
        if t == m.InitialTime:
            return m.PowerOutputStorage[s, t] >= m.StoragePowerOutputOnT0[s] - m.ScaledNominalRampDownLimitStorageOutput[s]
        else:
            return m.PowerOutputStorage[s, t] >= m.PowerOutputStorage[s, t-1] - m.ScaledNominalRampDownLimitStorageOutput[s]

    model.EnforceStorageOutputRampDownRates = Constraint(model.Storage, model.TimePeriods, rule=enforce_ramp_down_rates_power_output_storage_rule)

    def enforce_ramp_up_rates_power_input_storage_rule(m, s, t):
        if value(m.ScaledNominalRampUpLimitStorageInput[s]) >= \
                value(m.MaximumPowerInputStorage[s]-m.MinimumPowerInputStorage[s]):
            return Constraint.Skip
        if t == m.InitialTime:
            return m.PowerInputStorage[s, t] <= m.StoragePowerInputOnT0[s] + m.ScaledNominalRampUpLimitStorageInput[s]
        else:
            return m.PowerInputStorage[s, t] <= m.PowerInputStorage[s, t-1] + m.ScaledNominalRampUpLimitStorageInput[s]

    model.EnforceStorageInputRampUpRates = Constraint(model.Storage, model.TimePeriods, rule=enforce_ramp_up_rates_power_input_storage_rule)

    def enforce_ramp_down_rates_power_input_storage_rule(m, s, t):
        if value(m.ScaledNominalRampDownLimitStorageInput[s]) >= \
                value(m.MaximumPowerInputStorage[s]-m.MinimumPowerInputStorage[s]):
            return Constraint.Skip
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
        return m.SocStorage[s, value(m.NumTimePeriods)] >= m.EndPointSocStorage[s]
    model.EnforceEndPointSocStorage = Constraint(model.Storage, rule=storage_end_point_soc_rule)

    def storage_cost_rule(m, s, t):
        return m.ChargeCost[s]*m.PowerInputStorage[s,t]*m.TimePeriodLengthHours + \
                m.DischargeCost[s]*m.PowerOutputStorage[s,t]*m.TimePeriodLengthHours
    model.StorageCost = Expression(model.Storage, model.TimePeriods, rule=storage_cost_rule)

    return
## end storage_services


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
                                                  'status_vars': None,
                                                  'power_vars': None,
                                                 })
def ancillary_services(model):
    '''
    Defines ancillary services: regulation, spinning reserve, nonspinning reserve, operational reserve, flexible ramp
    ## NOTE: As in most markets, the value of ancillary services from high to low is regulation, spinning reserve, nonspinning reserve, and supplemental reserve.
    ##       We allow for a higher-quality ancillary service to be subtituted for a lower-quality one
    ##       Flexible ramp is treated differently, again as it is in most markets. There is no bid for flexible ramp, and it is priced at opportunity cost
    '''
    md = model.model_data

    system = md.data['system']
    elements = md.data['elements']

    TimeMapper = uc_time_helper(model.TimePeriods)


    ## list of possible ancillary services coming
    ## from model_data
    ancillary_service_list = [ 'spinning_reserve_requirement',
                               'non_spinning_reserve_requirement',
                               'regulation_up_requirement',
                               'regulation_down_requirement',
                               'supplemental_reserve_requirement',
                               'flexible_ramp_up_requirement',
                               'flexible_ramp_down_requirement',
                             ]

    if 'zone' not in elements:
        elements['zone'] = dict()
    if 'area' not in elements:
        elements['area'] = dict()
    
    ## check and see if each one of these services appears anywhere in model_data
    def _check_for_requirement( requirement ):
        if requirement in system:
            return True
        for zone in elements['zone'].values():
            if requirement in zone:
                return True
        for area in elements['area'].values():
            if requirement in area:
                return True
        return False

    ## flags for if ancillary services appear
    add_spinning_reserve =  _check_for_requirement('spinning_reserve_requirement')
    add_non_spinning_reserve = _check_for_requirement('non_spinning_reserve_requirement')
    add_regulation_reserve = (_check_for_requirement('regulation_up_requirement') or 
                            _check_for_requirement('regulation_down_requirement'))
    add_supplemental_reserve = _check_for_requirement('supplemental_reserve_requirement')
    add_flexi_ramp_reserve = (_check_for_requirement('flexible_ramp_up_requirement') or
                            _check_for_requirement('flexible_ramp_down_requirement'))

    
    ## check here and break if there's nothing to do
    no_reserves = not (add_spinning_reserve or add_non_spinning_reserve or add_regulation_reserve or add_supplemental_reserve or add_flexi_ramp_reserve)

    ## add a flag for which brach we took here
    if no_reserves:
        model.nonbasic_reserves = False
        model.regulation_service = None
        model.spinning_reserve = None
        model.non_spinning_reserve = None
        model.supplemental_reserve = None
        model.flexible_ramping = None
        return

    model.nonbasic_reserves = True

    ## check this here to avoid exceptions when the model has no ancillary services
    if model.status_vars not in ['garver_3bin_vars','garver_2bin_vars', 'garver_3bin_relaxed_stop_vars', 'ALS_state_transition_vars']:
        raise Exception('Exception adding ancillary_services! ancillary_services requires one of: garver_3bin_vars, garver_2bin_vars, garver_3bin_relaxed_stop_vars, ALS_state_transition_vars, to be used for the status_vars.')

    ## set some penalties by default based on the other model penalties
    default_reg_pen = value(model.LoadMismatchPenalty+model.ReserveShortfallPenalty)/2.
    ## set these penalties in relation to each other, from higher quality service to lower
    model.RegulationPenalty = Param(within=NonNegativeReals,
                                    initialize=system.get('regulation_penalty_price', default_reg_pen))

    default_spin_pen = value(model.RegulationPenalty+model.ReserveShortfallPenalty)/2.
    model.SpinningReservePenalty = Param(within=NonNegativeReals, 
                                         initialize=system.get('spinning_reserve_penalty_price', default_spin_pen))

    default_nspin_pen = value(model.SpinningReservePenalty+model.ReserveShortfallPenalty)/2.
    model.NonSpinningReservePenalty = Param(within=NonNegativeReals,
                                            initialize=system.get('non_spinning_reserve_penalty_price', default_nspin_pen))

    default_supp_pen = value(model.NonSpinningReservePenalty+model.ReserveShortfallPenalty)/2.
    model.SupplementalReservePenalty = Param(within=NonNegativeReals,
                                             initialize=system.get('supplemental_reserve_penalty_price', default_supp_pen))

    default_flex_pen = value(model.NonSpinningReservePenalty+model.SpinningReservePenalty)/2.
    model.FlexRampPenalty = Param(within=NonNegativeReals,
                                  initialize=system.get('flexible_ramp_penalty_price', default_flex_pen))

    thermal_gen_attrs = md.attributes(element_type='generator', generator_type='thermal')
    
    def zone_initializer_builder(reserve_checker):
        def init_reserve_zone(m):
            for an, area in elements['area'].items():
                if reserve_checker(area):
                    yield 'area_'+an
            for zn, zone in elements['zone'].items():
                if reserve_checker(zone):
                    yield 'zone_'+zn
        return init_reserve_zone

    zone_attrs = md.attributes(element_type='zone')
    area_attrs = md.attributes(element_type='area')

    def zone_requirement_getter(reserve_product):
        if reserve_product in zone_attrs:
            zone_r_time = TimeMapper(zone_attrs[reserve_product])
        if reserve_product in area_attrs:
            area_r_time = TimeMapper(area_attrs[reserve_product])
        def get_attribute(m, az, t):
            az_n = str(az)
            if az_n[:5] == 'zone_':
                z_n = az_n[5:]
                if (z_n,t) in zone_r_time:
                    return zone_r_time[z_n,t]
                else:
                    return 0.0
            elif az_n[:5] == 'area_':
                a_n = az_n[5:]
                if (a_n,t) in area_r_time:
                    return area_r_time[a_n,t]
                else:
                    return 0.0
            else:
                raise Exception('Unexpected case in instance of zone_requirement_getter')
        return get_attribute
    
    def gens_in_reserve_zone_getter(gen_attrs_subset=None):
        if gen_attrs_subset is None:
            gen_attrs = thermal_gen_attrs
        else:
            gen_attrs = gen_attrs_subset
        def get_gens_in_reserve_zone(m, az):
            az_n = str(az)
            if az_n[:5] == 'zone_':
                z_n = az_n[5:]
                for g in gen_attrs['names']:
                    if g in gen_attrs['area'] and gen_attrs['zone'][g] == z_n:
                        yield g
            elif az_n[:5] == 'area_':
                a_n = az_n[5:]
                for g in gen_attrs['names']:
                    if g in gen_attrs['area'] and gen_attrs['area'][g] == a_n:
                        yield g
            else:
                raise Exception('Unexpected case in instance of gens_in_reserve_zone_getter')
        return get_gens_in_reserve_zone

    ## these need to be added by high-quality to low-quality,
    ## except flexiramp which is it's own thing
    if add_regulation_reserve:
        regulation_services(model, zone_initializer_builder, zone_requirement_getter, gens_in_reserve_zone_getter)
    else:
        model.regulation_service = None

    if add_spinning_reserve:
        spinning_reserves(model, zone_initializer_builder, zone_requirement_getter, gens_in_reserve_zone_getter, thermal_gen_attrs)
    else:
        model.spinning_reserve = None

    if add_non_spinning_reserve:
        non_spinning_reserves(model, zone_initializer_builder, zone_requirement_getter, gens_in_reserve_zone_getter)
    else:
        model.non_spinning_reserve = None

    if add_supplemental_reserve:
        supplemental_reserves(model, zone_initializer_builder, zone_requirement_getter, gens_in_reserve_zone_getter, thermal_gen_attrs)
    else:
        model.supplemental_reserve = None

    if add_flexi_ramp_reserve:
        flexible_ramping_reserves(model, zone_initializer_builder, zone_requirement_getter, gens_in_reserve_zone_getter, thermal_gen_attrs)
    else:
        model.flexible_ramping = None

    ## Ancillary service capacity limits (enhance for ramping, start-up/shutdown)

    def ancillary_service_capacity_limit_upper(m, g, t):
        reg = (bool(m.regulation_service) and (g in m.AGC_Generators))
        return m.MaximumPowerAvailable[g,t] \
                    + (m.FlexUpProvided[g,t] if add_flexi_ramp_reserve else 0.) \
                    + (m.RegulationReserveUp[g,t] if reg else 0.) \
                    + (m.SpinningReserveDispatched[g,t] if add_spinning_reserve else 0.) \
                    + (m.SupplementalSpinReserveDispatched[g,t] if add_supplemental_reserve else 0.) \
                <= m.MaximumPowerOutput[g,t]*m.UnitOn[g,t] \
                    - ((m.MaximumPowerOutput[g,t] - m.RegulationHighLimit[g,t])*m.RegulationOn[g,t] if reg else 0.)
    model.AncillaryServiceCapacityLimitUpper = Constraint(model.ThermalGenerators, model.TimePeriods, rule=ancillary_service_capacity_limit_upper)

    def ancillary_service_capacity_limit_lower(m, g, t):
        if not (bool(m.flexible_ramping) or bool(m.regulation_service)):
            return Constraint.Feasible
        reg = (bool(m.regulation_service) and (g in m.AGC_Generators))
        return m.PowerGeneratedAboveMinimum[g,t] \
                    - (m.FlexDnProvided[g,t] if add_flexi_ramp_reserve else 0.) \
                    - (m.RegulationReserveDn[g,t] if reg else 0.) \
                >= \
                    ((m.RegulationLowLimit[g,t] - m.MinimumPowerOutput[g,t])*m.RegulationOn[g,t] if reg else 0.)
    model.AncillaryServiceCapacityLimitLower = Constraint(model.ThermalGenerators, model.TimePeriods, rule=ancillary_service_capacity_limit_lower)

    ## NOTE: ScaledNominalRampUpLimit/ScaledNominalRampDownLimit and ScaledStartupRampLimit/ScaledShutdownRampLimit
    ##       are not appropriate in the ramp sharing constraints that follow.
    ##       In particular, we need to possibly allow these to be larger than (MaximumPowerGenerated -
    ##       MinimumPowerGenerated), which these ramp limts do not allow for tightness and less error checking
    ##       in the base UC/ED constrants

    def as_ramp_up(m,g):
        return m.NominalRampUpLimit[g]*m.TimePeriodLengthHours
    model.AS_ScaledNominalRampUpLimit = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=as_ramp_up)

    def as_ramp_down(m,g):
        return m.NominalRampDownLimit[g]*m.TimePeriodLengthHours
    model.AS_ScaledNominalRampDownLimit = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=as_ramp_down)

    ##TODO: FIXME: REVISIT AFTER RAMPING CONSTRAINTS
    def as_startup_ramp(m,g,t):
        return (m.StartupRampLimit[g,t] - m.MinimumPowerOutput[g,t])*m.TimePeriodLengthHours + m.MinimumPowerOutput[g,t]
    model.AS_ScaledStartupRamp = Param(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, initialize=as_startup_ramp)

    def as_shutdown_ramp(m,g,t):
        return (m.ShutdownRampLimit[g,t] - m.MinimumPowerOutput[g,t])*m.TimePeriodLengthHours + m.MinimumPowerOutput[g,t]
    model.AS_ScaledShutdownRamp = Param(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, initialize=as_shutdown_ramp)

    def as_shutdown_ramp_t0(m,g):
        return (m.ShutdownRampLimitT0[g] - m.MinimumPowerOutputT0[g])*m.TimePeriodLengthHours + m.MinimumPowerOutputT0[g]
    model.AS_ScaledShutdownRampT0 = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=as_shutdown_ramp_t0)

    ## These are formulated similarly to the damci-kurt ramp limits
    def ancillary_service_ramp_up_limit(m,g,t):
        reg = (bool(m.regulation_service) and (g in m.AGC_Generators))
        if t == m.InitialTime:
            return m.MaximumPowerAvailableAboveMinimum[g, t] - m.PowerGeneratedT0[g]\
                     + ((m.TimePeriodLengthMinutes/m.RegulationMinutes)*m.RegulationReserveUp[g,t] if reg else 0.) \
                     + ((m.TimePeriodLengthMinutes/m.SpinningReserveMinutes)*m.SpinningReserveDispatched[g,t] if add_spinning_reserve else 0.) \
                     + ((m.TimePeriodLengthMinutes/m.FlexRampMinutes)*m.FlexUpProvided[g,t] if add_flexi_ramp_reserve else 0.) \
                     + ((m.TimePeriodLengthMinutes/m.SupplementalReserveMinutes)*m.SupplementalSpinReserveDispatched[g,t] if add_supplemental_reserve else 0.) \
                  <= \
                    (m.AS_ScaledNominalRampUpLimit[g] + 0 - m.MinimumPowerOutput[g,t])*m.UnitOn[g,t] + \
    		    (m.AS_ScaledStartupRamp[g,t] - 0 - m.AS_ScaledNominalRampUpLimit[g])*m.UnitStart[g,t] 
        else: ## average the regulation and spin over the two time periods, which is what is done in CAISO
            return m.MaximumPowerAvailableAboveMinimum[g, t] - m.PowerGeneratedAboveMinimum[g, t-1] \
                     + ((m.TimePeriodLengthMinutes/m.RegulationMinutes)*(m.RegulationReserveUp[g,t]+m.RegulationReserveUp[g,t-1])/2. if reg else 0.) \
                     + ((m.TimePeriodLengthMinutes/m.SpinningReserveMinutes)*(m.SpinningReserveDispatched[g,t]+m.SpinningReserveDispatched[g,t-1])/2. if add_spinning_reserve else 0.) \
                     + ((m.TimePeriodLengthMinutes/m.FlexRampMinutes)*m.FlexUpProvided[g,t] if add_flexi_ramp_reserve else 0.) \
                     + ((m.TimePeriodLengthMinutes/m.SupplementalReserveMinutes)*(m.SupplementalSpinReserveDispatched[g,t]+m.SupplementalSpinReserveDispatched[g,t-1])/2. if add_supplemental_reserve else 0.) \
                  <= \
                    (m.AS_ScaledNominalRampUpLimit[g] + m.MinimumPowerOutput[g,t-1] - m.MinimumPowerOutput[g,t])*m.UnitOn[g,t] + \
    		    (m.AS_ScaledStartupRamp[g,t] - m.MinimumPowerOutput[g,t-1] - m.AS_ScaledNominalRampUpLimit[g])*m.UnitStart[g,t] 
    model.AncillaryServiceRampUpLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=ancillary_service_ramp_up_limit)

    ## NOTE: for the regulation and flexible down services, these subtract from power generated at t, so they get added here
    def ancillary_service_ramp_dn_limit(m,g,t):
        if not (bool(m.flexible_ramping) or bool(m.regulation_service)):
            return Constraint.Feasible
        reg = (bool(m.regulation_service) and (g in m.AGC_Generators))
        if t == m.InitialTime:
            if not m.enforce_t1_ramp_rates:
                return Constraint.Skip
            else:
                return m.PowerGeneratedT0[g] - m.PowerGeneratedAboveMinimum[g, t] \
                          + ((m.TimePeriodLengthMinutes/m.RegulationMinutes)*m.RegulationReserveDn[g,t] if reg else 0.) \
                          + ((m.TimePeriodLengthMinutes/m.FlexRampMinutes)*m.FlexDnProvided[g,t] if add_flexi_ramp_reserve else 0.) \
                       <= \
                         (m.AS_ScaledNominalRampDownLimit[g] + m.MinimumPowerOutput[g,t] - 0)*m.UnitOnT0[g] + \
                         (m.AS_ScaledShutdownRampT0[g] - m.MinimumPowerOutput[g,t] - m.AS_ScaledNominalRampDownLimit[g])*m.UnitStop[g,t]
        else:
            return m.PowerGeneratedAboveMinimum[g, t-1] - m.PowerGeneratedAboveMinimum[g, t] \
                     + ((m.TimePeriodLengthMinutes/m.RegulationMinutes)*(m.RegulationReserveDn[g,t]+m.RegulationReserveDn[g,t-1])/2. if reg else 0.) \
                     + ((m.TimePeriodLengthMinutes/m.FlexRampMinutes)*m.FlexDnProvided[g,t] if add_flexi_ramp_reserve else 0.) \
                  <= \
                    (m.AS_ScaledNominalRampDownLimit[g] + m.MinimumPowerOutput[g,t] - m.MinimumPowerOutput[g,t-1])*m.UnitOn[g,t-1] + \
                    (m.AS_ScaledShutdownRamp[g,t-1] - m.MinimumPowerOutput[g,t] - m.AS_ScaledNominalRampDownLimit[g])*m.UnitStop[g,t]
    model.AncillaryServiceRampDnLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=ancillary_service_ramp_dn_limit)


@add_model_attr('regulation_service', requires = {'data_loader': None,
                                                  'status_vars': None,
                                                 })
def regulation_services(model, zone_initializer_builder, zone_requirement_getter, gens_in_reserve_zone_getter):

    md = model.model_data

    system = md.data['system']

    TimeMapper = uc_time_helper(model.TimePeriods)

    def _check_reg(e_dict):
        return ( ('regulation_up_requirement' in e_dict) \
                  or ('regulation_down_requirement' in e_dict) )
    model.RegulationZones = Set(initialize=zone_initializer_builder(_check_reg))

    ## begin regulation

    #################################
    # Regulation ancillary services #
    #################################

    agc_gen_attrs = md.attributes(element_type='generator', generator_type='thermal', agc_capable=True)
    
    model.AGC_Generators = Set(within=model.ThermalGenerators, initialize=agc_gen_attrs['names'])
    
    model.AGC_GeneratorsInRegulationZone = Set(model.RegulationZones, initialize=gens_in_reserve_zone_getter(agc_gen_attrs))

    model.RegulationMinutes = Param(within=PositiveReals, default=5.)
    
    # When units are selected for regulation, their limits are bounded by the RegulationHighLimit and RegulationLowLimit
    # I'll refer to it as the "regulation band"
    def regulation_high_limit_validator(m, v, g, t):
        return v <= value(m.MaximumPowerOutput[g,t])
    model.RegulationHighLimit = Param(model.AGC_Generators, model.TimePeriods, within=NonNegativeReals, validate=regulation_high_limit_validator, initialize=TimeMapper(agc_gen_attrs['p_max_agc']))
    
    def regulation_low_limit_validator(m, v, g, t):
        return (v <= value(m.RegulationHighLimit[g,t]) and v >= value(m.MinimumPowerOutput[g,t]))
    model.RegulationLowLimit = Param(model.AGC_Generators, model.TimePeriods, within=NonNegativeReals, validate=regulation_low_limit_validator, initialize=TimeMapper(agc_gen_attrs['p_min_agc']))
    
    # Regulation capacity is calculated as the min of "regulation band" and RegulationMinutes*AutomaticResponseRate
    model.AutomaticResponseRate = Param(model.AGC_Generators, within=NonNegativeReals, initialize=agc_gen_attrs['ramp_agc'])
    
    def calculate_regulation_capability_rule(m, g, t):
        temp1 = value(m.RegulationMinutes * m.AutomaticResponseRate[g])
        temp2 = value(m.RegulationHighLimit[g,t] - m.RegulationLowLimit[g,t])/2.
        if temp1 > temp2:
            return temp2
        else:
            return temp1
    
    model.RegulationUpCapability = Param(model.AGC_Generators, model.TimePeriods, within=NonNegativeReals, initialize=calculate_regulation_capability_rule)
    model.RegulationDnCapability = Param(model.AGC_Generators, model.TimePeriods, within=NonNegativeReals, initialize=calculate_regulation_capability_rule)

    model.ZonalRegulationUpRequirement = Param(model.RegulationZones, model.TimePeriods, within=NonNegativeReals, 
                                                    initialize=zone_requirement_getter('regulation_up_requirement'))
    
    model.SystemRegulationUpRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0, initialize=TimeMapper(system.get('regulation_up_requirement', dict())))

    model.ZonalRegulationDnRequirement = Param(model.RegulationZones, model.TimePeriods, within=NonNegativeReals,
                                                    initialize=zone_requirement_getter('regulation_down_requirement'))

    model.SystemRegulationDnRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0, initialize=TimeMapper(system.get('regulation_down_requirement', dict())))

    def validate_fixed_reg(m,v,g,t):
        if (v is not None) and (value(m.FixedCommitment[g,t]) is not None):
            return v <= value(m.FixedCommitment[g,t])
        else:
            return True
    model.FixedRegulation = Param(model.AGC_Generators, model.TimePeriods, initialize=TimeMapper(agc_gen_attrs.get('fixed_regulation', dict())),
                                    default=None, within=model.FixedCommitmentTypes, validate=validate_fixed_reg)

    def zonal_up_bounds(m, rz, t):
        return (0, m.ZonalRegulationUpRequirement[rz,t])
    model.ZonalRegulationUpShortfall = Var(model.RegulationZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_up_bounds)
    def zonal_dn_bounds(m, rz, t):
        return (0, m.ZonalRegulationDnRequirement[rz,t])
    model.ZonalRegulationDnShortfall = Var(model.RegulationZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_dn_bounds)

    def system_up_bounds(m, t):
        return (0, m.SystemRegulationUpRequirement[t])
    model.SystemRegulationUpShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_up_bounds)
    def system_dn_bounds(m, t):
        return (0, m.SystemRegulationDnRequirement[t])
    model.SystemRegulationDnShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_dn_bounds)
    
    # regulation cost for
    model.RegulationOfferFixedCost = Param(model.AGC_Generators, model.TimePeriods, within=NonNegativeReals, default=0.0, initialize=TimeMapper(agc_gen_attrs.get('agc_fixed_cost', dict())))
    model.RegulationOfferMarginalCost = Param(model.AGC_Generators, model.TimePeriods, within=NonNegativeReals, default=0.0, initialize=TimeMapper(agc_gen_attrs.get('agc_marginal_cost', dict())))

    if _is_relaxed(model):
        model.RegulationOn = Var(model.AGC_Generators, model.TimePeriods, within=UnitInterval)
    else:
        model.RegulationOn = Var(model.AGC_Generators, model.TimePeriods, within=Binary)

    def reg_up_bounds(m,g,t):
        return (0, m.RegulationUpCapability[g,t])
    def reg_dn_bounds(m,g,t):
        return (0, m.RegulationDnCapability[g,t])
    model.RegulationReserveUp = Var(model.AGC_Generators, model.TimePeriods, within=NonNegativeReals, bounds=reg_up_bounds)
    model.RegulationReserveDn = Var(model.AGC_Generators, model.TimePeriods, within=NonNegativeReals, bounds=reg_dn_bounds)

    def enforce_fixed_regulation_rule(m,g,t):
        if value(m.FixedRegulation[g,t]) is not None:
            m.RegulationOn[g,t].value = value(m.FixedRegulation[g,t])
            m.RegulationOn[g,t].fix()
    model.EnforceFixedRegulation = BuildAction(model.AGC_Generators, model.TimePeriods, rule=enforce_fixed_regulation_rule)

    # a generator can provide regulation only when it's on
    def provide_regulation_when_unit_on_rule(m, g, t):
        return m.RegulationOn[g, t] <= m.UnitOn[g, t]
    model.EnforceRegulationOnWhenUnitOn = Constraint(model.AGC_Generators, model.TimePeriods, rule=provide_regulation_when_unit_on_rule)

    def reg_up_rule(m,g,t):
        reg_up_limit = min(value(m.RegulationUpCapability[g,t]), value(m.NominalRampUpLimit[g]/60.*m.RegulationMinutes))
        return m.RegulationReserveUp[g,t] <= reg_up_limit*m.RegulationOn[g,t]
    model.EnforceRegulationUpBound = Constraint(model.AGC_Generators, model.TimePeriods, rule=reg_up_rule)

    def reg_dn_rule(m,g,t):
        reg_dn_limit = min(value(m.RegulationDnCapability[g,t]), value(m.NominalRampDownLimit[g]/60.*m.RegulationMinutes))
        return m.RegulationReserveDn[g,t] <= reg_dn_limit*m.RegulationOn[g,t]
    model.EnforceRegulationDnBound = Constraint(model.AGC_Generators, model.TimePeriods, rule=reg_dn_rule)

    def zonal_reg_up_provided(m,rz,t):
        return sum(m.RegulationReserveUp[g,t] for g in m.AGC_GeneratorsInRegulationZone[rz]) + m.ZonalRegulationUpShortfall[rz,t] 
    model.ZonalRegulationUpProvided = Expression(model.RegulationZones, model.TimePeriods, rule=zonal_reg_up_provided)

    def enforce_zonal_reg_up_requirement_rule(m, rz, t):
        return  m.ZonalRegulationUpProvided[rz,t] >= m.ZonalRegulationUpRequirement[rz,t]
    model.EnforceZonalRegulationUpRequirements = Constraint(model.RegulationZones, model.TimePeriods, rule=enforce_zonal_reg_up_requirement_rule)

    def enforce_zonal_reg_dn_requirement_rule(m, rz, t):
        return sum(m.RegulationReserveDn[g,t] for g in m.AGC_GeneratorsInRegulationZone[rz]) + \
                m.ZonalRegulationDnShortfall[rz,t] >= m.ZonalRegulationDnRequirement[rz,t]
    model.EnforceZonalRegulationDnRequirements = Constraint(model.RegulationZones, model.TimePeriods, rule=enforce_zonal_reg_dn_requirement_rule)

    ## NOTE: making sure not to double count the shortfall
    def system_reg_up_provided(m,t):
        return sum(m.RegulationReserveUp[g,t] for g in m.AGC_Generators) + \
                m.SystemRegulationUpShortfall[t] + sum(m.ZonalRegulationUpShortfall[rz,t] for rz in m.RegulationZones) 
    model.SystemRegulationUpProvided = Expression(model.TimePeriods, rule=system_reg_up_provided)

    def enforce_system_regulation_up_requirement_rule(m, t):
        return m.SystemRegulationUpProvided[t] >= m.SystemRegulationUpRequirement[t]
    model.EnforceSystemRegulationUpRequirement = Constraint(model.TimePeriods, rule=enforce_system_regulation_up_requirement_rule)

    def enforce_system_regulation_dn_requirement_rule(m, t):
        return sum(m.RegulationReserveDn[g,t] for g in m.AGC_Generators) + \
                m.SystemRegulationDnShortfall[t] + sum(m.ZonalRegulationDnShortfall[rz,t] for rz in m.RegulationZones) \
                >= m.SystemRegulationDnRequirement[t]
    model.EnforceSystemRegulationDnRequirement = Constraint(model.TimePeriods, rule=enforce_system_regulation_dn_requirement_rule)

    def regulation_cost_commitment(m,g,t):
        return m.RegulationOfferFixedCost[g,t] * m.RegulationOn[g, t]*m.TimePeriodLengthHours
    model.RegulationCostCommitment = Expression(model.AGC_Generators, model.TimePeriods, rule=regulation_cost_commitment)

    def regulation_cost_generation(m,g,t):
        return m.RegulationOfferMarginalCost[g,t]*m.TimePeriodLengthHours*(m.RegulationReserveUp[g,t] + m.RegulationReserveDn[g,t])
    model.RegulationCostGeneration = Expression(model.AGC_Generators, model.TimePeriods, rule=regulation_cost_generation)

    def regulation_cost_slacks(m,t):
        return m.TimePeriodLengthHours*m.RegulationPenalty*(
                        m.SystemRegulationUpShortfall[t] + m.SystemRegulationDnShortfall[t] \
                      + sum(m.ZonalRegulationUpShortfall[rz,t] for rz in m.RegulationZones) \
                      + sum(m.ZonalRegulationDnShortfall[rz,t] for rz in m.RegulationZones) \
                      )
    model.RegulationCostPenalty = Expression(model.TimePeriods, rule=regulation_cost_slacks)

    ## end regulation_services

@add_model_attr('spinning_reserve', requires = {'data_loader': None,
                                                'status_vars': None,
                                               })
def spinning_reserves(model, zone_initializer_builder, zone_requirement_getter, gens_in_reserve_zone_getter, thermal_gen_attrs):

    md = model.model_data

    system = md.data['system']

    TimeMapper = uc_time_helper(model.TimePeriods)

    def _check_spin(e_dict):
        return 'spinning_reserve_requirement' in e_dict

    model.SpinningReserveZones = Set(initialize=zone_initializer_builder(_check_spin))

    model.ThermalGeneratorsInSpinningReserveZone = Set(model.SpinningReserveZones, initialize=gens_in_reserve_zone_getter())
    ## begin spinning reserve

    # spinning reserve response time
    model.SpinningReserveMinutes = Param(within=PositiveReals, default=10.) # in minutes, varies among ISOs

    # limit,  cost of spinning reserves
    # NOTE: This is here in case the user wants to limit this beyond the ramping limits
    model.SpinningReserveCapability = Param(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, default=float('inf'),
                                                initialize=TimeMapper(thermal_gen_attrs.get('spinning_capacity', dict())))
    model.SpinningReservePrice = Param(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, default=0.0, initialize=TimeMapper(thermal_gen_attrs.get('spinning_cost', dict())))
    
    # spinning reserve requirements
    model.ZonalSpinningReserveRequirement = Param(model.SpinningReserveZones, model.TimePeriods, within=NonNegativeReals,
                                                        initialize=zone_requirement_getter('spinning_reserve_requirement'))
    model.SystemSpinningReserveRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0, initialize=TimeMapper(system.get('spinning_reserve_requirement', dict())))

    def zonal_spin_bounds(m,rz,t):
        return (0, m.ZonalSpinningReserveRequirement[rz,t])
    model.ZonalSpinningReserveShortfall = Var(model.SpinningReserveZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_spin_bounds)
    def system_spin_bounds(m,t):
        return (0, m.SystemSpinningReserveRequirement[t])
    model.SystemSpinningReserveShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_spin_bounds)

    # spinning reserve
    def spin_bounds(m,g,t):
        return (0,m.SpinningReserveCapability[g,t])
    model.SpinningReserveDispatched = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=spin_bounds)

    regup_reserves = bool(model.regulation_service)
    def spinning_reserve_available(m, g, t):
        spin_limit = min(value(m.SpinningReserveCapability[g,t]), value(m.NominalRampUpLimit[g]/60.*m.SpinningReserveMinutes))
        if regup_reserves and g in m.AGC_Generators:
            return m.RegulationReserveUp[g,t] + m.SpinningReserveDispatched[g, t] <= spin_limit*m.UnitOn[g,t]
        else:
            return m.SpinningReserveDispatched[g, t] <= spin_limit*m.UnitOn[g,t]
    model.SpinningReserveAvailableConstr = Constraint(model.ThermalGenerators, model.TimePeriods, rule=spinning_reserve_available)

    def zonal_spinning_reserve_provided(m, rz, t):
        return sum(m.SpinningReserveDispatched[g, t] for g in m.ThermalGeneratorsInSpinningReserveZone[rz])\
                + m.ZonalSpinningReserveShortfall[rz,t]
    model.ZonalSpinningReserveProvided = Expression(model.SpinningReserveZones, model.TimePeriods, rule=zonal_spinning_reserve_provided)

    def enforce_zonal_spinning_reserve_requirement(m, rz, t):
        if regup_reserves and (rz in m.RegulationZones):
            return  m.ZonalSpinningReserveProvided[rz,t] + m.ZonalRegulationUpProvided[rz,t] \
                    >= m.ZonalSpinningReserveRequirement[rz, t] + m.ZonalRegulationUpRequirement[rz, t]
        else:
            return  m.ZonalSpinningReserveProvided[rz,t] >= m.ZonalSpinningReserveRequirement[rz, t]
    model.EnforceZonalSpinningReserveRequirement = Constraint(model.SpinningReserveZones, model.TimePeriods, rule=enforce_zonal_spinning_reserve_requirement)

    def system_spinning_reserve_provided(m,t):
        return sum(m.SpinningReserveDispatched[g,t] for g in m.ThermalGenerators) \
                + sum(m.ZonalSpinningReserveShortfall[rz,t] for rz in m.SpinningReserveZones) \
                + m.SystemSpinningReserveShortfall[t]
    model.SystemSpinningReserveProvided = Expression(model.TimePeriods, rule=system_spinning_reserve_provided)

    def enforce_system_spinning_reserve_requirement(m, t):
        if regup_reserves:
            return m.SystemSpinningReserveProvided[t] + m.SystemRegulationUpProvided[t] \
                    >= m.SystemSpinningReserveRequirement[t] + m.SystemRegulationUpRequirement[t]
        else:
            return m.SystemSpinningReserveProvided[t] >= m.SystemSpinningReserveRequirement[t]
    model.EnforceSystemSpinningReserveRequirement = Constraint(model.TimePeriods, rule=enforce_system_spinning_reserve_requirement)

    def compute_spinning_reserve_cost(m, g, t):
        return m.SpinningReserveDispatched[g, t] * m.SpinningReservePrice[g,t] * m.TimePeriodLengthHours
    model.SpinningReserveCostGeneration = Expression(model.ThermalGenerators, model.TimePeriods, rule=compute_spinning_reserve_cost)

    def spinning_reserve_cost_slacks(m,t):
        return m.TimePeriodLengthHours*m.SpinningReservePenalty*(
                      m.SystemSpinningReserveShortfall[t] \
                    + sum(m.ZonalSpinningReserveShortfall[rz,t] for rz in m.SpinningReserveZones)
                    )
    model.SpinningReserveCostPenalty = Expression(model.TimePeriods, rule=spinning_reserve_cost_slacks)

    ## end spinning reserves


@add_model_attr('non_spinning_reserve', requires = {'data_loader': None,
                                                  'status_vars': None,
                                                   })
def non_spinning_reserves(model, zone_initializer_builder, zone_requirement_getter, gens_in_reserve_zone_getter):

    md = model.model_data

    system = md.data['system']

    TimeMapper = uc_time_helper(model.TimePeriods)

    def _check_nspin(e_dict):
        return 'non_spinning_reserve_requirement' in e_dict
    model.NonSpinReserveZones = Set(initialize=zone_initializer_builder(_check_nspin))

    ## begin non-spinning reserves
    nspin_gen_attrs = md.attributes(element_type='generator', generator_type='thermal', fast_start=True)

    model.NonSpinGenerators = Set(within=model.ThermalGenerators, initialize=nspin_gen_attrs['names'])

    model.NonSpinGeneratorsInNonSpinZone = Set(model.NonSpinReserveZones, initialize=gens_in_reserve_zone_getter(nspin_gen_attrs))
    
    # Non-spinning reserves are assumed to be fast -- Supplemental reserves are slow (30 minutes)

    def validate_nonspin_bid(m,v,g,t):
        return v <= value(m.MaximumPowerOutput[g,t])
    model.NonSpinningReserveCapability = Param(model.NonSpinGenerators, model.TimePeriods, within=NonNegativeReals, default=0.0, validate=validate_nonspin_bid,
                                                    initialize=TimeMapper(nspin_gen_attrs['non_spinning_capacity']))
    model.NonSpinningReservePrice = Param(model.NonSpinGenerators, model.TimePeriods, within=NonNegativeReals, default=0.0, initialize=TimeMapper(nspin_gen_attrs.get('non_spinning_cost', dict())))
    
    model.ZonalNonSpinningReserveRequirement = Param(model.NonSpinReserveZones, model.TimePeriods, within=NonNegativeReals, default=0.0,
                                                        initialize=zone_requirement_getter('non_spinning_reserve_requirement'))
    model.SystemNonSpinningReserveRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0, 
                                                        initialize=TimeMapper(system.get('non_spinning_reserve_requirement', dict())))

    def zonal_fast_bounds(m,rz,t):
        return (0, m.ZonalNonSpinningReserveRequirement[rz,t])
    model.ZonalNonSpinningReserveShortfall = Var(model.NonSpinReserveZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_fast_bounds)
    def system_fast_bounds(m,t):
        return (0, m.SystemNonSpinningReserveRequirement[t])
    model.SystemNonSpinningReserveShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_fast_bounds)

    def nspin_bounds(m,g,t):
        return (0,m.NonSpinningReserveCapability[g,t])
    model.NonSpinningReserveDispatched = Var(model.NonSpinGenerators, model.TimePeriods, within=NonNegativeReals, bounds=nspin_bounds)

    # non-spinning reserve
    def calculate_non_spinning_reserve_limit_rule(m, g, t):
        return m.NonSpinningReserveDispatched[g, t] <= m.NonSpinningReserveCapability[g,t] * (1 - m.UnitOn[g, t])
    model.CalculateNonSpinningReserveLimit = Constraint(model.NonSpinGenerators, model.TimePeriods, rule=calculate_non_spinning_reserve_limit_rule)

    def nspin_zonal_reserves_provided(m,rz,t):
        return sum(m.NonSpinningReserveDispatched[g,t] for g in m.NonSpinGeneratorsInNonSpinZone[rz]) \
                    + m.ZonalNonSpinningReserveShortfall[rz,t]
    model.NonSpinningZonalReservesProvided = Expression(model.NonSpinReserveZones, model.TimePeriods, rule=nspin_zonal_reserves_provided)

    spin_reserves = bool(model.spinning_reserve)
    regup_reserves = bool(model.regulation_service)

    def enforce_zonal_non_spinning_reserve_rule(m, rz, t):
        zonal_spin_reserves = (spin_reserves and rz in m.SpinningReserveZones)
        zonal_regup_reserves = (regup_reserves and rz in m.RegulationZones)
        return m.NonSpinningZonalReservesProvided[rz,t] \
                + (m.ZonalSpinningReserveProvided[rz,t] if zonal_spin_reserves else 0.) \
                + (m.ZonalRegulationUpProvided[rz,t] if zonal_regup_reserves else 0.) \
               >= m.ZonalNonSpinningReserveRequirement[rz, t] \
                + (m.ZonalSpinningReserveRequirement[rz,t] if zonal_spin_reserves else 0.) \
                + (m.ZonalRegulationUpRequirement[rz,t] if zonal_regup_reserves else 0.)
    model.EnforceNonSpinningZonalReserveRequirement = Constraint(model.NonSpinReserveZones, model.TimePeriods, rule=enforce_zonal_non_spinning_reserve_rule)

    def nspin_reserves_provided(m,t):
        return sum(m.NonSpinningReserveDispatched[g,t] for g in m.NonSpinGenerators) \
                + sum(m.ZonalNonSpinningReserveShortfall[rz,t] for rz in m.NonSpinReserveZones) \
                + m.SystemNonSpinningReserveShortfall[t]
    model.SystemNonSpinningReserveProvided = Expression(model.TimePeriods, rule=nspin_reserves_provided)

    def enforce_system_non_spinning_reserve_requirement(m, t):
        return m.SystemNonSpinningReserveProvided[t] \
                  + (m.SystemSpinningReserveProvided[t] if spin_reserves else 0.) \
                  + (m.SystemRegulationUpProvided[t] if regup_reserves else 0.) \
                >= m.SystemNonSpinningReserveRequirement[t] \
                  + (m.SystemSpinningReserveRequirement[t] if spin_reserves else 0.) \
                  + (m.SystemRegulationUpRequirement[t] if regup_reserves else 0.)
    model.EnforceSystemNonSpinningReserveRequirement = Constraint(model.TimePeriods, rule=enforce_system_non_spinning_reserve_requirement)

    def calculate_non_spinning_reserve_cost(m, g, t):
        return m.NonSpinningReserveDispatched[g, t] * m.NonSpinningReservePrice[g,t] * m.TimePeriodLengthHours
    model.NonSpinningReserveCostGeneration = Expression(model.NonSpinGenerators, model.TimePeriods, rule=calculate_non_spinning_reserve_cost)

    def non_spinning_reserve_cost_penalty(m,t):
        return m.TimePeriodLengthHours*m.NonSpinningReservePenalty*(
                        m.SystemNonSpinningReserveShortfall[t] \
                      + sum(m.ZonalNonSpinningReserveShortfall[rz,t] for rz in m.NonSpinReserveZones)
                      )
    model.NonSpinningReserveCostPenalty = Expression(model.TimePeriods, rule=non_spinning_reserve_cost_penalty)

    ## end non-spinning reserve


@add_model_attr('supplemental_reserve', requires = {'data_loader': None,
                                                  'status_vars': None,
                                                 })
def supplemental_reserves(model, zone_initializer_builder, zone_requirement_getter, gens_in_reserve_zone_getter, thermal_gen_attrs):

    md = model.model_data

    system = md.data['system']

    TimeMapper = uc_time_helper(model.TimePeriods)

    def _check_supplemental(e_dict):
        return 'supplemental_reserve_requirement' in e_dict
    model.SupplementalReserveZones = Set(initialize=zone_initializer_builder(_check_supplemental))

    supplemental_nspin_gen_attrs = md.attributes(element_type='generator', generator_type='thermal', supplemental_start=True)

    model.SupplementalNonSpinGenerators = Set(within=model.ThermalGenerators, initialize=supplemental_nspin_gen_attrs['names'])

    model.GeneratorsInSupplementalReserveZone = Set(model.SupplementalReserveZones, initialize=gens_in_reserve_zone_getter())

    ## begin supplemental reserve

    # Thirty-minute supplemental reserves, for generators which can start in 30 minutes
    def validate_nonspin_bid(m,v,g,t):
        return v <= value(m.MaximumPowerOutput[g,t])
    model.SupplementalReserveCapabilityNonSpin = Param(model.SupplementalNonSpinGenerators, model.TimePeriods, within=NonNegativeReals, default=0.0, 
                                                        validate=validate_nonspin_bid, initialize=TimeMapper(supplemental_nspin_gen_attrs.get('supplemental_non_spinning_capacity', dict())))

    ## NOTE: this param is here if the user wants to limit this beyond the nominal ramping limits
    model.SupplementalReserveCapabilitySpin = Param(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, default=float('inf'),
                                                        initialize=TimeMapper(thermal_gen_attrs.get('supplemental_spinning_capacity', dict())))

    model.SupplementalReservePrice = Param(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, default=0.0, initialize=TimeMapper(thermal_gen_attrs.get('supplemental_cost', dict())))
    model.SupplementalReserveMinutes = Param(within=PositiveReals, default=30.)

    # Supplemental reserve requirement

    model.ZonalSupplementalReserveRequirement = Param(model.SupplementalReserveZones, model.TimePeriods, within=NonNegativeReals, default=0.0,
                                                        initialize=zone_requirement_getter('supplemental_reserve_requirement'))
    model.SystemSupplementalReserveRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0, 
                                                        initialize=TimeMapper(system.get('supplemental_reserve_requirement', dict())))

    def zonal_op_bounds(m,rz,t):
        return (0, m.ZonalSupplementalReserveRequirement[rz,t])
    model.ZonalSupplementalReserveShortfall = Var(model.SupplementalReserveZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_op_bounds)
    def system_op_bounds(m,t):
        return (0, m.SystemSupplementalReserveRequirement[t])
    model.SystemSupplementalReserveShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_op_bounds)
    
    def op_bounds_nspin(m,g,t):
        return (0,m.SupplementalReserveCapabilityNonSpin[g,t])
    model.SupplementalNonSpinReserveDispatched = Var(model.SupplementalNonSpinGenerators, model.TimePeriods, within=NonNegativeReals, bounds=op_bounds_nspin)

    def op_bounds_spin(m,g,t):
        return (0,m.SupplementalReserveCapabilitySpin[g,t])
    model.SupplementalSpinReserveDispatched = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, bounds=op_bounds_spin)

    nspin_reserves = bool(model.non_spinning_reserve)
    # thirty-minute supplemental reserve, for units which are off
    def calculate_supplemental_reserve_limit_rule_nonspin(m, g, t):
        if nspin_reserves and g in m.NonSpinGenerators:
            return m.SupplementalNonSpinReserveDispatched[g, t] + m.NonSpinningReserveDispatched[g, t] <= m.SupplementalReserveCapabilityNonSpin[g,t] * (1 - m.UnitOn[g, t])
        else:
            return m.SupplementalNonSpinReserveDispatched[g, t] <= m.SupplementalReserveCapabilityNonSpin[g,t] * (1 - m.UnitOn[g, t])
    model.CalculateSupplementalReserveLimits = Constraint(model.SupplementalNonSpinGenerators, model.TimePeriods, rule=calculate_supplemental_reserve_limit_rule_nonspin)

    regup_reserves = bool(model.regulation_service)
    spin_reserves = bool(model.spinning_reserve)
    def calculate_supplemental_reserve_limit_rule_spin(m, g, t):
        spin_limit = min(value(m.SupplementalReserveCapabilitySpin[g,t]), value(m.NominalRampUpLimit[g]/60.*m.SupplementalReserveMinutes))
        regup = (regup_reserves and (g in m.AGC_Generators))
        return m.SupplementalSpinReserveDispatched[g,t] \
                 + (m.SpinningReserveDispatched[g, t] if spin_reserves else 0.) \
                 + (m.RegulationReserveUp[g,t] if regup else 0.) \
               <= spin_limit*m.UnitOn[g,t]
    model.CalculateSupplementalReserveLimitsSpin = Constraint(model.ThermalGenerators, model.TimePeriods, rule=calculate_supplemental_reserve_limit_rule_spin)

    def supplemental_reserve_expr_rule(m, g, t):
        if g in m.SupplementalNonSpinGenerators:
            return m.SupplementalNonSpinReserveDispatched[g,t] + m.SupplementalSpinReserveDispatched[g,t]
        else:
            return m.SupplementalSpinReserveDispatched[g,t]
    model.SupplementalReserveDispatched = Expression(model.ThermalGenerators, model.TimePeriods, rule=supplemental_reserve_expr_rule)

    def operational_zonal_reserves_provided(m,rz,t):
        return sum(m.SupplementalReserveDispatched[g,t] for g in m.GeneratorsInSupplementalReserveZone[rz]) + m.ZonalSupplementalReserveShortfall[rz,t]
    model.SupplementalZonalReservesProvided = Expression(model.SupplementalReserveZones, model.TimePeriods, rule=operational_zonal_reserves_provided)

    def enforce_zonal_supplemental_reserve_requirement_rule(m, rz, t):
        reg_up = (regup_reserves and (rz in m.RegulationZones))
        nspin = (nspin_reserves and (rz in m.NonSpinReserveZones))
        spin = (spin_reserves and (rz in m.SpinningReserveZones))
        return m.SupplementalZonalReservesProvided[rz,t] \
                  + (m.NonSpinningZonalReservesProvided[rz,t] if nspin else 0.) \
                  + (m.ZonalSpinningReserveProvided[rz,t] if spin else 0.) \
                  + (m.ZonalRegulationUpRequirement[rz,t] if reg_up else 0.) \
                >= m.ZonalSupplementalReserveRequirement[rz,t] \
                  + (m.ZonalNonSpinningReserveRequirement[rz,t] if nspin else 0.) \
                  + (m.ZonalSpinningReserveRequirement[rz,t] if spin else 0.)\
                  + (m.ZonalRegulationUpRequirement[rz,t] if reg_up else 0.)
    model.EnforceZonalSupplementalReserveRequirement = Constraint(model.SupplementalReserveZones, model.TimePeriods, rule=enforce_zonal_supplemental_reserve_requirement_rule)

    def operational_reserves_provided(m,t):
        return sum(m.SupplementalReserveDispatched[g,t] for g in m.ThermalGenerators) \
                + sum(m.ZonalSupplementalReserveShortfall[rz,t] for rz in m.SupplementalReserveZones) \
                + m.SystemSupplementalReserveShortfall[t]
    model.SystemSupplementalReserveProvided = Expression(model.TimePeriods, rule=operational_reserves_provided)

    def enforce_system_supplemental_reserve_requirement(m, t):
        return m.SystemSupplementalReserveProvided[t] \
                    + (m.SystemNonSpinningReserveProvided[t] if nspin_reserves else 0.) \
                    + (m.SystemSpinningReserveProvided[t] if spin_reserves else 0.) \
                    + (m.SystemRegulationUpRequirement[t] if regup_reserves else 0.)\
                >= m.SystemSupplementalReserveRequirement[t] \
                    + (m.SystemNonSpinningReserveRequirement[t] if nspin_reserves else 0.) \
                    + (m.SystemSpinningReserveRequirement[t] if spin_reserves else 0.)\
                    + (m.SystemRegulationUpRequirement[t] if regup_reserves else 0.)

    model.EnforceSystemSupplementalReserveRequirement = Constraint(model.TimePeriods, rule=enforce_system_supplemental_reserve_requirement)

    def calculate_supplemental_reserve_cost_rule(m, g, t):
        return m.SupplementalReserveDispatched[g, t] * m.SupplementalReservePrice[g,t] * m.TimePeriodLengthHours
    model.SupplementalReserveCostGeneration = Expression(model.ThermalGenerators, model.TimePeriods, rule=calculate_supplemental_reserve_cost_rule)

    def supplemental_reserve_cost_penalty(m,t):
        return m.TimePeriodLengthHours*m.SupplementalReservePenalty*(
                        m.SystemSupplementalReserveShortfall[t] \
                      + sum(m.ZonalSupplementalReserveShortfall[rz,t] for rz in m.SupplementalReserveZones)
                      )
    model.SupplementalReserveCostPenalty = Expression(model.TimePeriods, rule=supplemental_reserve_cost_penalty)

    ## end supplemental reserve

@add_model_attr('flexible_ramping', requires = {'data_loader': None,
                                                  'status_vars': None,
                                                 })
def flexible_ramping_reserves(model, zone_initializer_builder, zone_requirement_getter, gens_in_reserve_zone_getter, thermal_gen_attrs):

    md = model.model_data

    system = md.data['system']

    TimeMapper = uc_time_helper(model.TimePeriods)

    def _check_flex(e_dict):
        return ( ('flexible_ramp_up_requirement' in e_dict) \
                  or ('flexible_ramp_down_requirement' in e_dict) )
    model.FlexRampZones = Set(initialize=zone_initializer_builder(_check_flex))

    model.ThermalGeneratorsInFlexRampZone = Set(model.FlexRampZones, initialize=gens_in_reserve_zone_getter())

    ## begin flexible_ramp
    model.FlexRampMinutes = Param(within=PositiveReals, default=20.)

    model.ZonalFlexUpRequirement = Param(model.FlexRampZones, model.TimePeriods, within=NonNegativeReals, default=0.0,
                                            initialize=zone_requirement_getter('flexible_ramp_up_requirement'))
    model.ZonalFlexDnRequirement = Param(model.FlexRampZones, model.TimePeriods, within=NonNegativeReals, default=0.0,
                                            initialize=zone_requirement_getter('flexible_ramp_down_requirement'))

    model.SystemFlexUpRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0, initialize=TimeMapper(system.get('flexible_ramp_up_requirement', dict())))
    model.SystemFlexDnRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0, initialize=TimeMapper(system.get('flexible_ramp_down_requirement', dict())))

    def zonal_flex_up_bounds(m, rz, t):
        return (0, m.ZonalFlexUpRequirement[rz,t])
    model.ZonalFlexUpShortfall = Var(model.FlexRampZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_flex_up_bounds)
    def zonal_flex_dn_bounds(m, rz, t):
        return (0, m.ZonalFlexDnRequirement[rz,t])
    model.ZonalFlexDnShortfall = Var(model.FlexRampZones, model.TimePeriods, within=NonNegativeReals, bounds=zonal_flex_dn_bounds)

    def system_flex_up_bounds(m, t):
        return (0, m.SystemFlexUpRequirement[t])
    model.SystemFlexUpShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_flex_up_bounds)
    def system_flex_dn_bounds(m, t):
        return (0, m.SystemFlexDnRequirement[t])
    model.SystemFlexDnShortfall = Var(model.TimePeriods, within=NonNegativeReals, bounds=system_flex_dn_bounds)


    model.FlexUpProvided = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)
    model.FlexDnProvided = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)

    def flex_up_limit_rule(m, g, t):
        return m.FlexUpProvided[g,t] <= m.FlexRampMinutes*(m.NominalRampUpLimit[g]/60.)*m.UnitOn[g,t]
    model.FlexUpLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=flex_up_limit_rule)

    def flex_down_limit_rule(m, g, t):
        return m.FlexDnProvided[g,t] <= m.FlexRampMinutes*(m.NominalRampDownLimit[g]/60.)*m.UnitOn[g,t]
    model.FlexDnLimit = Constraint(model.ThermalGenerators, model.TimePeriods, rule=flex_down_limit_rule)

    def zonal_flex_up_requirement_rule(m, rz, t):
        return sum(m.FlexUpProvided[g,t] for g in m.ThermalGeneratorsInFlexRampZone[rz]) \
                    + m.ZonalFlexUpShortfall[rz,t] >= m.ZonalFlexUpRequirement[rz,t]
    model.ZonalFlexUpRequirementConstr = Constraint(model.FlexRampZones, model.TimePeriods, rule=zonal_flex_up_requirement_rule)

    def zonal_flex_dn_requirement_rule(m, rz, t):
        return sum(m.FlexDnProvided[g,t] for g in m.ThermalGeneratorsInFlexRampZone[rz]) \
                    + m.ZonalFlexDnShortfall[rz,t] >= m.ZonalFlexDnRequirement[rz,t]
    model.ZonalFlexDnRequirementConstr = Constraint(model.FlexRampZones, model.TimePeriods, rule=zonal_flex_dn_requirement_rule)

    def system_flex_up_requirement_rule(m, t):
        return sum(m.FlexUpProvided[g,t] for g in m.ThermalGenerators) \
                 + sum(m.ZonalFlexUpShortfall[rz,t] for rz in m.FlexRampZones) \
                 + m.SystemFlexUpShortfall[t] \
                 >= m.SystemFlexUpRequirement[t]
    model.SystemFlexUpRequirementConstr = Constraint(model.TimePeriods, rule=system_flex_up_requirement_rule)

    def system_flex_dn_requirement_rule(m, t):
        return sum(m.FlexDnProvided[g,t] for g in m.ThermalGenerators) \
                 + sum(m.ZonalFlexDnShortfall[rz,t] for rz in m.FlexRampZones) \
                 + m.SystemFlexDnShortfall[t] \
                 >= m.SystemFlexDnRequirement[t]
    model.SystemFlexDnRequirementConstr = Constraint(model.TimePeriods, rule=system_flex_dn_requirement_rule)
    
    def flex_ramp_penalty_cost(m, t):
        return m.TimePeriodLengthHours*m.FlexRampPenalty*(
                        m.SystemFlexUpShortfall[t] + m.SystemFlexDnShortfall[t] \
                      + sum(m.ZonalFlexUpShortfall[rz,t]+m.ZonalFlexDnShortfall[rz,t] for rz in m.FlexRampZones)
                      )
    model.FlexibleRampingCostPenalty = Expression(model.TimePeriods, rule = flex_ramp_penalty_cost)


    ## end flexible ramp

