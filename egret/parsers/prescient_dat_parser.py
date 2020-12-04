#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## loads and validates input unit commitment data
from pyomo.environ import *
import os.path
import egret.data.model_data as md

def create_ModelData(dat_file):
    '''
    Create a ModelData object from a prescient dat file

    Parameters
    ----------
    dat_file : str
        Path to prescient *.dat file

    Returns
    -------
    egret.model_data.ModelData
        Returns a ModelData object with the dat file data specified
    '''
    return md.ModelData(create_model_data_dict(dat_file))

def create_model_data_dict(dat_file):
    model = get_uc_model()
    params = model.create_instance(dat_file)
    return create_model_data_dict_params(params)

def get_uc_model():
    uc_model = AbstractModel()
    setup_abstract_model(uc_model)
    return uc_model

def create_model_data_dict_params(params, keep_names=False):

    md_dict = md.ModelData.empty_model_data_dict()

    elements = md_dict['elements']
    system = md_dict['system']

    system['time_keys'] = list(str(t) for t in params.TimePeriods)
    system['time_period_length_minutes'] = value(params.TimePeriodLengthMinutes)

    system['load_mismatch_cost'] = value(params.LoadMismatchPenalty)
    system['reserve_shortfall_cost'] = value(params.ReserveShortfallPenalty)

    ## These UC param files have the baseMVA factored out
    system['baseMVA'] = 1.

    bus_dict = dict() 
    set_reference = True
    gen_bus_dict = dict()
    renewable_gen_bus_dict = dict()
    storage_bus_dict = dict()
    for b in params.Buses:
        if set_reference:
            system['reference_bus'] = b
            system['reference_bus_angle'] = 0.0
            set_reference = False
        bus_dict[b] = {'base_kv' : params.BusKV[b]}
        for g in params.ThermalGeneratorsAtBus[b]:
            gen_bus_dict[g] = b
        for n in params.NondispatchableGeneratorsAtBus[b]:
            renewable_gen_bus_dict[n] = b
        for s in params.StorageAtBus[b]:
            storage_bus_dict[s] = b
    elements['bus'] = bus_dict

    load_dict = dict()
    for b in params.Buses:
        l_d = { 'bus' : b, 
                'in_service': True,
                'p_load':
                        {'data_type':'time_series',
                            'values': [params.Demand[b,t] for t in params.TimePeriods ]
                        }
               }
        load_dict[b] = l_d
    elements['load'] = load_dict

    reserve_dict = [ value(params.ReserveRequirement[t]) for t in params.TimePeriods ]
    system['reserve_requirement'] = { 'data_type':'time_series', 'values': reserve_dict }

    branch_dict = dict()
    for l in params.TransmissionLines:
        b_d = { 'from_bus' : params.BusFrom[l],
                'to_bus' : params.BusTo[l],
                'reactance' : params.Impedence[l],
                'rating_long_term' : params.ThermalLimit[l],
                'rating_short_term' : params.ThermalLimit[l],
                'rating_emergency' : params.ThermalLimit[l],
                'in_service' : True,
                'branch_type' : 'line',
                'angle_diff_min': -90,
                'angle_diff_max': 90,
                }
        branch_dict[l] = b_d
    elements['branch'] = branch_dict

    interface_dict = dict()
    for i in params.Interfaces:
        i_d = { 'interface_lines' : list(params.InterfaceLines[l]),
                'interface_from_limit': params.InterfaceFromLimit[l],
                'interface_to_limit': params.InterfaceToLimit[l],
              }
        interface_dict[i] = i_d
    elements['interface'] = interface_dict

    zone_dict = dict()
    for z in params.ReserveZones:
        reserve_dict = [ params.ZonalReserveRequirement[z,t] for t in params.TimePeriods ]
        z_d = { 'reserve_requirement' : {'data_type': 'time_series', 'values' : reserve_dict } }
        zone_dict[z] = z_d
    elements['zone'] = zone_dict

    gen_dict = dict()

    for g in params.ThermalGenerators:
        g_d = { 'generator_type':'thermal', }
        g_d['bus'] = gen_bus_dict[g]
        g_d['fuel'] = params.ThermalGeneratorType[g]
        g_d['fast_start'] = (g in params.QuickStartGenerators)
        g_d['fixed_commitment'] = (1 if params.MustRun[g] else None)
        g_d['in_service'] = True
        g_d['zone'] = params.ReserveZoneLocation[g]
        g_d['failure_rate'] = params.FailureProbability[g]
        g_d['p_min'] = params.MinimumPowerOutput[g]
        g_d['p_max'] = params.MaximumPowerOutput[g]
        g_d['ramp_up_60min'] = params.NominalRampUpLimit[g]
        g_d['ramp_down_60min'] = params.NominalRampDownLimit[g]
        g_d['startup_capacity'] = params.StartupRampLimit[g]
        g_d['shutdown_capacity'] = params.ShutdownRampLimit[g]
        g_d['min_up_time'] = params.MinimumUpTime[g]
        g_d['min_down_time'] = params.MinimumDownTime[g]
        g_d['initial_status'] = params.UnitOnT0State[g]
        g_d['initial_p_output'] = params.PowerGeneratedT0[g]
        g_d['startup_cost'] = list(zip(params.StartupLags[g],params.StartupCosts[g]))
        g_d['shutdown_cost'] = params.ShutdownFixedCost[g]
        p_cost = {'data_type' : 'cost_curve' }
        if value(params.PiecewiseType) == "NoPiecewise":
            p_cost['cost_curve_type'] = 'polynomial'
            p_cost['values'] = { 0 : params.FuelCost[g]*params.ProductionCostA0[g],
                                 1 : params.FuelCost[g]*params.ProductionCostA1[g],
                                 2 : params.FuelCost[g]*params.ProductionCostA2[g],
                               }
        else:
            p_cost['cost_curve_type'] = 'piecewise'
            p_cost['values'] = list(zip(params.CostPiecewisePoints[g], 
                                        (params.FuelCost[g]*val for val in params.CostPiecewiseValues[g])))
        g_d['p_cost'] =  p_cost

        ## NOTE: generators need unique names
        if keep_names:
            if g in gen_dict:
                raise RuntimeError("Nonunique generator names")
            gen_dict[g] = g_d
        else:
            gen_dict[g+'_t'] = g_d
        
    for g in params.AllNondispatchableGenerators:
        g_d = { 'generator_type':'renewable', }
        g_d['bus'] = renewable_gen_bus_dict[g]
        g_d['in_service'] = True
        g_d['fuel'] = params.NondispatchableGeneratorType[g]
        g_d['p_min'] = { 'data_type':'time_series', 
                            'values': [ params.MinNondispatchablePower[g,t] for t in params.TimePeriods ]
                       }
        g_d['p_max'] = { 'data_type':'time_series', 
                            'values': [ params.MaxNondispatchablePower[g,t] for t in params.TimePeriods ]
                       }
        ## NOTE: generators need unique names
        if keep_names:
            if g in gen_dict:
                raise RuntimeError("Nonunique generator names")
            gen_dict[g] = g_d
        else:
            gen_dict[g+'_r'] = g_d

    elements['generator'] = gen_dict

    storage_dict = {}
    for s in params.Storage:
        s_d = dict()
        s_d['bus'] = storage_bus_dict[s]
        s_d['min_discharge_rate'] = params.MinimumPowerOutputStorage[s]
        s_d['max_discharge_rate'] = params.MaximumPowerOutputStorage[s]
        s_d['min_charge_rate'] = params.MinimumPowerInputStorage[s]
        s_d['max_charge_rate'] = params.MaximumPowerInputStorage[s]
        s_d['ramp_up_output_60min'] = params.NominalRampUpLimitStorageOutput[s]
        s_d['ramp_down_output_60min'] = params.NominalRampDownLimitStorageOutput[s]
        s_d['ramp_up_input_60min'] = params.NominalRampUpLimitStorageInput[s]
        s_d['ramp_down_input_60min'] = params.NominalRampDownLimitStorageInput[s]
        s_d['energy_capacity'] = params.MaximumEnergyStorage[s]
        s_d['minimum_state_of_charge'] = params.MinimumSocStorage[s]
        s_d['charge_efficiency'] = params.InputEfficiencyEnergy[s]
        s_d['discharge_efficiency'] = params.OutputEfficiencyEnergy[s]
        s_d['retention_rate_60min'] = params.RetentionRate[s]
        s_d['initial_state_of_charge'] = params.StorageSocOnT0[s]
        s_d['initial_discharge_rate'] = params.StoragePowerOutputOnT0[s]
        s_d['initial_charge_rate'] = params.StoragePowerInputOnT0[s]

        storage_dict[s] = s_d
    elements['storage'] = storage_dict

    return md_dict


def _verify_must_run_t0_state_consistency(model):
    # ensure that the must-run flag and the t0 state are consistent. in partcular, make
    # sure that the unit has satisifed its minimum down time condition if UnitOnT0 is negative.
    
    def verify_must_run_t0_state_consistency_rule(m, g):
        if value(m.MustRun[g]):
            t0_state = value(m.UnitOnT0State[g])
            if t0_state < 0:
                min_down_time = value(m.MinimumDownTime[g])
                if abs(t0_state) < min_down_time:
                    print("DATA ERROR: The generator %s has been flagged as must-run, but its T0 state=%d is inconsistent with its minimum down time=%d" % (g, t0_state, min_down_time))
                    return False
        return True
    
    model.VerifyMustRunT0StateConsistency = BuildAction(model.ThermalGenerators, rule=verify_must_run_t0_state_consistency_rule)
    
def _populate_reserve_requirements(model):
    def populate_reserve_requirements_rule(m):
       reserve_factor = value(m.ReserveFactor)
       if reserve_factor > 0.0:
          for t in m.TimePeriods:
             demand = sum(value(m.Demand[b,t]) for b in sorted(m.Buses))
             m.ReserveRequirement[t] = (reserve_factor * demand)
    
    model.PopulateReserveRequirements = BuildAction(rule=populate_reserve_requirements_rule)

def setup_abstract_model(model):
    
    '''
    This adds an AbstractModel shell for dat file parsing
    '''
    warn_neg_load = False
    #
    # Parameters
    #
    
    ##############################################
    # string indentifiers for the set of busses. #
    ##############################################
    
    model.Buses = Set()
    
    ###################
    #   Load Zones    #
    ###################
    #Aggregated loads are distributed in the system based on load coefficient values
    
    model.Zones = Set(initialize=['SingleZone'])
    
    def buildBusZone(m):
        an_element = next(m.Zones.__iter__())
        if len(m.Zones) == 1 and an_element == 'SingleZone':
            for b in m.Buses:
                m.BusZone[b] = an_element
        else:
            print("Multiple buses is not supported by buildBusZone in ReferenceModel.py -- someone should fix that!")
            exit(1)
    
    model.BusZone = Param(model.Buses, mutable=True, within=Any)
    model.BuildBusZone = BuildAction(rule=buildBusZone)
    
    model.LoadCoefficient = Param(model.Buses, default=0.0)
    
    def total_load_coefficient_per_zone(m, z):
        return sum(m.LoadCoefficient[b] for b in m.Buses if str(value(m.BusZone[b]))==str(z))
    model.TotalLoadCoefficientPerZone = Param(model.Zones, initialize=total_load_coefficient_per_zone)
    
    def load_factors_per_bus(m,b):
        if (m.TotalLoadCoefficientPerZone[value(m.BusZone[b])] != 0.0):
            return m.LoadCoefficient[b]/m.TotalLoadCoefficientPerZone[value(m.BusZone[b])]
        else:
            return 0.0
    model.LoadFactor = Param(model.Buses, initialize=load_factors_per_bus, within=NonNegativeReals)

    model.BusKV = Param(model.Buses, default=1000.)
    
    ################################
    
    model.StageSet = Set(ordered=True) 
    
    # IMPORTANT: The stage set must be non-empty - otherwise, zero costs result.
    def check_stage_set(m):
       return (len(m.StageSet) != 0)
    model.CheckStageSet = BuildCheck(rule=check_stage_set)

    ## for backwards capatability (for now)
    model.TimePeriodLength = Param(default=1, within=PositiveReals)
    def time_period_length_validator(m):
        assert(m.TimePeriodLength == 1)
    model.TimePeriodLengthIsOne = BuildAction(rule=time_period_length_validator)
    
    ## IN HOURS, assert athat this must be a positive number
    model.TimePeriodLengthHours = Param(default=1.0, within=PositiveReals)

    ## in minutes, assert that this must be a positive integer
    model.TimePeriodLengthMinutes = Param(default=60, within=PositiveIntegers)

    ## sync the two time period lengths depending on the user's specification
    def harmonize_times(m):
        ## the user can only specify a non-default for 
        ## one of the time period lengths
        assert( (value(m.TimePeriodLengthHours) == 1.0) or (value(m.TimePeriodLengthMinutes) == 60) )
        if value(m.TimePeriodLengthHours) != 1.0:
            m.TimePeriodLengthMinutes = int(round(value(m.TimePeriodLengthHours)*60))
        if value(m.TimePeriodLengthMinutes) != 60:
            m.TimePeriodLengthHours = value(m.TimePeriodLengthMinutes)/60.
    
    model.HarmonizeTimes = BuildAction(rule=harmonize_times)

    model.NumTimePeriods = Param(within=PositiveIntegers, )
    
    model.InitialTime = Param(within=PositiveIntegers, default=1)
    model.TimePeriods = RangeSet(model.InitialTime, model.NumTimePeriods)
    
    # the following sets must must come from the data files or from an initialization function that uses 
    # a parameter that tells you when the stages end (and that thing needs to come from the data files)
    model.CommitmentTimeInStage = Set(model.StageSet, within=model.TimePeriods) 
    model.GenerationTimeInStage = Set(model.StageSet, within=model.TimePeriods)
    
    ##############################################
    # Network definition (S)
    ##############################################

    ## for older .dat files
    model.NumTransmissionLines = Param(default=0)
    def num_transimission_lines_validator(m):
        assert(m.NumTransmissionLines == 0)
    model.NumTransmissionLinesIsZero = BuildAction(rule=num_transimission_lines_validator)
    
    model.TransmissionLines = Set()
    
    model.BusFrom = Param(model.TransmissionLines, within=Any)
    model.BusTo   = Param(model.TransmissionLines, within=Any)

    model.Impedence = Param(model.TransmissionLines, within=NonNegativeReals)

    model.ThermalLimit = Param(model.TransmissionLines) # max flow across the line

    ## Interfaces
    ## NOTE: Lines in iterfaces should be all go "from" the
    ##       other network "to" the modeled network
    model.Interfaces = Set()

    model.InterfaceLines = Set(model.Interfaces, within=model.TransmissionLines)
    model.InterfaceFromLimit = Param(model.Interfaces, within=NonNegativeReals)
    model.InterfaceToLimit = Param(model.Interfaces, within=NonNegativeReals)
    
    ##########################################################
    # string indentifiers for the set of thermal generators. #
    # and their locations. (S)                               #
    ##########################################################
    
    model.ThermalGenerators = Set()
    model.ThermalGeneratorsAtBus = Set(model.Buses)
    
    # thermal generator types must be specified as 'N', 'C', 'G', and 'H',
    # with the obvious interpretation.
    # TBD - eventually add a validator.
    
    model.ThermalGeneratorType = Param(model.ThermalGenerators, within=Any, default='C')
    
    def verify_thermal_generator_buses_rule(m, g):
       for b in m.Buses:
          if g in m.ThermalGeneratorsAtBus[b]:
             return 
       print("DATA ERROR: No bus assigned for thermal generator=%s" % g)
       assert(False)
    
    model.VerifyThermalGeneratorBuses = BuildAction(model.ThermalGenerators, rule=verify_thermal_generator_buses_rule)
    
    model.QuickStart = Param(model.ThermalGenerators, within=Boolean, default=False)

    def init_quick_start_generators(m):
        return [g for g in m.ThermalGenerators if value(m.QuickStart[g]) == 1]
    model.QuickStartGenerators = Set(within=model.ThermalGenerators, initialize=init_quick_start_generators)
    
    # optionally force a unit to be on.
    model.MustRun = Param(model.ThermalGenerators, within=Boolean, default=False)
    
    def nd_gen_init(m,b):
        return []
    model.NondispatchableGeneratorsAtBus = Set(model.Buses, initialize=nd_gen_init)
    
    def nondispatchable_generator_init(m):
        for b in m.Buses:
            for gen in m.NondispatchableGeneratorsAtBus[b]:
                yield gen
    model.AllNondispatchableGenerators = Set(initialize=nondispatchable_generator_init)

    model.NondispatchableGeneratorType = Param(model.AllNondispatchableGenerators, within=Any, default='W')
    
    ######################
    #   Reserve Zones    #
    ######################
    
    # Generators are grouped in zones to provide zonal reserve requirements. #
    # All generators can contribute to global reserve requirements           #
    
    model.ReserveZones = Set()
    model.ZonalReserveRequirement = Param(model.ReserveZones, model.TimePeriods, default=0.0, within=NonNegativeReals)
    model.ReserveZoneLocation = Param(model.ThermalGenerators, default='None', within=Any)
    
    #################################################################
    # the global system demand, for each time period. units are MW. #
    # demand as at busses (S) so total demand is derived            #
    #################################################################
    
    # at the moment, we allow for negative demand. this is probably
    # not a good idea, as "stuff" representing negative demand - including
    # renewables, interchange schedules, etc. - should probably be modeled
    # explicitly.
    
    # Demand can also be given by Zones
    
    model.DemandPerZone = Param(model.Zones, model.TimePeriods, default=0.0, )
    
    # Convert demand by zone to demand by bus
    def demand_per_bus_from_demand_per_zone(m,b,t):
        return m.DemandPerZone[value(m.BusZone[b]), t] * m.LoadFactor[b]
    model.Demand = Param(model.Buses, model.TimePeriods, initialize=demand_per_bus_from_demand_per_zone, )
    
    def calculate_total_demand(m, t):
        return sum(value(m.Demand[b,t]) for b in sorted(m.Buses))
    model.TotalDemand = Param(model.TimePeriods, initialize=calculate_total_demand)
    
    # at this point, a user probably wants to see if they have negative demand.
    def warn_about_negative_demand_rule(m, b, t):
       this_demand = value(m.Demand[b,t])
       if this_demand < 0.0:
          print("***WARNING: The demand at bus="+str(b)+" for time period="+str(t)+" is negative - value="+str(this_demand)+"; model="+str(m.name)+".")
    
    if warn_neg_load:
        model.WarnAboutNegativeDemand = BuildAction(model.Buses, model.TimePeriods, rule=warn_about_negative_demand_rule)
    
    ##################################################################
    # the global system reserve, for each time period. units are MW. #
    # NOTE: We don't have per-bus / zonal reserve requirements. they #
    #       would be easy to add. (dlw oct 2013: this comment is incorrect, I think)                                   #
    ##################################################################
    
    # we provide two mechanisms to specify reserve requirements. the
    # first is a scaling factor relative to demand, on a per time 
    # period basis. the second is an explicit parameter that specifies
    # the reserver requirement on a per-time-period basis. if the 
    # reserve requirement factor is > 0, then it is used to populate
    # the reserve requirements. otherwise, the user-supplied reserve
    # requirements are used.
    
    model.ReserveFactor = Param(within=Reals, default=-1.0, )
    
    model.ReserveRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0, mutable=True )
    
    _populate_reserve_requirements(model)
    
    ##############################################################
    # failure probability for each generator, in any given hour. #
    # not used within the model itself at present, but rather    #
    # used by scripts that read / manipulate the model.          #
    ##############################################################
    
    def probability_failure_validator(m, v, g):
       return v >= 0.0 and v <= 1.0
    
    model.FailureProbability = Param(model.ThermalGenerators, validate=probability_failure_validator, default=0.0)
    
    #####################################################################################
    # a binary indicator as to whether or not each generator is on-line during a given  #
    # time period. intended to represent a sampled realization of the generator failure #
    # probability distributions. strictly speaking, we interpret this parameter value   #
    # as indicating whether or not the generator is contributing (injecting) power to   #
    # the PowerBalance constraint. this parameter is not intended to be used in the     #
    # context of ramping or time up/down constraints.                                   # 
    #####################################################################################
    
    model.GeneratorForcedOutage = Param(model.ThermalGenerators * model.TimePeriods, within=Binary, default=False)
    
    ####################################################################################
    # minimum and maximum generation levels, for each thermal generator. units are MW. #
    # could easily be specified on a per-time period basis, but are not currently.     #
    ####################################################################################
    
    # you can enter generator limits either once for the generator or for each period (or just take 0)
    
    model.MinimumPowerOutput = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0)
    
    def maximum_power_output_validator(m, v, g):
       return v >= value(m.MinimumPowerOutput[g])
    
    model.MaximumPowerOutput = Param(model.ThermalGenerators, within=NonNegativeReals, validate=maximum_power_output_validator, default=0.0)
    
    # wind is similar, but max and min will be equal for non-dispatchable wind
    
    model.MinNondispatchablePower = Param(model.AllNondispatchableGenerators, model.TimePeriods, within=NonNegativeReals, default=0.0, )
    
    def maximum_nd_output_validator(m, v, g, t):
       return v >= value(m.MinNondispatchablePower[g,t])
    
    model.MaxNondispatchablePower = Param(model.AllNondispatchableGenerators, model.TimePeriods, within=NonNegativeReals, default=0.0, validate=maximum_nd_output_validator)
    
    #################################################
    # generator ramp up/down rates. units are MW/h. #
    # IMPORTANT: Generator ramp limits can exceed   #
    # the maximum power output, because it is the   #
    # ramp limit over an hour. If the unit can      #
    # fully ramp in less than an hour, then this    #
    # will occur.                                   #
    #################################################
    
    # limits for normal time periods
    model.NominalRampUpLimit = Param(model.ThermalGenerators, within=NonNegativeReals, )
    model.NominalRampDownLimit = Param(model.ThermalGenerators, within=NonNegativeReals, )
    
    # limits for time periods in which generators are brought on or off-line.
    # must be no less than the generator minimum output.
    def ramp_limit_validator(m, v, g):
       return v >= m.MinimumPowerOutput[g]

    ## These defaults follow what is in most market manuals
    ## We scale this for the time period below
    def startup_ramp_default(m, g):
        return m.MinimumPowerOutput[g]+m.NominalRampUpLimit[g]/2.
    def shutdown_ramp_default(m, g):
        return m.MinimumPowerOutput[g]+m.NominalRampDownLimit[g]/2.

    model.StartupRampLimit = Param(model.ThermalGenerators, within=NonNegativeReals, default=startup_ramp_default, validate=ramp_limit_validator, )
    model.ShutdownRampLimit = Param(model.ThermalGenerators, within=NonNegativeReals, default=shutdown_ramp_default,  validate=ramp_limit_validator, )
    
    ##########################################################################################################
    # the minimum number of time periods that a generator must be on-line (off-line) once brought up (down). #
    ##########################################################################################################
    
    model.MinimumUpTime = Param(model.ThermalGenerators, within=NonNegativeIntegers, default=0)
    model.MinimumDownTime = Param(model.ThermalGenerators, within=NonNegativeIntegers, default=0)
    
    
    #############################################
    # unit on state at t=0 (initial condition). #
    #############################################
    
    # if positive, the number of hours prior to (and including) t=0 that the unit has been on.
    # if negative, the number of hours prior to (and including) t=0 that the unit has been off.
    # the value cannot be 0, by definition.
    
    def t0_state_nonzero_validator(m, v, g):
        return v != 0
    
    model.UnitOnT0State = Param(model.ThermalGenerators, within=Reals, validate=t0_state_nonzero_validator, )
    
    def t0_unit_on_rule(m, g):
        return int(value(m.UnitOnT0State[g]) >= 1)
    
    model.UnitOnT0 = Param(model.ThermalGenerators, within=Binary, initialize=t0_unit_on_rule, )
    
    _verify_must_run_t0_state_consistency(model)

    ####################################################################
    # generator power output at t=0 (initial condition). units are MW. #
    ####################################################################
    
    def between_limits_validator(m, v, g):
       status = (v <= (value(m.MaximumPowerOutput[g]) * value(m.UnitOnT0[g]))  and v >= (value(m.MinimumPowerOutput[g]) * value(m.UnitOnT0[g])))
       if status == False:
          print("Failed to validate PowerGeneratedT0 value for g="+g+"; new value="+str(v)+", UnitOnT0="+str(value(m.UnitOnT0[g])))
       return v <= (value(m.MaximumPowerOutput[g]) * value(m.UnitOnT0[g]))  and v >= (value(m.MinimumPowerOutput[g]) * value(m.UnitOnT0[g]))
    model.PowerGeneratedT0 = Param(model.ThermalGenerators, within=NonNegativeReals, validate=between_limits_validator, )
    
    
    ###############################################
    # startup cost parameters for each generator. #
    ###############################################
    
    # startup costs are conceptually expressed as pairs (x, y), where x represents the number of hours that a unit has been off and y represents
    # the cost associated with starting up the unit after being off for x hours. these are broken into two distinct ordered sets, as follows.
    
    def startup_lags_init_rule(m, g):
       return [value(m.MinimumDownTime[g])] 
    model.StartupLags = Set(model.ThermalGenerators, within=NonNegativeIntegers, ordered=True, initialize=startup_lags_init_rule) # units are hours / time periods.
    
    def startup_costs_init_rule(m, g):
       return [0.0] 
    
    model.StartupCosts = Set(model.ThermalGenerators, within=NonNegativeReals, ordered=True, initialize=startup_costs_init_rule) # units are $.
    
    # startup lags must be monotonically increasing...
    def validate_startup_lags_rule(m, g):
       startup_lags = list(m.StartupLags[g])
    
       if len(startup_lags) == 0:
          print("DATA ERROR: The number of startup lags for thermal generator="+str(g)+" must be >= 1.")
          assert(False)
    
       if startup_lags[0] != value(m.MinimumDownTime[g]):
          print("DATA ERROR: The first startup lag for thermal generator="+str(g)+" must be equal the minimum down time="+str(value(m.MinimumDownTime[g]))+".")
          assert(False)      
    
       for i in range(0, len(startup_lags)-1):
          if startup_lags[i] >= startup_lags[i+1]:
             print("DATA ERROR: Startup lags for thermal generator="+str(g)+" must be monotonically increasing.")
             assert(False)
    
    model.ValidateStartupLags = BuildAction(model.ThermalGenerators, rule=validate_startup_lags_rule)
    
    # while startup costs must be monotonically non-decreasing!
    def validate_startup_costs_rule(m, g):
       startup_costs = list(m.StartupCosts[g])
       for i in range(0, len(startup_costs)-2):
          if startup_costs[i] > startup_costs[i+1]:
             print("DATA ERROR: Startup costs for thermal generator="+str(g)+" must be monotonically non-decreasing.")
             assert(False)
    
    model.ValidateStartupCosts = BuildAction(model.ThermalGenerators, rule=validate_startup_costs_rule)
    
    def validate_startup_lag_cost_cardinalities(m, g):
       if len(m.StartupLags[g]) != len(m.StartupCosts[g]):
          print("DATA ERROR: The number of startup lag entries ("+str(len(m.StartupLags[g]))+") for thermal generator="+str(g)+" must equal the number of startup cost entries ("+str(len(m.StartupCosts[g]))+")")
          assert(False)
    
    model.ValidateStartupLagCostCardinalities = BuildAction(model.ThermalGenerators, rule=validate_startup_lag_cost_cardinalities)
    
    # for purposes of defining constraints, it is useful to have a set to index the various startup costs parameters.
    # entries are 1-based indices, because they are used as indicies into Pyomo sets - which use 1-based indexing.
    
    def startup_cost_indices_init_rule(m, g):
       return range(1, len(m.StartupLags[g])+1)
    
    model.StartupCostIndices = Set(model.ThermalGenerators, within=NonNegativeIntegers, initialize=startup_cost_indices_init_rule)

    ##################################################################################
    # shutdown cost for each generator. in the literature, these are often set to 0. #
    ##################################################################################
    
    model.ShutdownFixedCost = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0) # units are $.
    
    ## BEGIN PRODUCTION COST
    ## NOTE: For better or worse, we handle scaling this to the time period length in the objective function.
    ##       In particular, this is done in objective.py.

    ##################################################################################################################
    # production cost coefficients (for the quadratic) a0=constant, a1=linear coefficient, a2=quadratic coefficient. #
    ##################################################################################################################
    
    model.ProductionCostA0 = Param(model.ThermalGenerators, default=0.0) # units are $/hr (or whatever the time unit is).
    model.ProductionCostA1 = Param(model.ThermalGenerators, default=0.0) # units are $/MWhr.
    model.ProductionCostA2 = Param(model.ThermalGenerators, default=0.0) # units are $/(MWhr^2).
    
    # the parameters below are populated if cost curves are specified as linearized heat rate increment segments.
    #
    # CostPiecewisePoints represents the power output levels defining the segment boundaries.
    # these *must* include the minimum and maximum power output points - a validation check
    # if performed below.
    # 
    # CostPiecewiseValues are the absolute heat rates / costs associated with the corresponding 
    # power output levels. the precise interpretation of whether a value is a heat rate or a cost
    # depends on the value of the FuelCost parameter, specified below.
    
    
    def piecewise_type_validator(m, v):
       return (v == "NoPiecewise") or (v == "Absolute")
    
    def piecewise_type_init(m):
        boo = False
        for g in m.ThermalGenerators:
            if not (m.ProductionCostA0[g] == 0.0 and m.ProductionCostA1[g] == 0.0 and m.ProductionCostA2[g] == 0.0):
                boo = True
                break
        if boo:
            return "NoPiecewise"
        else:
            return "Absolute"
    
    model.PiecewiseType = Param(validate=piecewise_type_validator,initialize=piecewise_type_init, within=Any)
    
    def piecewise_init(m, g):
        return []
    
    model.CostPiecewisePoints = Set(model.ThermalGenerators, initialize=piecewise_init, ordered=True, within=NonNegativeReals)
    model.CostPiecewiseValues = Set(model.ThermalGenerators, initialize=piecewise_init, ordered=True, within=NonNegativeReals)
    
    # a check to ensure that the cost piecewise point parameter was correctly populated.
    # these are global checks, which cannot be performed as a set validation (which 
    # operates on a single element at a time).
    
    def validate_cost_piecewise_points_and_values_rule(m, g):
        
        if value(m.PiecewiseType) == "NoPiecewise":
            # if there isn't any piecewise data specified, we shouldn't find any.
            if len(m.CostPiecewisePoints[g]) > 0:
                print("DATA ERROR: The PiecewiseType parameter was set to NoPiecewise, but piecewise point data was specified!")
                return False
            # if there isn't anything to validate and we didn't expect piecewise 
            # points, we can safely skip the remaining validation steps.
            return True
        else:
            # if the user said there was going to be piecewise data and none was 
            # supplied, they should be notified as to such.
            if len(m.CostPiecewisePoints[g]) == 0:
                print("DATA ERROR: The PiecewiseType parameter was set to something other than NoPiecewise, but no piecewise point data was specified!")
                return False
    
       # per the requirement below, there have to be at least two piecewise points if there are any.
    
        min_output = value(m.MinimumPowerOutput[g])
        max_output = value(m.MaximumPowerOutput[g])   
    
        points = list(m.CostPiecewisePoints[g])
    
        if min_output not in points:
            print("DATA ERROR: Cost piecewise points for generator g="+str(g)+" must contain the minimum output level="+str(min_output))
            return False
    
        if max_output not in points:
            print("DATA ERROR: Cost piecewise points for generator g="+str(g)+" must contain the maximum output level="+str(max_output))
            return False

        return True
    
    model.ValidateCostPiecewisePointsAndValues = BuildCheck(model.ThermalGenerators, rule=validate_cost_piecewise_points_and_values_rule)
    
    model.FuelCost = Param(model.ThermalGenerators, default=1.0, within=Reals) 
    
    # Minimum production cost (needed because Piecewise constraint on ProductionCost 
    # has to have lower bound of 0, so the unit can cost 0 when off -- this is added
    # back in to the objective if a unit is on
    def minimum_production_cost(m, g):
        if len(m.CostPiecewisePoints[g]) > 1:
            return m.CostPiecewiseValues[g].first() * m.FuelCost[g]
        else:
            return  m.FuelCost[g] * \
                   (m.ProductionCostA0[g] + \
                    m.ProductionCostA1[g] * m.MinimumPowerOutput[g] + \
                    m.ProductionCostA2[g] * (m.MinimumPowerOutput[g]**2))
    model.MinimumProductionCost = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=minimum_production_cost, )
    
    ##############################################################################################
    # number of pieces in the linearization of each generator's quadratic cost production curve. #
    ##############################################################################################
    
    model.NumGeneratorCostCurvePieces = Param(within=PositiveIntegers, default=2, )

    ModeratelyBigPenalty = 1e3
    
    model.ReserveShortfallPenalty = Param(within=NonNegativeReals, default=ModeratelyBigPenalty, )

    #########################################
    # penalty costs for constraint violation #
    #########################################
    
    BigPenalty = 1e4
    
    model.LoadMismatchPenalty = Param(within=NonNegativeReals, default=BigPenalty, )

    ## END PRODUCTION COST CALCULATIONS

    #
    # STORAGE parameters
    #
    
    
    model.Storage = Set()
    model.StorageAtBus = Set(model.Buses)
    
    def verify_storage_buses_rule(m, s):
        for b in m.Buses:
            if s in m.StorageAtBus[b]:
                return
        print("DATA ERROR: No bus assigned for storage element=%s" % s)
        assert(False)
    
    model.VerifyStorageBuses = BuildAction(model.Storage, rule=verify_storage_buses_rule)
    
    ####################################################################################
    # minimum and maximum power ratings, for each storage unit. units are MW.          #
    # could easily be specified on a per-time period basis, but are not currently.     #
    ####################################################################################
    
    # Storage power output >0 when discharging
    
    model.MinimumPowerOutputStorage = Param(model.Storage, within=NonNegativeReals, default=0.0)
    
    def maximum_power_output_validator_storage(m, v, s):
        return v >= value(m.MinimumPowerOutputStorage[s])
    
    model.MaximumPowerOutputStorage = Param(model.Storage, within=NonNegativeReals, validate=maximum_power_output_validator_storage, default=0.0)
    
    #Storage power input >0 when charging
    
    model.MinimumPowerInputStorage = Param(model.Storage, within=NonNegativeReals, default=0.0)
    
    def maximum_power_input_validator_storage(m, v, s):
        return v >= value(m.MinimumPowerInputStorage[s])
    
    model.MaximumPowerInputStorage = Param(model.Storage, within=NonNegativeReals, validate=maximum_power_input_validator_storage, default=0.0)
    
    ###############################################
    # storage ramp up/down rates. units are MW/h. #
    ###############################################
    
    # ramp rate limits when discharging
    model.NominalRampUpLimitStorageOutput    = Param(model.Storage, within=NonNegativeReals)
    model.NominalRampDownLimitStorageOutput  = Param(model.Storage, within=NonNegativeReals)
    
    # ramp rate limits when charging
    model.NominalRampUpLimitStorageInput     = Param(model.Storage, within=NonNegativeReals)
    model.NominalRampDownLimitStorageInput   = Param(model.Storage, within=NonNegativeReals)
    
    ####################################################################################
    # minimum state of charge (SOC) and maximum energy ratings, for each storage unit. #
    # units are MWh for energy rating and p.u. (i.e. [0,1]) for SOC     #
    ####################################################################################
    
    # you enter storage energy ratings once for each storage unit
    
    model.MaximumEnergyStorage = Param(model.Storage, within=NonNegativeReals, default=0.0)
    model.MinimumSocStorage = Param(model.Storage, within=PercentFraction, default=0.0)
    
    ################################################################################
    # round trip efficiency for each storage unit given as a fraction (i.e. [0,1]) #
    ################################################################################
    
    model.InputEfficiencyEnergy  = Param(model.Storage, within=PercentFraction, default=1.0)
    model.OutputEfficiencyEnergy = Param(model.Storage, within=PercentFraction, default=1.0)
    model.RetentionRate          = Param(model.Storage, within=PercentFraction, default=1.0) ## assumed to be %/hr

    
    ########################################################################
    # end-point SOC for each storage unit. units are in p.u. (i.e. [0,1])  #
    ########################################################################
    
    # end-point values are the SOC targets at the final time period. With no end-point constraints
    # storage units will always be empty at the final time period.
    
    model.EndPointSocStorage = Param(model.Storage, within=PercentFraction, default=0.5)
    
    ############################################################
    # storage initial conditions: SOC, power output and input  #
    ############################################################
    
    def t0_storage_power_input_validator(m, v, s):
        return (v >= value(m.MinimumPowerInputStorage[s])) and (v <= value(m.MaximumPowerInputStorage[s]))
    
    def t0_storage_power_output_validator(m, v, s):
        return (v >= value(m.MinimumPowerInputStorage[s])) and (v <= value(m.MaximumPowerInputStorage[s]))
    
    model.StoragePowerOutputOnT0 = Param(model.Storage, within=NonNegativeReals, validate=t0_storage_power_output_validator, default=0.0)
    model.StoragePowerInputOnT0  = Param(model.Storage, within=NonNegativeReals, validate=t0_storage_power_input_validator, default=0.0)
    model.StorageSocOnT0         = Param(model.Storage, within=PercentFraction, default=0.5)

    return model
