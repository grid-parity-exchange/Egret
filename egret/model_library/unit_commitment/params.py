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
import math
from egret.data.model_data import map_items, zip_items
from egret.model_library.transmission import tx_utils
    
from .uc_utils import add_model_attr, build_uc_time_mapping

component_name = 'data_loader'

def _verify_must_run_t0_state_consistency(model):
    # ensure that the must-run flag and the t0 state are consistent. in partcular, make
    # sure that the unit has satisifed its minimum down time condition if UnitOnT0 is negative.
    
    def verify_must_run_t0_state_consistency_rule(m, g):
        t0_state = value(m.UnitOnT0State[g]) / value(m.TimePeriodLengthHours)
        if t0_state < 0:
            min_down_time = value(m.ScaledMinimumDownTime[g])
            if abs(t0_state) < min_down_time:
                for t in range(m.TimePeriods.first(), value(m.InitialTimePeriodsOffLine[g])+m.TimePeriods.first()):
                    fixed_commitment = value(m.FixedCommitment[g,t])
                    if (fixed_commitment is not None) and (fixed_commitment == 1):
                        print("DATA ERROR: The generator %s has been flagged as must-run at time %d, but its T0 state=%d is inconsistent with its minimum down time=%d" % (g, t, t0_state, min_down_time))
                        return False
        else: # t0_state > 0
            min_up_time = value(m.ScaledMinimumUpTime[g])
            if abs(t0_state) < min_up_time:
                for t in range(m.TimePeriods.first(), value(m.InitialTimePeriodsOnLine[g])+m.TimePeriods.first()):
                    fixed_commitment = value(m.FixedCommitment[g,t])
                    if (fixed_commitment is not None) and (fixed_commitment == 0):
                        print("DATA ERROR: The generator %s has been flagged as off at time %d, but its T0 state=%d is inconsistent with its minimum up time=%d" % (g, t, t0_state, min_down_time))
                        return False
        return True
    
    model.VerifyMustRunT0StateConsistency = BuildAction(model.ThermalGenerators, rule=verify_must_run_t0_state_consistency_rule)
    

def _add_initial_time_periods_on_off_line(model):
    #######################################################################################
    # the number of time periods that a generator must initally on-line (off-line) due to #
    # its minimum up time (down time) constraint.                                         #
    #######################################################################################
    
    def initial_time_periods_online_rule(m, g):
       if not value(m.UnitOnT0[g]):
          return 0
       else:
          return int(min(value(m.NumTimePeriods),
                 round(max(0, value(m.MinimumUpTime[g]) - value(m.UnitOnT0State[g])) / value(m.TimePeriodLengthHours))))
    
    model.InitialTimePeriodsOnLine = Param(model.ThermalGenerators, within=NonNegativeIntegers, initialize=initial_time_periods_online_rule, mutable=True)
    
    def initial_time_periods_offline_rule(m, g):
       if value(m.UnitOnT0[g]):
          return 0
       else:
          return int(min(value(m.NumTimePeriods),
                 round(max(0, value(m.MinimumDownTime[g]) + value(m.UnitOnT0State[g])) / value(m.TimePeriodLengthHours)))) # m.UnitOnT0State is negative if unit is off
    
    model.InitialTimePeriodsOffLine = Param(model.ThermalGenerators, within=NonNegativeIntegers, initialize=initial_time_periods_offline_rule, mutable=True)

@add_model_attr(component_name)
def load_params(model, model_data):
    
    '''
    This loads unit commitment params from a GridModel object
    '''
    warn_neg_load = False

    md = model_data
    model.model_data = model_data

    system = md.data['system']
    elements = md.data['elements']

    time_keys = system['time_indices']
    TimeMapper = build_uc_time_mapping(time_keys)

    
    ## insert potentially missing keys
    if 'branch' not in elements:
        elements['branch'] = dict()
    if 'interface' not in elements:
        elements['interface'] = dict()
    if 'storage' not in elements:
        elements['storage'] = dict()

    ## NOTE: generator, bus, and load should be in here for a well-defined problem

    loads = dict(md.elements(element_type='load'))
    thermal_gens = dict(md.elements(element_type='generator', generator_type='thermal'))
    renewable_gens = dict(md.elements(element_type='generator', generator_type='renewable'))
    buses = dict(md.elements(element_type='bus'))
    shunts = dict(md.elements(element_type='shunt'))
    branches = dict(md.elements(element_type='branch'))
    storage = dict(md.elements(element_type='storage'))

    thermal_gen_attrs = md.attributes(element_type='generator', generator_type='thermal')
    renewable_gen_attrs = md.attributes(element_type='generator', generator_type='renewable')
    bus_attrs = md.attributes(element_type='bus')
    branch_attrs = md.attributes(element_type='branch')
    load_attrs = md.attributes(element_type='load')
    interface_attrs = md.attributes(element_type='interface')
    storage_attrs = md.attributes(element_type='storage')


    inlet_branches_by_bus, outlet_branches_by_bus = \
        tx_utils.inlet_outlet_branches_by_bus(branches, buses)
    thermal_gens_by_bus = tx_utils.gens_by_bus(buses, thermal_gens)
    renewable_gens_by_bus = tx_utils.gens_by_bus(buses, renewable_gens)
    storage_by_bus = tx_utils.gens_by_bus(buses, storage)

    ### get the fixed shunts at the buses
    bus_bs_fixed_shunts, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    ## attach some of these to the model object for ease/speed later
    #model._loads = loads
    model._buses = buses
    model._branches = branches
    model._shunts = shunts
    model._bus_gs_fixed_shunts = bus_gs_fixed_shunts
    #model._bus_bs_fixed_shunts = bus_bs_fixed_shunts
    #model._TimeMapper = TimeMapper

    #
    # Parameters
    #
    
    ##############################################
    # string indentifiers for the set of busses. #
    ##############################################
    
    model.Buses = Set(initialize=bus_attrs['names'])
    
    if 'reference_bus' in system and system['reference_bus'] in model.Buses:
        reference_bus = system['reference_bus']
    else:
        reference_bus = list(sorted(m.Buses))[0]

    model.ReferenceBus = Param(within=model.Buses, initialize=reference_bus)

    if 'reference_bus_angle' in system:
        ref_angle = system['reference_bus_angle']
    else:
        ref_angle = 0.

    model.ReferenceBusAngle = Param(within=Reals, initialize=ref_angle)
    
    ################################
    
    ## in minutes, assert that this must be a positive integer
    model.TimePeriodLengthMinutes = Param(default=60, within=PositiveIntegers, initialize=system['time_period_length_minutes'])

    ## IN HOURS, assert athat this must be a positive number
    model.TimePeriodLengthHours = Param(default=value(model.TimePeriodLengthMinutes)/60., within=PositiveReals)

    model.NumTimePeriods = Param(within=PositiveIntegers, initialize=len(system['time_indices']))
    
    model.InitialTime = Param(within=PositiveIntegers, default=1)
    model.TimePeriods = RangeSet(model.InitialTime, model.NumTimePeriods)
    
    ## For now, hard code these. Probably need to be able to specify in model_data
    model.StageSet = Set(ordered=True, initialize=['Stage_1', 'Stage_2']) 

    # the following sets must must come from the data files or from an initialization function that uses 
    # a parameter that tells you when the stages end (and that thing needs to come from the data files)

    model.CommitmentTimeInStage = Set(model.StageSet, within=model.TimePeriods,
                                        initialize={'Stage_1':model.TimePeriods, 'Stage_2': list() } ) 
    model.GenerationTimeInStage = Set(model.StageSet, within=model.TimePeriods,
                                        initialize={'Stage_1': list(), 'Stage_2': model.TimePeriods } )
    
    ##############################################
    # Network definition (S)
    ##############################################
    
    model.TransmissionLines = Set(initialize=branch_attrs['names'])
    
    model.BusFrom = Param(model.TransmissionLines, initialize=branch_attrs.get('from_bus'))
    model.BusTo   = Param(model.TransmissionLines, initialize=branch_attrs.get('to_bus'))

    model.LinesTo = Set(model.Buses, initialize=inlet_branches_by_bus)
    model.LinesFrom = Set(model.Buses, initialize=outlet_branches_by_bus)

    model.Impedence = Param(model.TransmissionLines, within=NonNegativeReals, initialize=branch_attrs.get('reactance'))

    model.ThermalLimit = Param(model.TransmissionLines, initialize=branch_attrs.get('rating_long_term')) # max flow across the line

    model.LineOutOfService = Param(model.TransmissionLines, model.TimePeriods, within=Boolean, default=False,
                                    initialize=TimeMapper(branch_attrs.get('planned_outage')))

    ## Interfaces
    ## NOTE: Lines in iterfaces should be all go "from" the
    ##       other network "to" the modeled network
    model.Interfaces = Set(initialize=interface_attrs['names'])

    model.InterfaceLines = Set(model.Interfaces, within=model.TransmissionLines, initialize=interface_attrs.get('lines'))
    model.InterfaceFromLimit = Param(model.Interfaces, within=NonNegativeReals, initialize=interface_attrs.get('interface_from_limit'))
    model.InterfaceToLimit = Param(model.Interfaces, within=NonNegativeReals, initialize=interface_attrs.get('interface_to_limit'))
    
    ##########################################################
    # string indentifiers for the set of thermal generators. #
    # and their locations. (S)                               #
    ##########################################################
    
    model.ThermalGenerators = Set(initialize=thermal_gen_attrs['names'])
    model.ThermalGeneratorsAtBus = Set(model.Buses, initialize=thermal_gens_by_bus)
    
    model.ThermalGeneratorType = Param(model.ThermalGenerators, within=Any, default='C', initialize=thermal_gen_attrs.get('fuel'))
    
    def verify_thermal_generator_buses_rule(m, g):
       for b in m.Buses:
          if g in m.ThermalGeneratorsAtBus[b]:
             return 
       print("DATA ERROR: No bus assigned for thermal generator=%s" % g)
       assert(False)
    
    model.VerifyThermalGeneratorBuses = BuildAction(model.ThermalGenerators, rule=verify_thermal_generator_buses_rule)
    
    model.QuickStart = Param(model.ThermalGenerators, within=Boolean, default=False, initialize=thermal_gen_attrs.get('quickstart_capable'))
    
    def init_quick_start_generators(m):
        return [g for g in m.ThermalGenerators if value(m.QuickStart[g]) == 1]
    
    model.QuickStartGenerators = Set(within=model.ThermalGenerators, initialize=init_quick_start_generators)
    
    # optionally force a unit to be on/off
    model.FixedCommitmentTypes = Set(initialize=[0,1,None])
    model.FixedCommitment = Param(model.ThermalGenerators,
                                  model.TimePeriods,
                                  within=model.FixedCommitmentTypes,
                                  default=None,
                                  initialize=TimeMapper(thermal_gen_attrs.get('fixed_commitment')),)
    
    model.NondispatchableGeneratorsAtBus = Set(model.Buses, initialize=renewable_gens_by_bus)
    
    model.AllNondispatchableGenerators = Set(initialize=renewable_gen_attrs['names'])

    model.NondispatchableGeneratorType = Param(model.AllNondispatchableGenerators, within=Any, default='W', 
                                                initialize=renewable_gen_attrs.get('fuel'))
    
    
    #################################################################
    # the global system demand, for each time period. units are MW. #
    # demand as at busses (S) so total demand is derived            #
    #################################################################
    
    # at the moment, we allow for negative demand. this is probably
    # not a good idea, as "stuff" representing negative demand - including
    # renewables, interchange schedules, etc. - should probably be modeled
    # explicitly.

    bus_loads = { (b,t) : 0 for b in bus_attrs['names'] for t in model.TimePeriods}

    for lname, load in loads.items():
        bus = load['bus']
        load_time = TimeMapper(load['p_load'])
        load_in_service = TimeMapper(load['in_service'])
        for t in model.TimePeriods:
            bus_loads[bus, t] += load_in_service(None,t)*load_time(None,t)
    model.Demand = Param(model.Buses, model.TimePeriods, initialize=bus_loads, mutable=True)
    
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
    ##################################################################

    reserve_requirement = system.get("reserve_requirement", 0.)
    model.ReserveRequirement = Param(model.TimePeriods, within=NonNegativeReals, 
                                        initialize=TimeMapper(reserve_requirement), mutable=True)
    
    
    
    ####################################################################################
    # minimum and maximum generation levels, for each thermal generator. units are MW. #
    # could easily be specified on a per-time period basis, but are not currently.     #
    ####################################################################################
    
    # you can enter generator limits either once for the generator or for each period (or just take 0)
    
    model.MinimumPowerOutput = Param(model.ThermalGenerators, within=NonNegativeReals, 
                                        initialize=thermal_gen_attrs['p_min'],
                                        default=0.0)
    
    def maximum_power_output_validator(m, v, g):
       return v >= value(m.MinimumPowerOutput[g])
    
    model.MaximumPowerOutput = Param(model.ThermalGenerators, within=NonNegativeReals, 
                                        validate=maximum_power_output_validator, 
                                        initialize=thermal_gen_attrs['p_max'],
                                        default=0.0)

    model.MinimumReactivePowerOutput = Param(model.ThermalGenerators, within=Reals,
                                                initialize=thermal_gen_attrs.get('q_min'),
                                                default=0.0)

    def maximum_reactive_output_validator(m, v, g):
        return v >= value(m.MinimumReactivePowerOutput[g])

    model.MaximumReactivePowerOutput = Param(model.ThermalGenerators, within=Reals,
                                                initialize=thermal_gen_attrs.get('q_max'),
                                                default=0.0)
    
    # wind is similar, but max and min will be equal for non-dispatchable wind
    
    model.MinNondispatchablePower = Param(model.AllNondispatchableGenerators,
                                            model.TimePeriods, 
                                            within=NonNegativeReals,
                                            default=0.0,
                                            mutable=True,
                                            initialize=TimeMapper(renewable_gen_attrs.get('p_min')))
    
    def maximum_nd_output_validator(m, v, g, t):
       return v >= value(m.MinNondispatchablePower[g,t])
    
    model.MaxNondispatchablePower = Param(model.AllNondispatchableGenerators,
                                            model.TimePeriods,
                                            within=NonNegativeReals,
                                            default=0.0,
                                            mutable=True,
                                            validate=maximum_nd_output_validator,
                                            initialize=TimeMapper(renewable_gen_attrs.get('p_max')))
    
    #################################################
    # generator ramp up/down rates. units are MW/h. #
    # IMPORTANT: Generator ramp limits can exceed   #
    # the maximum power output, because it is the   #
    # ramp limit over an hour. If the unit can      #
    # fully ramp in less than an hour, then this    #
    # will occur.                                   #
    #################################################
    
    # limits for normal time periods
    model.NominalRampUpLimit = Param(model.ThermalGenerators,
                                        within=NonNegativeReals,
                                        mutable=True,
                                        initialize=thermal_gen_attrs['ramp_up_60min'])
    model.NominalRampDownLimit = Param(model.ThermalGenerators,
                                        within=NonNegativeReals,
                                        mutable=True,
                                        initialize=thermal_gen_attrs['ramp_down_60min'])
    
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

    model.StartupRampLimit = Param(model.ThermalGenerators, 
                                    within=NonNegativeReals,
                                    default=startup_ramp_default,
                                    validate=ramp_limit_validator,
                                    mutable=True,
                                    initialize=thermal_gen_attrs.get('startup_capacity'))
    model.ShutdownRampLimit = Param(model.ThermalGenerators, 
                                    within=NonNegativeReals,
                                    default=shutdown_ramp_default, 
                                    validate=ramp_limit_validator,
                                    mutable=True,
                                    initialize=thermal_gen_attrs.get('shutdown_capacity'))
    
    ## These get used in the basic UC constraints, which implicity assume RU, RD <= Pmax
    def scale_ramp_up(m, g):
        temp = m.NominalRampUpLimit[g] * m.TimePeriodLengthHours
        if value(temp) > value(m.MaximumPowerOutput[g]):
            return m.MaximumPowerOutput[g]
        else:
            return temp
    model.ScaledNominalRampUpLimit = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=scale_ramp_up, mutable=True)
    
    def scale_ramp_down(m, g):
        temp = m.NominalRampDownLimit[g] * m.TimePeriodLengthHours
        if value(temp) > value(m.MaximumPowerOutput[g]):
            return m.MaximumPowerOutput[g]
        else:
            return temp
    model.ScaledNominalRampDownLimit = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=scale_ramp_down, mutable=True)
    
    def scale_startup_limit(m, g):
        ## temp now has the "running room" over Pmin. This will be scaled for the time period length, 
        ## most market models do not have this notion, so this is set-up so that the defaults
        ## will be scaled as they would be in most market models
        temp = (m.StartupRampLimit[g] - m.MinimumPowerOutput[g])*m.TimePeriodLengthHours
        if value(temp) > value(m.MaximumPowerOutput[g] - m.MinimumPowerOutput[g]):
            return m.MaximumPowerOutput[g]
        else:
            return temp + m.MinimumPowerOutput[g]
    model.ScaledStartupRampLimit = Param(model.ThermalGenerators, within=NonNegativeReals, validate=ramp_limit_validator, initialize=scale_startup_limit, mutable=True)
    
    def scale_shutdown_limit(m, g):
        ## temp now has the "running room" over Pmin. This will be scaled for the time period length
        ## most market models do not have this notion, so this is set-up so that the defaults
        ## will be scaled as they would be in most market models
        temp = (m.ShutdownRampLimit[g] - m.MinimumPowerOutput[g])*m.TimePeriodLengthHours
        if value(temp) > value(m.MaximumPowerOutput[g] - m.MinimumPowerOutput[g]):
            return m.MaximumPowerOutput[g]
        else:
            return temp + m.MinimumPowerOutput[g]
    model.ScaledShutdownRampLimit = Param(model.ThermalGenerators, within=NonNegativeReals, validate=ramp_limit_validator, initialize=scale_shutdown_limit, mutable=True)
    
    
    ##########################################################################################################
    # the minimum number of time periods that a generator must be on-line (off-line) once brought up (down). #
    ##########################################################################################################
    
    model.MinimumUpTime = Param(model.ThermalGenerators,
                                    within=NonNegativeReals,
                                    default=0,
                                    initialize=thermal_gen_attrs['min_up_time'])
    model.MinimumDownTime = Param(model.ThermalGenerators,
                                    within=NonNegativeReals,
                                    default=0,
                                    initialize=thermal_gen_attrs['min_down_time'])
    
    ## Assert that MUT and MDT are at least 1 in the time units of the model.
    ## Otherwise, turn on/offs may not be enforced correctly.
    def scale_min_uptime(m, g):
        scaled_up_time = int(round(m.MinimumUpTime[g] / m.TimePeriodLengthHours))
        return min(max(scaled_up_time,1), value(m.NumTimePeriods))
    model.ScaledMinimumUpTime = Param(model.ThermalGenerators, within=NonNegativeIntegers, initialize=scale_min_uptime)
    
    def scale_min_downtime(m, g):
        scaled_down_time = int(round(m.MinimumDownTime[g] / m.TimePeriodLengthHours))
        return min(max(scaled_down_time,1), value(m.NumTimePeriods))
    model.ScaledMinimumDownTime = Param(model.ThermalGenerators, within=NonNegativeIntegers, initialize=scale_min_downtime)
    
    #############################################
    # unit on state at t=0 (initial condition). #
    #############################################
    
    # if positive, the number of hours prior to (and including) t=0 that the unit has been on.
    # if negative, the number of hours prior to (and including) t=0 that the unit has been off.
    # the value cannot be 0, by definition.
    
    def t0_state_nonzero_validator(m, v, g):
        return v != 0
    
    model.UnitOnT0State = Param(model.ThermalGenerators,
                                within=Reals,
                                validate=t0_state_nonzero_validator,
                                mutable=True,
                                initialize=thermal_gen_attrs['initial_status'])
    
    def t0_unit_on_rule(m, g):
        return int(value(m.UnitOnT0State[g]) >= 1)
    
    model.UnitOnT0 = Param(model.ThermalGenerators,
                            within=Binary,
                            initialize=t0_unit_on_rule,
                            mutable=True)
    
    _add_initial_time_periods_on_off_line(model)
    _verify_must_run_t0_state_consistency(model)
    
    ####################################################################
    # generator power output at t=0 (initial condition). units are MW. #
    ####################################################################
    
    def between_limits_validator(m, v, g):
       status = (v <= (value(m.MaximumPowerOutput[g]) * value(m.UnitOnT0[g]))  and v >= (value(m.MinimumPowerOutput[g]) * value(m.UnitOnT0[g])))
       if status == False:
          print("Failed to validate PowerGeneratedT0 value for g="+g+"; new value="+str(v)+", UnitOnT0="+str(value(m.UnitOnT0[g])))
       return v <= (value(m.MaximumPowerOutput[g]) * value(m.UnitOnT0[g]))  and v >= (value(m.MinimumPowerOutput[g]) * value(m.UnitOnT0[g]))
    model.PowerGeneratedT0 = Param(model.ThermalGenerators, 
                                    within=NonNegativeReals, 
                                    validate=between_limits_validator, 
                                    mutable=True,
                                    initialize=thermal_gen_attrs['initial_p_output'])
    
    
    ###############################################
    # startup cost parameters for each generator. #
    ###############################################
    
    # startup costs are conceptually expressed as pairs (x, y), where x represents the number of hours that a unit has been off and y represents
    # the cost associated with starting up the unit after being off for x hours. these are broken into two distinct ordered sets, as follows.
    
    def startup_lags_init_rule(m, g):
        startup_cost = thermal_gens[g].get('startup_cost')
        startup_fuel = thermal_gens[g].get('startup_fuel')
        if startup_cost is not None and startup_fuel is not None:
            print("WARNING: found startup_fuel for generator {}, ignoring startup_cost".format(g))
        if startup_fuel is None and startup_cost is None:
            return [value(m.MinimumDownTime[g])] 
        elif startup_cost is None:
            return [i[0] for i in startup_fuel]
        else:
            return [i[0] for i in startup_cost]
    model.StartupLags = Set(model.ThermalGenerators, within=NonNegativeReals, ordered=True, initialize=startup_lags_init_rule) # units are hours / time periods.
    
    def startup_costs_init_rule(m, g):
        startup_cost = thermal_gens[g].get('startup_cost')
        startup_fuel = thermal_gens[g].get('startup_fuel')
        fixed_startup_cost = thermal_gens[g].get('non_fuel_startup_cost')
        if fixed_startup_cost is None:
            fixed_startup_cost = 0.
        if startup_fuel is None and startup_cost is None:
            return [fixed_startup_cost]
        elif startup_cost is None:
            fuel_cost = thermal_gens[g].get('fuel_cost')
            if fuel_cost is None:
                raise Exception("No fuel cost for generator {}, but data is provided for fuel tracking".format(g))
            return [fixed_startup_cost+fuel_cost*i[1] for i in startup_fuel]
        else:
            return [fixed_startup_cost+i[1] for i in startup_cost]
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
    
    ## scale the startup lags
    ## Again, assert that this must be at least one in the time units of the model
    def scaled_startup_lags_rule(m, g):
        return [ max(int(round(this_lag / m.TimePeriodLengthHours)),1) for this_lag in m.StartupLags[g] ]
    model.ScaledStartupLags = Set(model.ThermalGenerators, within=NonNegativeIntegers, ordered=True, initialize=scaled_startup_lags_rule)

    ##################################################################################
    # shutdown cost for each generator. in the literature, these are often set to 0. #
    ##################################################################################
    
    model.ShutdownFixedCost = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0, initialize=thermal_gen_attrs.get('shutdown_cost', dict())) # units are $.
    
    ## BEGIN PRODUCTION COST
    ## NOTE: For better or worse, we handle scaling this to the time period length in the objective function.
    ##       In particular, this is done in objective.py.

    ##################################################################################################################
    # production cost coefficients (for the quadratic) a0=constant, a1=linear coefficient, a2=quadratic coefficient. #
    ##################################################################################################################

    def _init_A(coeff):
        def _init(m,g):
            cost = thermal_gens[g].get('p_cost')
            if cost is None:
                return 0.
            elif cost['data_type'] != 'cost_curve':
                raise Exception("p_cost must be of data_type cost_curve.")
            elif cost['cost_curve_type'] == 'piecewise':
                return 0.
            elif cost['cost_curve_type'] == 'polynomial':
                values = cost['values']
                if set(values.keys()) <= {0,1,2}:
                    if coeff in values:
                        return values[coeff]
                    else:
                        return 0.
                else:
                    raise Exception("Polynomial cost curves must be quatric.")
            else:
                raise Exception("Unexpected cost_curve_type")
        return _init
    
    model.ProductionCostA0 = Param(model.ThermalGenerators, default=0.0, initialize=_init_A(0)) # units are $/hr (or whatever the time unit is).
    model.ProductionCostA1 = Param(model.ThermalGenerators, default=0.0, initialize=_init_A(1)) # units are $/MWhr.
    model.ProductionCostA2 = Param(model.ThermalGenerators, default=0.0, initialize=_init_A(2)) # units are $/(MWhr^2).
    
    # the parameters below are populated if cost curves are specified as linearized heat rate increment segments.
    #
    # CostPiecewisePoints represents the power output levels defining the segment boundaries.
    # these *must* include the minimum and maximum power output points - a validation check
    # if performed below.
    # 
    # CostPiecewiseValues are the absolute costs associated with the corresponding
    # power output levels.
    
    # there are many ways to interpret the cost piecewise point/value data, when translating into
    # an actual piecewise construct for the model. this interpretation is controlled by the following
    # string parameter, whose legal values are: NoPiecewise (no data provided) and Absolute.
    # NoPiecewise means that we're using quadraic cost curves, and will 
    # construct the piecewise data ourselves directly from that cost curve. 
    
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
    
    model.PiecewiseType = Param(validate=piecewise_type_validator,initialize=piecewise_type_init, mutable=True)  #irios: default="Absolute" initialize=piecewise_type_init

    def get_piecewise_tuple_index(g, tuple_index):
        cost = thermal_gens[g].get('p_cost')
        fuel = thermal_gens[g].get('p_fuel')
        fuel_cost = thermal_gens[g].get('fuel_cost')
        fixed_no_load = thermal_gens[g].get('non_fuel_no_load_cost')
        if cost is None and fuel is None and fixed_no_load is None:
            return list()
        if fixed_no_load is None or (tuple_index == 0): ## don't add for cost piecewise points
            fixed_no_load = 0.
        if cost is not None and fuel is not None:
            print("WARNING: ignoring provided p_cost and using fuel cost data from p_fuel for generator {}".format(g))
        if fuel is None:
            if cost['data_type'] != 'cost_curve':
                raise Exception("p_cost must be of data_type cost_curve.")
            elif cost['cost_curve_type'] == 'polynomial':
                return list()
            elif cost['cost_curve_type'] == 'piecewise':
                return (fixed_no_load + i[tuple_index] for i in cost['values'])
            else:
                raise Exception("Unexpected cost_curve type")
        else:
            if fuel['data_type'] != 'fuel_curve':
                raise Exception("p_cost must be of data_type fuel_curve for generator {}".format(g))
            if fuel_cost is None:
                raise Exception("must supply fuel costs for generator {} with p_fuel".format(g))
            if tuple_index == 0: ## don't multiply for cost piecewise points
                return (i[tuple_index] for i in fuel['values'])
            return (fixed_no_load + fuel_cost*i[tuple_index] for i in fuel['values'])


    
    def piecewise_points_init(m, g):
        return get_piecewise_tuple_index(g, 0)

    def piecewise_values_init(m, g):
        return get_piecewise_tuple_index(g, 1)
    
    model.CostPiecewisePoints = Set(model.ThermalGenerators, initialize=piecewise_points_init, ordered=True, within=NonNegativeReals)
    model.CostPiecewiseValues = Set(model.ThermalGenerators, initialize=piecewise_values_init, ordered=True, within=NonNegativeReals)
    
    # a check to ensure that the cost piecewise point parameter was correctly populated.
    # these are global checks, which cannot be performed as a set validation (which 
    # operates on a single element at a time).
    
    # irios: When the check fails, I add the missing PiecewisePoints and Values in order to make it work.
    # I did this in an arbitrary way, but you can change it. In particular, I erased those values which are not 
    # between the minimum power output and the maximum power output. Also, I added those values if they are not in
    # the input data. Finally, I added values (0 and this_generator_piecewise_values[-1] + 1) to end with the same 
    # number of points and values.
    
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
            raise Exception("Cost piecewise points for generator g="+str(g)+
                            " must contain the minimum output level="+str(min_output))
    
        if max_output not in points:
            raise Exception("Cost piecewise points for generator g="+str(g)+
                            " must contain the maximum output level="+str(max_output))
        return True

    model.ValidateCostPiecewisePoints = BuildCheck(model.ThermalGenerators, rule=validate_cost_piecewise_points_and_values_rule)
    
    # Minimum production cost (needed because Piecewise constraint on ProductionCost 
    # has to have lower bound of 0, so the unit can cost 0 when off -- this is added
    # back in to the objective if a unit is on
    def minimum_production_cost(m, g):
        if len(m.CostPiecewisePoints[g]) > 1:
            return m.CostPiecewiseValues[g].first()
        else:
            return (m.ProductionCostA0[g] + \
                    m.ProductionCostA1[g] * m.MinimumPowerOutput[g] + \
                    m.ProductionCostA2[g] * (m.MinimumPowerOutput[g]**2))
    model.MinimumProductionCost = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=minimum_production_cost, mutable=True)
    
    ##############################################################################################
    # number of pieces in the linearization of each generator's quadratic cost production curve. #
    ##############################################################################################
    
    model.NumGeneratorCostCurvePieces = Param(within=PositiveIntegers, default=2, mutable=True)


    #######################################################################
    # points for piecewise linearization of power generation cost curves. #
    #######################################################################
    
    # BK -- changed to reflect that the generator's power output variable is always above minimum in the ME model
    #       this simplifies things quite a bit..
    
    # maps a (generator, time-index) pair to a list of points defining the piecewise cost linearization breakpoints.
    # the time index is redundant, but required - in the current implementation of the Piecewise construct, the 
    # breakpoints must be indexed the same as the Piecewise construct itself.
    
    # the points are expected to be on the interval [0, maxpower-minpower], and must contain both endpoints. 
    # power generated can always be 0, and piecewise expects the entire variable domain to be represented.
    model.PowerGenerationPiecewisePoints = {}
    
    # NOTE: the values are relative to the minimum production cost, i.e., the values represent
    # incremental costs relative to the minimum production cost.
    
    model.PowerGenerationPiecewiseValues = {}
    
    def power_generation_piecewise_points_rule(m, g, t):
    
        # factor out the fuel cost here, as the piecewise approximation is scaled by fuel cost
        # elsewhere in the model (i.e., in the Piecewise construct below).
        minimum_production_cost = value(m.MinimumProductionCost[g])
    
        # minimum output
        minimum_power_output = value(m.MinimumPowerOutput[g])
        
        piecewise_type = value(m.PiecewiseType)
    
        if piecewise_type == "Absolute":
    
           piecewise_values = list(m.CostPiecewiseValues[g])
           piecewise_points = list(m.CostPiecewisePoints[g])
           m.PowerGenerationPiecewiseValues[g,t] = {}
           m.PowerGenerationPiecewisePoints[g,t] = [] 
           for i in range(len(piecewise_points)):
              this_point = piecewise_points[i] - minimum_power_output
              m.PowerGenerationPiecewisePoints[g,t].append(this_point)
              m.PowerGenerationPiecewiseValues[g,t][this_point] = piecewise_values[i] - minimum_production_cost
    
        else: # piecewise_type == "NoPiecewise"
    
           if value(m.ProductionCostA2[g]) == 0:
              # If cost is linear, we only need two points -- (0,0) and (MaxOutput-MinOutput, MaxCost-MinCost))
              min_power = value(m.MinimumPowerOutput[g])
              max_power = value(m.MaximumPowerOutput[g])
              if min_power == max_power:
                 m.PowerGenerationPiecewisePoints[g, t] = [0.0]
              else:
                 m.PowerGenerationPiecewisePoints[g, t] = [0.0, max_power-min_power]
    
              m.PowerGenerationPiecewiseValues[g,t] = {}
    
              m.PowerGenerationPiecewiseValues[g,t][0.0] = 0.0
    
              if min_power != max_power:
                 m.PowerGenerationPiecewiseValues[g,t][max_power-min_power] = \
                     value(m.ProductionCostA0[g]) + \
                     value(m.ProductionCostA1[g]) * max_power \
                     - minimum_production_cost
    
           else:
               min_power = value(m.MinimumPowerOutput[g])
               max_power = value(m.MaximumPowerOutput[g])
               n = value(m.NumGeneratorCostCurvePieces)
               width = (max_power - min_power) / float(n)
               if width == 0:
                   m.PowerGenerationPiecewisePoints[g, t] = [0]
               else:
                   m.PowerGenerationPiecewisePoints[g, t] = []
                   m.PowerGenerationPiecewisePoints[g, t].extend([0 + i*width for i in range(0,n+1)])
                   # NOTE: due to numerical precision limitations, the last point in the x-domain
                   #       of the generation piecewise cost curve may not be precisely equal to the 
                   #       maximum power output level of the generator. this can cause Piecewise to
                   #       sqawk, as it would like the upper bound of the variable to be represented
                   #       in the domain. so, we will make it so.
                   m.PowerGenerationPiecewisePoints[g, t][-1] = max_power - min_power
               m.PowerGenerationPiecewiseValues[g,t] = {}
               for i in range(len(m.PowerGenerationPiecewisePoints[g, t])):
                   m.PowerGenerationPiecewiseValues[g,t][m.PowerGenerationPiecewisePoints[g,t][i]] = \
                              value(m.ProductionCostA0[g]) + \
                              value(m.ProductionCostA1[g]) * (m.PowerGenerationPiecewisePoints[g, t][i] + min_power) + \
                              value(m.ProductionCostA2[g]) * (m.PowerGenerationPiecewisePoints[g, t][i] + min_power)**2 \
                              - minimum_production_cost
               assert(m.PowerGenerationPiecewisePoints[g, t][0] == 0)
        
        # validate the computed points, independent of the method used to generate them.
        # nothing should be negative, and the costs should be monotonically non-decreasing.
        for i in range(0, len(m.PowerGenerationPiecewisePoints[g, t])):
           this_level = m.PowerGenerationPiecewisePoints[g, t][i]
           assert this_level >= 0.0
    
    model.CreatePowerGenerationPiecewisePoints = BuildAction(model.ThermalGenerators * model.TimePeriods, rule=power_generation_piecewise_points_rule)

    ModeratelyBigPenalty = 1e3*system['baseMVA']
    
    model.ReserveShortfallPenalty = Param(within=NonNegativeReals, default=ModeratelyBigPenalty, mutable=True, initialize=system.get('reserve_shortfall_cost'))

    #########################################
    # penalty costs for constraint violation #
    #########################################
    
    BigPenalty = 1e4*system['baseMVA']
    
    model.LoadMismatchPenalty = Param(within=NonNegativeReals, default=BigPenalty, mutable=True, initialize=system.get('load_mismatch_cost'))
    model.LoadMismatchPenaltyReactive = Param(within=NonNegativeReals, default=BigPenalty/2., mutable=True, initialize=system.get('q_load_mismatch_cost'))

    ## END PRODUCTION COST CALCULATIONS

    ## FUEL-SUPPLY Sets

    def fuel_supply_gens_init(m):
        if 'fuel_supply' not in elements and ('fuel_supply' in thermal_gen_attrs or 'aux_fuel_supply' in thermal_gen_attrs):
            print('WARNING: Some generators have \'fuel_supply\' marked, but no fuel supply was found on ModelData.data[\'system\']')
            return iter(())
        if 'fuel_supply' in elements and ('fuel_supply' not in thermal_gen_attrs and 'aux_fuel_supply' not in thermal_gen_attrs):
            print('WARNING: fuel_supply in ModelData.data["elements"], but no generators are attached to any fuel supply')
            return iter(())
        if 'fuel_supply' not in thermal_gen_attrs:
            thermal_gen_attrs['fuel_supply'] = dict()
        if 'aux_fuel_supply' not in thermal_gen_attrs:
            thermal_gen_attrs['aux_fuel_supply'] = dict()
        gen_set = set(thermal_gen_attrs['fuel_supply'].keys())
        gen_set.update(thermal_gen_attrs['aux_fuel_supply'].keys())
        return gen_set

    def gen_cost_fuel_validator(m,g):
        if 'p_fuel' in thermal_gen_attrs and g in thermal_gen_attrs['p_fuel']:
            pass
        else:
            print('ERROR: All fuel-constrained generators must have "p_fuel" attribute which tracks their fuel consumption')
            print('ERROR: Could not find such an attribute for generator {}'.format(g))
            return False
        return True

    model.FuelSupplyGenerators = Set(within=model.ThermalGenerators, initialize=fuel_supply_gens_init, validate=gen_cost_fuel_validator)

    ## DUAL-FUEL Sets

    def dual_fuel_init(m):
        for g, g_dict in thermal_gens.items():
            if 'aux_fuel_capable' in g_dict and g_dict['aux_fuel_capable']:
                yield g
    model.DualFuelGenerators = Set(within=model.ThermalGenerators, initialize=dual_fuel_init)

    ## This set is for modeling elements that are exhanged
    ## in whole for the dual-fuel model
    model.SingleFuelGenerators = model.ThermalGenerators - model.DualFuelGenerators

    #
    # STORAGE parameters
    #
    


    
    model.Storage = Set(initialize=storage_attrs['names'])
    model.StorageAtBus = Set(model.Buses, initialize=storage_by_bus)
    
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
    
    model.MinimumPowerOutputStorage = Param(model.Storage, within=NonNegativeReals,
                                            default=0.0, initialize=storage_attrs.get('min_discharge_rate'))
    
    def maximum_power_output_validator_storage(m, v, s):
        return v >= value(m.MinimumPowerOutputStorage[s])
    
    model.MaximumPowerOutputStorage = Param(model.Storage, within=NonNegativeReals,
                                            validate=maximum_power_output_validator_storage, default=0.0,
                                            initialize=storage_attrs.get('max_discharge_rate'))
    
    #Storage power input >0 when charging
    
    model.MinimumPowerInputStorage = Param(model.Storage, within=NonNegativeReals,
                                            default=0.0, initialize=storage_attrs.get('min_charge_rate'))
    
    def maximum_power_input_validator_storage(m, v, s):
        return v >= value(m.MinimumPowerInputStorage[s])
    
    model.MaximumPowerInputStorage = Param(model.Storage, within=NonNegativeReals,
                                            validate=maximum_power_input_validator_storage, default=0.0,
                                            initialize=storage_attrs.get('max_charge_rate'))
    
    ###############################################
    # storage ramp up/down rates. units are MW/h. #
    ###############################################
    
    # ramp rate limits when discharging
    model.NominalRampUpLimitStorageOutput    = Param(model.Storage, within=NonNegativeReals,
                                                        initialize=storage_attrs.get('ramp_up_output_60min'))
    model.NominalRampDownLimitStorageOutput  = Param(model.Storage, within=NonNegativeReals,
                                                        initialize=storage_attrs.get('ramp_down_output_60min'))
    
    # ramp rate limits when charging
    model.NominalRampUpLimitStorageInput     = Param(model.Storage, within=NonNegativeReals,
                                                        initialize=storage_attrs.get('ramp_up_input_60min'))
    model.NominalRampDownLimitStorageInput   = Param(model.Storage, within=NonNegativeReals,
                                                        initialize=storage_attrs.get('ramp_down_input_60min'))
    
    def scale_storage_ramp_up_out(m, s):
        return m.NominalRampUpLimitStorageOutput[s] * m.TimePeriodLengthHours
    model.ScaledNominalRampUpLimitStorageOutput = Param(model.Storage, within=NonNegativeReals, initialize=scale_storage_ramp_up_out)
    
    def scale_storage_ramp_down_out(m, s):
        return m.NominalRampDownLimitStorageOutput[s] * m.TimePeriodLengthHours
    model.ScaledNominalRampDownLimitStorageOutput = Param(model.Storage, within=NonNegativeReals, initialize=scale_storage_ramp_down_out)
    
    def scale_storage_ramp_up_in(m, s):
        return m.NominalRampUpLimitStorageInput[s] * m.TimePeriodLengthHours
    model.ScaledNominalRampUpLimitStorageInput = Param(model.Storage, within=NonNegativeReals, initialize=scale_storage_ramp_up_in)
    
    def scale_storage_ramp_down_in(m, s):
        return m.NominalRampDownLimitStorageInput[s] * m.TimePeriodLengthHours
    model.ScaledNominalRampDownLimitStorageInput = Param(model.Storage, within=NonNegativeReals, initialize=scale_storage_ramp_down_in)
    
    ####################################################################################
    # minimum state of charge (SOC) and maximum energy ratings, for each storage unit. #
    # units are MWh for energy rating and p.u. (i.e. [0,1]) for SOC     #
    ####################################################################################
    
    # you enter storage energy ratings once for each storage unit
    
    model.MaximumEnergyStorage = Param(model.Storage, within=NonNegativeReals, default=0.0,
                                        initialize=storage_attrs.get('energy_capacity'))
    model.MinimumSocStorage = Param(model.Storage, within=PercentFraction, default=0.0,
                                        initialize=storage_attrs.get('minimum_state_of_charge'))
    
    ################################################################################
    # round trip efficiency for each storage unit given as a fraction (i.e. [0,1]) #
    ################################################################################
    
    model.InputEfficiencyEnergy  = Param(model.Storage, within=PercentFraction, default=1.0,
                                            initialize=storage_attrs.get('charge_efficiency'))
    model.OutputEfficiencyEnergy = Param(model.Storage, within=PercentFraction, default=1.0,
                                            initialize=storage_attrs.get('discharge_efficienty'))
    model.RetentionRate          = Param(model.Storage, within=PercentFraction, default=1.0,
                                            initialize=storage_attrs.get('retention_rate_60min')) ## assumed to be %/hr

    model.ChargeCost = Param(model.Storage, within=Reals, default=0.0, initialize=storage_attrs.get('charge_cost'))
    model.DischargeCost = Param(model.Storage, within=Reals, default=0.0, initialize=storage_attrs.get('discharge_cost'))

    ## this will be multiplied by itself 1/m.TimePeriodLengthHours times, so this is the scaling to
    ## get us back to %/hr
    def scaled_retention_rate(m,s):
        return value(m.RetentionRate[s])**value(m.TimePeriodLengthHours)
    model.ScaledRetentionRate = Param(model.Storage, within=PercentFraction, initialize=scaled_retention_rate)
    
    ########################################################################
    # end-point SOC for each storage unit. units are in p.u. (i.e. [0,1])  #
    ########################################################################
    
    # end-point values are the SOC targets at the final time period. With no end-point constraints
    # storage units will always be empty at the final time period.
    
    model.EndPointSocStorage = Param(model.Storage, within=PercentFraction, default=0.5,
                                        initialize=storage_attrs.get('initial_state_of_charge'))
    
    ############################################################
    # storage initial conditions: SOC, power output and input  #
    ############################################################
    
    def t0_storage_power_input_validator(m, v, s):
        return (v >= value(m.MinimumPowerInputStorage[s])) and (v <= value(m.MaximumPowerInputStorage[s]))
    
    def t0_storage_power_output_validator(m, v, s):
        return (v >= value(m.MinimumPowerInputStorage[s])) and (v <= value(m.MaximumPowerInputStorage[s]))
    
    model.StoragePowerOutputOnT0 = Param(model.Storage, within=NonNegativeReals,
                                            validate=t0_storage_power_output_validator,
                                            default=0.0,
                                            initialize=storage_attrs.get('initial_discharge_rate'))
    model.StoragePowerInputOnT0  = Param(model.Storage, within=NonNegativeReals,
                                            validate=t0_storage_power_input_validator,
                                            default=0.0,
                                            initialize=storage_attrs.get('initial_charge_rate'))
    model.StorageSocOnT0         = Param(model.Storage, within=PercentFraction,
                                            default=0.5, initialize=storage_attrs.get('initial_state_of_charge'))

    ##############################################################
    # failure probability for each generator, in any given hour. #
    # not used within the model itself at present, but rather    #
    # used by scripts that read / manipulate the model.          #
    ##############################################################
    
    def probability_failure_validator(m, v, g):
       return v >= 0.0 and v <= 1.0
    
    model.FailureProbability = Param(model.ThermalGenerators, validate=probability_failure_validator, default=0.0)
    
    return model
