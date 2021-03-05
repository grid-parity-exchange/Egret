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
from egret.data.data_utils import map_items, zip_items
from egret.model_library.transmission import tx_utils
from egret.common.log import logger
    
from .uc_utils import add_model_attr, uc_time_helper

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
                        print("DATA ERROR: The generator %s has been flagged as off at time %d, but its T0 state=%d is inconsistent with its minimum up time=%d" % (g, t, t0_state, min_up_time))
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

    time_keys = system['time_keys']
    
    ## insert potentially missing keys
    if 'branch' not in elements:
        elements['branch'] = dict()
    if 'interface' not in elements:
        elements['interface'] = dict()
    if 'storage' not in elements:
        elements['storage'] = dict()
    if 'dc_branch' not in elements:
        elements['dc_branch'] = dict()

    ## NOTE: generator, bus, and load should be in here for a well-defined problem

    loads = dict(md.elements(element_type='load'))
    thermal_gens = dict(md.elements(element_type='generator', generator_type='thermal'))
    renewable_gens = dict(md.elements(element_type='generator', generator_type='renewable'))
    buses = dict(md.elements(element_type='bus'))
    shunts = dict(md.elements(element_type='shunt'))
    branches = dict(md.elements(element_type='branch'))
    interfaces = dict(md.elements(element_type='interface'))
    contingencies = dict(md.elements(element_type='contingency'))
    storage = dict(md.elements(element_type='storage'))
    dc_branches = dict(md.elements(element_type='dc_branch'))

    thermal_gen_attrs = md.attributes(element_type='generator', generator_type='thermal')
    renewable_gen_attrs = md.attributes(element_type='generator', generator_type='renewable')
    bus_attrs = md.attributes(element_type='bus')
    branch_attrs = md.attributes(element_type='branch')
    load_attrs = md.attributes(element_type='load')
    interface_attrs = md.attributes(element_type='interface')
    storage_attrs = md.attributes(element_type='storage')
    dc_branch_attrs = md.attributes(element_type='dc_branch')


    inlet_branches_by_bus, outlet_branches_by_bus = \
        tx_utils.inlet_outlet_branches_by_bus(branches, buses)
    dc_inlet_branches_by_bus, dc_outlet_branches_by_bus = \
        tx_utils.inlet_outlet_branches_by_bus(dc_branches, buses)
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
    model._interfaces = interfaces
    model._contingencies = contingencies
    model._dc_branches = dc_branches
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
        reference_bus = list(sorted(model.Buses))[0]

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

    model.NumTimePeriods = Param(within=PositiveIntegers, initialize=len(system['time_keys']))
    
    model.InitialTime = Param(within=PositiveIntegers, default=1)
    model.TimePeriods = RangeSet(model.InitialTime, model.NumTimePeriods)

    TimeMapper = uc_time_helper(model.TimePeriods)
    
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
    model.HVDCLines = Set(initialize=dc_branch_attrs['names'])
    
    model.BusFrom = Param(model.TransmissionLines, within=model.Buses, initialize=branch_attrs.get('from_bus', dict()))
    model.BusTo   = Param(model.TransmissionLines, within=model.Buses, initialize=branch_attrs.get('to_bus', dict()))

    model.HVDCBusFrom = Param(model.HVDCLines, within=model.Buses, initialize=dc_branch_attrs.get('from_bus', dict()))
    model.HVDCBusTo   = Param(model.HVDCLines, within=model.Buses, initialize=dc_branch_attrs.get('to_bus', dict()))

    model.LinesTo = Set(model.Buses, initialize=inlet_branches_by_bus)
    model.LinesFrom = Set(model.Buses, initialize=outlet_branches_by_bus)

    model.HVDCLinesTo = Set(model.Buses, initialize=dc_inlet_branches_by_bus)
    model.HVDCLinesFrom = Set(model.Buses, initialize=dc_outlet_branches_by_bus)

    def _warn_neg_impedence(m, v, l):
        if v == 0.:
            logger.error(f"Found zero reactance for line {l}")
            return False
        elif v < 0.:
            # We allow negative reactance, as it just reverses the
            # direction of the line. But we do print a warning.
            logger.warning(f"WARNING: found negative reactance for line {l}")
            return True
        return True
    model.Impedence = Param(model.TransmissionLines, within=Reals, initialize=branch_attrs.get('reactance', dict()), validate=_warn_neg_impedence)

    model.ThermalLimit = Param(model.TransmissionLines, initialize=branch_attrs.get('rating_long_term', dict())) # max flow across the line
    model.HVDCThermalLimit = Param(model.HVDCLines, initialize=dc_branch_attrs.get('rating_long_term', dict())) # max flow across the line

    model.LineOutOfService = Param(model.TransmissionLines, model.TimePeriods, within=Boolean, default=False,
                                    initialize=TimeMapper(branch_attrs.get('planned_outage', dict())))

    model.HVDCLineOutOfService = Param(model.HVDCLines, model.TimePeriods, within=Boolean, default=False,
                                       initialize=TimeMapper(dc_branch_attrs.get('planned_outage', dict())))

    _branch_penalties = dict()
    _md_violation_penalties = branch_attrs.get('violation_penalty')
    if _md_violation_penalties is not None:
        for i, val in _md_violation_penalties.items():
            if val is not None:
                _branch_penalties[i] = val
                if val <= 0:
                    logger.warning("Branch {} has a non-positive penalty {}, this will cause its limits to be ignored!".format(i,val))

    model.BranchesWithSlack = Set(within=model.TransmissionLines, initialize=_branch_penalties.keys())

    model.BranchLimitPenalty = Param(model.BranchesWithSlack, within=NonNegativeReals, initialize=_branch_penalties)

    ## Interfaces
    model.Interfaces = Set(initialize=interface_attrs['names'])

    model.InterfaceLines = Set(model.Interfaces, within=model.TransmissionLines, initialize=interface_attrs.get('lines', dict()), ordered=True)
    model.InterfaceMinFlow = Param(model.Interfaces, within=Reals, initialize=interface_attrs.get('minimum_limit', dict()))
    model.InterfaceMaxFlow = Param(model.Interfaces, within=Reals, initialize=interface_attrs.get('maximum_limit', dict()))

    def check_min_less_max_interface_flow_limits(m):
        for i in m.Interfaces:
            if value(m.InterfaceMinFlow[i]) > value(m.InterfaceMaxFlow[i]):
                raise Exception("Interface {} has a minimum_limit which is greater than the maximum_limit".format(i))

    model.CheckInterfaceFlowLimits = BuildAction(rule=check_min_less_max_interface_flow_limits)

    def get_interface_line_pairs(m):
        for i in m.Interfaces:
            for l in m.InterfaceLines[i]:
                yield i,l
    model.InterfaceLinePairs = Set(initialize=get_interface_line_pairs, dimen=2)

    _interface_line_orientation_dict = dict()
    for i, interface in interfaces.items():
        for l, sign in zip(interface['lines'],interface['line_orientation']):
            _interface_line_orientation_dict[i,l] = sign

    model.InterfaceLineOrientation = Param(model.InterfaceLinePairs, initialize=_interface_line_orientation_dict, within=set([-1,0,1]))

    _interface_penalties = dict()
    _md_violation_penalties = interface_attrs.get('violation_penalty')
    if _md_violation_penalties is not None:
        for i, val in _md_violation_penalties.items():
            if val is not None:
                _interface_penalties[i] = val
                if val <= 0:
                    logger.warning("Interface {} has a non-positive penalty {}, this will cause its limits to be ignored!".format(i,val))

    model.InterfacesWithSlack = Set(within=model.Interfaces, initialize=_interface_penalties.keys())

    model.InterfaceLimitPenalty = Param(model.InterfacesWithSlack, within=NonNegativeReals, initialize=_interface_penalties)
  
    ##########################################################
    # string indentifiers for the set of thermal generators. #
    # and their locations. (S)                               #
    ##########################################################
    
    model.ThermalGenerators = Set(initialize=thermal_gen_attrs['names'])
    model.ThermalGeneratorsAtBus = Set(model.Buses, initialize=thermal_gens_by_bus)
    
    model.ThermalGeneratorType = Param(model.ThermalGenerators, within=Any, default='C', initialize=thermal_gen_attrs.get('fuel', dict()))
    
    def verify_thermal_generator_buses_rule(m, g):
       for b in m.Buses:
          if g in m.ThermalGeneratorsAtBus[b]:
             return 
       print("DATA ERROR: No bus assigned for thermal generator=%s" % g)
       assert(False)
    
    model.VerifyThermalGeneratorBuses = BuildAction(model.ThermalGenerators, rule=verify_thermal_generator_buses_rule)
    
    model.QuickStart = Param(model.ThermalGenerators, within=Boolean, default=False, initialize=thermal_gen_attrs.get('fast_start', dict()))
    
    def init_quick_start_generators(m):
        return [g for g in m.ThermalGenerators if value(m.QuickStart[g]) == 1]
    
    model.QuickStartGenerators = Set(within=model.ThermalGenerators, initialize=init_quick_start_generators)
    
    # optionally force a unit to be on/off
    model.FixedCommitmentTypes = Set(initialize=[0,1,None])
    model.FixedCommitment = Param(model.ThermalGenerators,
                                  model.TimePeriods,
                                  within=model.FixedCommitmentTypes,
                                  default=None,
                                  initialize=TimeMapper(thermal_gen_attrs.get('fixed_commitment', dict())),)
    
    model.NondispatchableGeneratorsAtBus = Set(model.Buses, initialize=renewable_gens_by_bus)
    
    model.AllNondispatchableGenerators = Set(initialize=renewable_gen_attrs['names'])

    model.NondispatchableGeneratorType = Param(model.AllNondispatchableGenerators, within=Any, default='W', 
                                                initialize=renewable_gen_attrs.get('fuel', dict()))
    
    
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
        load_time = TimeMapper(load['p_load'])
        bus = load['bus']
        if isinstance(bus, dict):
            assert bus['data_type'] == 'load_distribution_factor'
            for bn, multi in bus['values'].items():
                for t in model.TimePeriods:
                    bus_loads[bn, t] += multi*load_time[t]
        else:
            for t in model.TimePeriods:
                bus_loads[bus, t] += load_time[t]
    model.Demand = Param(model.Buses, model.TimePeriods, initialize=bus_loads, mutable=True)
    
    def calculate_total_demand(m, t):
        return sum(value(m.Demand[b,t]) for b in sorted(m.Buses))
    model.TotalDemand = Param(model.TimePeriods, initialize=calculate_total_demand)
    
    # at this point, a user probably wants to see if they have negative demand.
    def warn_about_negative_demand_rule(m, b, t):
       this_demand = value(m.Demand[b,t])
       if this_demand < 0.0:
          logger.warning("***WARNING: The demand at bus="+str(b)+" for time period="+str(t)+" is negative - value="+str(this_demand)+"; model="+str(m.name)+".")
    
    if warn_neg_load:
        model.WarnAboutNegativeDemand = BuildAction(model.Buses, model.TimePeriods, rule=warn_about_negative_demand_rule)
    
    ##################################################################
    # the global system reserve, for each time period. units are MW. #
    ##################################################################

    reserve_requirement = system.get("reserve_requirement", 0.)
    model.ReserveRequirement = Param(model.TimePeriods, within=NonNegativeReals, 
                                        initialize=TimeMapper(reserve_requirement), mutable=True)
    
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
    
    ####################################################################################
    # minimum and maximum generation levels, for each thermal generator. units are MW. #
    # could easily be specified on a per-time period basis, but are not currently.     #
    ####################################################################################
    
    # you can enter generator limits either once for the generator or for each period (or just take 0)
    
    model.MinimumPowerOutput = Param(model.ThermalGenerators, model.TimePeriods, 
                                        within=NonNegativeReals, 
                                        initialize=TimeMapper(thermal_gen_attrs['p_min']),
                                        default=0.0)
    
    def maximum_power_output_validator(m, v, g, t):
       return v >= value(m.MinimumPowerOutput[g,t])
    
    model.MaximumPowerOutput = Param(model.ThermalGenerators, model.TimePeriods, 
                                        within=NonNegativeReals, 
                                        validate=maximum_power_output_validator, 
                                        initialize=TimeMapper(thermal_gen_attrs['p_max']),
                                        default=0.0)

    model.MinimumReactivePowerOutput = Param(model.ThermalGenerators, model.TimePeriods,
                                                within=Reals,
                                                initialize=TimeMapper(thermal_gen_attrs.get('q_min', dict())),
                                                default=0.0)

    def maximum_reactive_output_validator(m, v, g, t):
        return v >= value(m.MinimumReactivePowerOutput[g,t])

    model.MaximumReactivePowerOutput = Param(model.ThermalGenerators, model.TimePeriods, 
                                                within=Reals,
                                                initialize=TimeMapper(thermal_gen_attrs.get('q_max', dict())),
                                                default=0.0)
    
    # wind is similar, but max and min will be equal for non-dispatchable wind
    
    model.MinNondispatchablePower = Param(model.AllNondispatchableGenerators,
                                            model.TimePeriods, 
                                            within=Reals, # more permissive; e.g. CSP
                                            default=0.0,
                                            mutable=True,
                                            initialize=TimeMapper(renewable_gen_attrs.get('p_min', dict())))
    
    def maximum_nd_output_validator(m, v, g, t):
       return v >= value(m.MinNondispatchablePower[g,t])
    
    model.MaxNondispatchablePower = Param(model.AllNondispatchableGenerators,
                                            model.TimePeriods,
                                            within=Reals, # more permissive; e.g. CSP
                                            default=0.0,
                                            mutable=True,
                                            validate=maximum_nd_output_validator,
                                            initialize=TimeMapper(renewable_gen_attrs.get('p_max', dict())))

    #################################################
    # generator ramp up/down rates. units are MW/h. #
    # IMPORTANT: Generator ramp limits can exceed   #
    # the maximum power output, because it is the   #
    # ramp limit over an hour. If the unit can      #
    # fully ramp in less than an hour, then this    #
    # will occur.                                   #
    #################################################

    ## be sure the generator can ramp
    ## between all the p_min/p_max values
    def ramp_up_validator(m, v, g):
        t1 = m.InitialTime
        for t in m.TimePeriods:
            if t == t1:
                continue
            diff = value(m.MinimumPowerOutput[g,t] - m.MaximumPowerOutput[g,t-1])
            if v*m.TimePeriodLengthHours < diff:
                logger.error('Generator {} has an infeasible ramp up between time periods {} and {}'.format(g,t-1,t))
                return False
        return True

    ## be sure the generator can ramp
    ## between all the p_min/p_max values
    def ramp_down_validator(m, v, g):
        t1 = m.InitialTime
        for t in m.TimePeriods:
            if t == t1:
                continue
            diff = value(m.MinimumPowerOutput[g,t-1] - m.MaximumPowerOutput[g,t])
            if v*m.TimePeriodLengthHours < diff:
                logger.error('Generator {} has an infeasible ramp down between time periods {} and {}'.format(g,t-1,t))
                return False
        return True

    # limits for normal time periods
    model.NominalRampUpLimit = Param(model.ThermalGenerators,
                                        within=NonNegativeReals,
                                        mutable=True,
                                        initialize=thermal_gen_attrs['ramp_up_60min'],
                                        validate=ramp_up_validator)
    model.NominalRampDownLimit = Param(model.ThermalGenerators,
                                        within=NonNegativeReals,
                                        mutable=True,
                                        initialize=thermal_gen_attrs['ramp_down_60min'],
                                        validate=ramp_down_validator)

    #############################################
    # unit on state at t=0 (initial condition). #
    #############################################
    
    # if positive, the number of hours prior to (and including) t=0 that the unit has been on.
    # if negative, the number of hours prior to (and including) t=0 that the unit has been off.
    # the value cannot be 0, by definition.
    
    def t0_state_nonzero_validator(m, v, g):
        return v != 0.
    
    model.UnitOnT0State = Param(model.ThermalGenerators,
                                within=Reals,
                                validate=t0_state_nonzero_validator,
                                mutable=True,
                                initialize=thermal_gen_attrs['initial_status'])
    
    def t0_unit_on_rule(m, g):
        return int(value(m.UnitOnT0State[g]) > 0.)
    
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
        t = m.TimePeriods.first() 

        if value(m.UnitOnT0[g]):
            v_less_max = v <= value(m.MaximumPowerOutput[g,t] + m.NominalRampDownLimit[g]*m.TimePeriodLengthHours)
            if not v_less_max:
                logger.error('Generator {} has more output at T0 than is feasible to ramp down to'.format(g))
                return False
            v_greater_min = v >= value(m.MinimumPowerOutput[g,t] - m.NominalRampUpLimit[g]*m.TimePeriodLengthHours)
            if not v_less_max:
                logger.error('Generator {} has less output at T0 than is feasible to ramp up to'.format(g))
                return False
            return True

        else:
            return v == 0.

    model.PowerGeneratedT0 = Param(model.ThermalGenerators, 
                                    within=NonNegativeReals, 
                                    validate=between_limits_validator, 
                                    mutable=True,
                                    initialize=thermal_gen_attrs['initial_p_output'])
    
    # limits for time periods in which generators are brought on or off-line.
    # must be no less than the generator minimum output.
    def ramp_limit_validator(m, v, g, t):
       return v >= m.MinimumPowerOutput[g,t]

    ## These defaults follow what is in most market manuals
    ## We scale this for the time period below
    def startup_ramp_default(m, g, t):
        return m.MinimumPowerOutput[g,t]+m.NominalRampUpLimit[g]/2.

    ## shutdown is based on the last period *on*
    def shutdown_ramp_default(m, g, t):
        return m.MinimumPowerOutput[g,t]+m.NominalRampDownLimit[g]/2.

    model.StartupRampLimit = Param(model.ThermalGenerators, 
                                    model.TimePeriods,
                                    within=NonNegativeReals,
                                    default=startup_ramp_default,
                                    validate=ramp_limit_validator,
                                    mutable=True,
                                    initialize=TimeMapper(thermal_gen_attrs.get('startup_capacity', dict())))
    model.ShutdownRampLimit = Param(model.ThermalGenerators, 
                                    model.TimePeriods,
                                    within=NonNegativeReals,
                                    default=shutdown_ramp_default, 
                                    validate=ramp_limit_validator,
                                    mutable=True,
                                    initialize=TimeMapper(thermal_gen_attrs.get('shutdown_capacity', dict())))
    
    ## These get used in the basic UC constraints, which implicity assume RU, RD <= Pmax
    ## Ramping constraints look backward, so these will accordingly as well
    ## NOTES: can't ramp up higher than the current pmax from the previous value
    ##        can't ramp down more than the pmax from the prior time period
    def scale_ramp_up(m, g, t):
        temp = m.NominalRampUpLimit[g] * m.TimePeriodLengthHours
        if value(temp) > value(m.MaximumPowerOutput[g,t]):
            return m.MaximumPowerOutput[g,t]
        else:
            return temp
    model.ScaledNominalRampUpLimit = Param(model.ThermalGenerators, model.TimePeriods,  within=NonNegativeReals, initialize=scale_ramp_up, mutable=True)
    
    def scale_ramp_down(m, g, t):
        temp = m.NominalRampDownLimit[g] * m.TimePeriodLengthHours
        if t == m.InitialTime:
            param = max(value(m.PowerGeneratedT0[g]), value(m.MaximumPowerOutput[g,t]))
        else:
            param = m.MaximumPowerOutput[g,t-1]
        if value(temp) > value(param):
            return param
        else:
            return temp
    model.ScaledNominalRampDownLimit = Param(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, initialize=scale_ramp_down, mutable=True)
    
    def scale_startup_limit(m, g, t):
        ## temp now has the "running room" over Pmin. This will be scaled for the time period length, 
        ## most market models do not have this notion, so this is set-up so that the defaults
        ## will be scaled as they would be in most market models
        temp = (m.StartupRampLimit[g,t] - m.MinimumPowerOutput[g,t])*m.TimePeriodLengthHours
        if value(temp) > value(m.MaximumPowerOutput[g,t] - m.MinimumPowerOutput[g,t]):
            return m.MaximumPowerOutput[g,t]
        else:
            return temp + m.MinimumPowerOutput[g,t]
    model.ScaledStartupRampLimit = Param(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, validate=ramp_limit_validator, initialize=scale_startup_limit, mutable=True)
    
    def scale_shutdown_limit(m, g, t):
        ## temp now has the "running room" over Pmin. This will be scaled for the time period length
        ## most market models do not have this notion, so this is set-up so that the defaults
        ## will be scaled as they would be in most market models
        temp = (m.ShutdownRampLimit[g,t] - m.MinimumPowerOutput[g,t])*m.TimePeriodLengthHours
        if value(temp) > value(m.MaximumPowerOutput[g,t] - m.MinimumPowerOutput[g,t]):
            return m.MaximumPowerOutput[g,t]
        else:
            return temp + m.MinimumPowerOutput[g,t]
    model.ScaledShutdownRampLimit = Param(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals, validate=ramp_limit_validator, initialize=scale_shutdown_limit, mutable=True)

    ## Some additional ramping parameters to 
    ## deal with shutdowns at time=1
    
    def _init_p_min_t0(m,g):
        if 'initial_p_min' in thermal_gen_attrs and \
                g in thermal_gen_attrs['initial_p_min']:
            return thermal_gen_attrs['initial_p_min'][g]
        else:
            return m.MinimumPowerOutput[g,m.InitialTime]

    model.MinimumPowerOutputT0 = Param(model.ThermalGenerators, 
                                    within=NonNegativeReals, 
                                    mutable=True,
                                    initialize=_init_p_min_t0)

    def _init_sd_t0(m,g):
        if 'initial_shutdown_capacity' in thermal_gen_attrs and\
                g in thermal_gen_attrs['initial_shutdown_capacity']:
            return thermal_gen_attrs['initial_shutdown_capacity'][g]
        return m.ShutdownRampLimit[g,m.InitialTime]

    model.ShutdownRampLimitT0 = Param(model.ThermalGenerators,
                                    within=NonNegativeReals,
                                    mutable=True,
                                    initialize=_init_sd_t0)
    
    def scale_shutdown_limit_t0(m, g):
        ## temp now has the "running room" over Pmin. This will be scaled for the time period length
        ## most market models do not have this notion, so this is set-up so that the defaults
        ## will be scaled as they would be in most market models
        temp = (m.ShutdownRampLimitT0[g] - m.MinimumPowerOutputT0[g])*m.TimePeriodLengthHours
        if value(temp) > value(m.PowerGeneratedT0[g] - m.MinimumPowerOutputT0[g]):
            return m.PowerGeneratedT0[g]
        else:
            return temp + m.MinimumPowerOutputT0[g]
    model.ScaledShutdownRampLimitT0 = Param(model.ThermalGenerators, within=NonNegativeReals, initialize=scale_shutdown_limit_t0, mutable=True)
    
    ###############################################
    # startup cost parameters for each generator. #
    ###############################################
    
    # startup costs are conceptually expressed as pairs (x, y), where x represents the number of hours that a unit has been off and y represents
    # the cost associated with starting up the unit after being off for x hours. these are broken into two distinct ordered sets, as follows.

    def _get_startup_lag(startup,default):
        try:
            iter(startup)
        except TypeError:
            return [default]
        else:
            return [i[0] for i in startup]
    
    def startup_lags_init_rule(m, g):
        startup_cost = thermal_gens[g].get('startup_cost')
        startup_fuel = thermal_gens[g].get('startup_fuel')
        if startup_cost is not None and startup_fuel is not None:
            logger.warning("WARNING: found startup_fuel for generator {}, ignoring startup_cost".format(g))
        if startup_fuel is None and startup_cost is None:
            return [value(m.MinimumDownTime[g])] 
        elif startup_cost is None:
            return _get_startup_lag(startup_fuel, value(m.MinimumDownTime[g]))
        else:
            return _get_startup_lag(startup_cost, value(m.MinimumDownTime[g]))
    model.StartupLags = Set(model.ThermalGenerators, within=NonNegativeReals, ordered=True, initialize=startup_lags_init_rule) # units are hours / time periods.

    def _get_startup_cost(startup, fixed_adder, multiplier):
        try:
            iter(startup)
        except TypeError:
            return [fixed_adder+multiplier*startup]
        else:
            return [fixed_adder+multiplier*i[1] for i in startup]
    
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
            return _get_startup_cost(startup_fuel, fixed_startup_cost, fuel_cost)
        else:
            return _get_startup_cost(startup_cost, fixed_startup_cost, 1.)
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
       startup_costs = m.StartupCosts[g]
       for i in range(1, len(startup_costs)-1):
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

    ## FUEL-SUPPLY Sets

    def fuel_supply_gens_init(m):
        if 'fuel_supply' not in elements and ('fuel_supply' in thermal_gen_attrs or 'aux_fuel_supply' in thermal_gen_attrs):
            logger.warning('WARNING: Some generators have \'fuel_supply\' marked, but no fuel supply was found on ModelData.data[\'system\']')
            return iter(())
        if 'fuel_supply' in elements and ('fuel_supply' not in thermal_gen_attrs and 'aux_fuel_supply' not in thermal_gen_attrs):
            logger.warning('WARNING: fuel_supply in ModelData.data["elements"], but no generators are attached to any fuel supply')
            return iter(())
        if 'fuel_supply' not in thermal_gen_attrs:
            thermal_gen_attrs['fuel_supply'] = dict()
        if 'aux_fuel_supply' not in thermal_gen_attrs:
            thermal_gen_attrs['aux_fuel_supply'] = dict()
        fuel_supply = thermal_gen_attrs['fuel_supply']
        for g in fuel_supply:
            yield g
        for g in thermal_gen_attrs['aux_fuel_supply']:
            if g not in fuel_supply:
                yield g

    def gen_cost_fuel_validator(m,g):
        # validators may get called once 
        # with key None for empty sets
        if g is None:
            return True
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
    
    ## BEGIN PRODUCTION COST
    ## NOTE: For better or worse, we handle scaling this to the time period length in the objective function.
    ##       In particular, this is done in objective.py.
    
    def _check_curve(m, g, curve, curve_type):

        for i, t in enumerate(m.TimePeriods):
            ## first, get a cost_curve out of time series
            if curve['data_type'] == 'time_series':
                curve_t = curve['values'][i]
            else:
                curve_t = curve 

            ## validate that what we have is a cost_curve
            if curve_t['data_type'] != curve_type:
                raise Exception("p_cost must be of data_type cost_curve.")

            ## get the values, check against something empty
            values = curve_t['values']
            if len(values) == 0:
                if curve_t == curve:
                    logger.warning("WARNING: Generator {} has no cost information associated with it".format(g))
                    return True
                else:
                    logger.warning("WARNING: Generator {} has no cost information associated with it at time {}".format(g,t))

            ## if we have a piecewise cost curve, ensure its convexity past p_min
            ## if no curve_type+'_type' is specified, we assume piecewise (for backwards 
            ## compatibility with no 'fuel_curve_type')
            if curve_type+'_type' not in curve_t or \
                    curve_t[curve_type+'_type'] == 'piecewise':
                p_min = value(m.MinimumPowerOutput[g,t])
                last_slope = None
                for (o1, c1), (o2, c2) in zip(values, values[1:]):
                    if o2 <= p_min or math.isclose(p_min, o2):
                        continue
                    if math.isclose(o2,o1):
                        if math.isclose(c2,c1):
                            continue
                        raise Exception("Piecewise {} must be convex above p_min. ".format(curve_type) + \
                                        "Found non-convex piecewise {} for generator {} at time {}".format(curve_type,g,t))
                    ## else p_min > o2
                    if last_slope is None:
                        last_slope = (c2-c1)/(o2-o1)
                        continue
                    this_slope = (c2-c1)/(o2-o1)
                    if this_slope < last_slope and not math.isclose(this_slope, last_slope):
                        raise Exception("Piecewise {} must be convex above p_min. ".format(curve_type) + \
                                        "Found non-convex piecewise {} for generator {} at time {}".format(curve_type,g,t))
                ## verify the last output value is at least p_max
                o_last = values[-1][0]
                if value(m.MaximumPowerOutput[g,t]) > o_last and \
                        not math.isclose(value(m.MaximumPowerOutput[g,t]), o_last):
                    raise Exception("{} does not contain p_max for generator {} at time {}".format(curve_type,g,t))

            ## if we have a quadratic cost curve, ensure its convexity
            elif curve_t[curve_type+'_type'] == 'polynomial':
                if not _check_curve.warn_piecewise_approx:
                    logger.warning("WARNING: Polynomial cost curves will be approximated using piecewise segments")
                    _check_curve.warn_piecewise_approx = True
                values = curve_t['values']
                if set(values.keys()) <= {0,1,2}:
                    if 2 in values and values[2] < 0:
                        raise Exception("Polynomial {}s must be convex. ".format(curve_type) + \
                                        "Found non-convex {} for generator {} at time {}.".format(curve_type,g,t))
                    if curve_t == curve: ## in this case, no need to check the other time periods
                        return
                else:
                    raise Exception("Polynomial {}s must be quatric. ".format(curve_type) + \
                                    "Found non-quatric {} for generator {} at time {}.".format(curve_type,g,t))
            else:
                raise Exception("Unexpected {}_type".format(curve_type))

    ## set "static" variable for this function
    _check_curve.warn_piecewise_approx = False
    
    def validate_cost_rule(m, g):
        gen_dict = thermal_gens[g]
        cost = gen_dict.get('p_cost')
        fuel = gen_dict.get('p_fuel')
        fuel_cost = gen_dict.get('fuel_cost')

        if cost is None and fuel is None:
            logger.warning("WARNING: Generator {} has no cost information associated with it".format(g))
            return True
        if cost is not None and fuel is not None:
            logger.warning("WARNING: ignoring provided p_cost and using fuel cost data from p_fuel for generator {}".format(g))

        ## look at p_cost through time
        if fuel is None:
            _check_curve(m, g, cost, 'cost_curve')
        else:
            if fuel_cost is None:
                raise Exception("Found fuel_curve but not fuel_cost for generator {}".format(g))
            _check_curve(m, g, fuel, 'fuel_curve')
            for i, t in enumerate(m.TimePeriods):
                if fuel_cost is dict:
                    if fuel_cost['data_type'] != 'time_series':
                        raise Exception("fuel_cost must be either numeric or time_series")
                    fuel_cost_t = fuel_cost['values'][i]
                else:
                    fuel_cost_t = fuel_cost
                if fuel_cost_t < 0:
                    raise Exception("fuel_cost must be non-negative, found negative fuel_cost for generator {}".format(g))
                if fuel_cost_t == fuel_cost:
                    break
        return True

    model.ValidateGeneratorCost = BuildCheck(model.ThermalGenerators, rule=validate_cost_rule)
    

    ##############################################################################################
    # number of pieces in the linearization of each generator's quadratic cost production curve. #
    ##############################################################################################
    ## TODO: option-drive with Egret, either globally or per-generator
    
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
    model.PowerGenerationPiecewiseCostValues = {}

    # NOTE; these values are relative to the minimum fuel conumption
    model.PowerGenerationPiecewiseFuelValues = {}
    
    _minimum_production_cost = {}
    _minimum_fuel_consumption = {}


    def _eliminate_piecewise_duplicates(input_func):
        if len(input_func) <= 1:
            return input_func
        new = [input_func[0]]
        for (o1, c1), (o2, c2) in zip(input_func, input_func[1:]):
            if not math.isclose(o1,o2) and not math.isclose(c1,c2):
                new.append((o2,c2))
        return new

    def _much_less_than(v1, v2):
        return v1 < v2 and not math.isclose(v1,v2)

    def _piecewise_adjustment_helper(m, p_min, p_max, input_func):

        minimum_val = 0.
        new_points = []
        new_vals = []

        input_func = _eliminate_piecewise_duplicates(input_func)

        set_p_min = False

        # NOTE: this implicitly inserts a (0.,0.)
        #       into every cost array
        prior_output, prior_cost = 0., 0. 

        for output, cost in input_func:
            ## catch this case
            if math.isclose(output, p_min) and math.isclose(output, p_max):
                new_points.append(0.)
                new_vals.append(0.)
                minimum_val = cost
                break

            ## output < p_min
            elif _much_less_than(output, p_min):
                pass

            ## p_min == output
            elif math.isclose(output, p_min):
                assert set_p_min is False
                new_points.append(0.)
                new_vals.append(0.)
                minimum_val = cost
                set_p_min = True

            ## p_min < output
            elif _much_less_than(p_min, output) and _much_less_than(output, p_max):
                if not set_p_min:
                    new_points.append(0.)
                    new_vals.append(0.)

                    price = ((cost-prior_cost)/(output-prior_output))
                    minimum_val = (p_min - prior_output) * price + prior_cost
                    
                    new_points.append( output - p_min )
                    new_vals.append( (output - p_min) * price )

                    set_p_min = True
                else:
                    new_points.append( output - p_min )
                    new_vals.append( cost - minimum_val )

            elif math.isclose(output, p_max) or _much_less_than(p_max, output):
                if not set_p_min:
                    new_points.append(0.)
                    new_vals.append(0.)

                    price = ((cost-prior_cost)/(output-prior_output))
                    minimum_val = (p_min - prior_output) * price + prior_cost
                    
                    new_points.append( p_max - p_min )

                    if math.isclose(output, p_max):
                        new_vals.append( cost - minimum_val )
                    else:
                        new_vals.append( (p_max - p_min) * price )
                    set_p_min = True

                else:
                    new_points.append( p_max - p_min )
                    if math.isclose(output, p_max):
                        new_vals.append( cost - minimum_val )
                    else:
                        price = ((cost-prior_cost)/(output-prior_output))
                        new_vals.append( (p_max - prior_output) * price + prior_cost - minimum_val )

                break

            else:
                raise Exception("Unexpected case in _piecewise_adjustment_helper, "
                                "p_min={}, p_max={}, output={}".format(p_min, p_max, output))
            
            prior_output, prior_cost = output, cost

        return new_points, new_vals, minimum_val

    def _polynomial_to_piecewise_helper(m, p_min, p_max, input_func):
        segment_max = value(m.NumGeneratorCostCurvePieces)

        for key in {0,1,2}:
            if key not in input_func:
                input_func[key] = 0.

        poly_func = lambda x : input_func[0] + input_func[1]*x + input_func[2]*x*x

        if p_min >= p_max:
            minimum_val = poly_func(p_min)
            new_points = [0.]
            new_vals = [0.]
            return new_points, new_vals, minimum_val

        elif input_func[2] == 0.: ## not actually quadratic 
            minimum_val = poly_func(p_min)
            new_points = [0., p_max - p_min]
            new_vals = [0., poly_func(p_max) - minimum_val]
            return new_points, new_vals, minimum_val

        ## actually quadratic
        width = (p_max - p_min)/float(segment_max)

        new_points = [i*width for i in range(0, segment_max+1)]

        ## replace the last with (p_max - p_min)
        new_points[-1] = p_max - p_min

        minimum_val = poly_func(p_min)
        new_vals = [ poly_func(pnt+p_min) - minimum_val for pnt in new_points ]

        return new_points, new_vals, minimum_val

    def _piecewise_helper(m, p_min, p_max, curve, curve_type):
        if curve_type not in curve or \
                curve[curve_type] == 'piecewise':
            return _piecewise_adjustment_helper(m, p_min, p_max, curve['values']) 
        else:
            assert curve[curve_type] == 'polynomial'
            return _polynomial_to_piecewise_helper(m, p_min, p_max, curve['values']) 

    
    def power_generation_piecewise_points_rule(m, g):

        ## NOTE: it is often (usually) the case that cost curves
        ##       are the same in every time period, This function
        ##       is optimized to avoid data redunancy and recalculation
        ##       for that case

        gen_dict = thermal_gens[g]

        fuel_curve = gen_dict.get('p_fuel')
        cost_curve = gen_dict.get('p_cost')
        fuel_cost = gen_dict.get('fuel_cost', 0.)
        no_load_cost = gen_dict.get('non_fuel_no_load_cost', 0.)

        if isinstance(fuel_cost,dict):
            fuel_costs = fuel_cost['values']
        else:
            fuel_costs = ( fuel_cost for t in m.TimePeriods )
        if isinstance(no_load_cost,dict):
            no_load_costs = no_load_cost['values']
        else:
            no_load_costs = ( no_load_cost for t in m.TimePeriods )

        _curve_cache = dict()

        if fuel_curve is not None:

            g_in_fuel_supply_generators = g in m.FuelSupplyGenerators
            g_in_single_fuel_generators = g in m.SingleFuelGenerators

            if isinstance(fuel_curve,dict) and fuel_curve['data_type'] == 'time_series':
                fuel_curves = fuel_curve['values']
                one_fuel_curve = False
            else:
                fuel_curves = ( fuel_curve for t in m.TimePeriods )
                one_fuel_curve = True

            for fuel_curve, fuel_cost, nlc, t in zip(fuel_curves, fuel_costs, no_load_costs, m.TimePeriods):
                p_min = value(m.MinimumPowerOutput[g,t])
                p_max = value(m.MaximumPowerOutput[g,t])

                if (p_min, p_max, fuel_cost, nlc) in _curve_cache:
                    curve = _curve_cache[p_min, p_max, fuel_cost, nlc]
                    if one_fuel_curve or curve['fuel_curve'] == fuel_curve:
                        m.PowerGenerationPiecewisePoints[g,t] = curve['points']
                        if g_in_fuel_supply_generators:
                            _minimum_fuel_consumption[g,t] = curve['min_fuel_consumption']
                            m.PowerGenerationPiecewiseFuelValues[g,t] = curve['fuel_values']
                        if g_in_single_fuel_generators:
                            _minimum_production_cost[g,t] = curve['min_production_cost']
                            m.PowerGenerationPiecewiseCostValues[g,t] = curve['cost_values']
                        continue
                    
                points, values, minimum_val = _piecewise_helper(m, p_min, p_max, fuel_curve, 'fuel_curve_type')
                
                curve = { 'points' : points }

                if not one_fuel_curve:
                    curve['fuel_curve'] = fuel_curve

                m.PowerGenerationPiecewisePoints[g,t] = points
                if g_in_fuel_supply_generators:
                    _minimum_fuel_consumption[g,t] = minimum_val
                    curve['min_fuel_consumption'] = minimum_val

                    m.PowerGenerationPiecewiseFuelValues[g,t] = values
                    curve['fuel_values'] = values

                if g_in_single_fuel_generators:
                    
                    min_production_cost = minimum_val*fuel_cost + no_load_cost
                    _minimum_production_cost[g,t] = min_production_cost
                    curve['min_production_cost'] = min_production_cost

                    cost_values = [ fuel_cost*val for val in values ]
                    m.PowerGenerationPiecewiseCostValues[g,t] = cost_values
                    curve['cost_values'] = cost_values

                _curve_cache[p_min, p_max, fuel_cost, nlc] = curve

            return ## we can assume below that we don't have a fuel curve

        if isinstance(cost_curve,dict) and cost_curve['data_type'] == 'time_series':
            cost_curves = cost_curve['values']
            one_cost_curve = False
        else:
            cost_curves = ( cost_curve for t in m.TimePeriods )
            one_cost_curve = True

        for cost_curve, nlc, t in zip(cost_curves, no_load_costs, m.TimePeriods):
            p_min = value(m.MinimumPowerOutput[g,t])
            p_max = value(m.MaximumPowerOutput[g,t])

            if (p_min, p_max, nlc) in _curve_cache:
                curve = _curve_cache[p_min, p_max, nlc]
                if one_cost_curve or curve['cost_curve'] == cost_curve:
                    m.PowerGenerationPiecewisePoints[g,t] = curve['points']
                    m.PowerGenerationPiecewiseCostValues[g,t] = curve['cost_values']
                    _minimum_production_cost[g,t] = curve['min_production']
                    continue

            if cost_curve is None:
                if p_min >= p_max: ## only one point
                    points = [0.]
                    values = [0.]
                else:
                    points = [0., p_max - p_min]
                    values = [0., 0.]
                min_production = nlc
            else:
                points, values, minimum_val = _piecewise_helper(m, p_min, p_max, cost_curve, 'cost_curve_type')
                min_production = minimum_val + nlc
    
            curve = {'points':points, 'cost_values':values, 'min_production':min_production}
            if not one_cost_curve:
                curve['cost_curve'] = cost_curve
            _curve_cache[p_min, p_max, nlc] = curve

            m.PowerGenerationPiecewisePoints[g,t] = points
            m.PowerGenerationPiecewiseCostValues[g,t] = values
            _minimum_production_cost[g,t] = min_production 

    model.CreatePowerGenerationPiecewisePoints = BuildAction(model.ThermalGenerators, rule=power_generation_piecewise_points_rule)

    # Minimum production cost (needed because Piecewise constraint on ProductionCost 
    # has to have lower bound of 0, so the unit can cost 0 when off -- this is added
    # back in to the objective if a unit is on

    model.MinimumProductionCost = Param(model.SingleFuelGenerators, model.TimePeriods, within=NonNegativeReals, initialize=_minimum_production_cost, mutable=True)

    model.MinimumFuelConsumption = Param(model.FuelSupplyGenerators, model.TimePeriods, within=NonNegativeReals, initialize=_minimum_fuel_consumption, mutable=True)

    ## END PRODUCTION COST CALCULATIONS

    #########################################
    # penalty costs for constraint violation #
    #########################################

    ModeratelyBigPenalty = 1e3*system['baseMVA']
    
    model.ReserveShortfallPenalty = Param(within=NonNegativeReals, default=ModeratelyBigPenalty, mutable=True, initialize=system.get('reserve_shortfall_cost', ModeratelyBigPenalty))
    
    BigPenalty = 1e4*system['baseMVA']
    
    model.LoadMismatchPenalty = Param(within=NonNegativeReals, mutable=True, initialize=system.get('load_mismatch_cost', BigPenalty))
    model.LoadMismatchPenaltyReactive = Param(within=NonNegativeReals, mutable=True, initialize=system.get('q_load_mismatch_cost', BigPenalty/2.))

    model.Contingencies = Set(initialize=contingencies.keys())

    # leaving this unindexed for now for simpility
    model.ContingencyLimitPenalty = Param(within=NonNegativeReals, initialize=system.get('contingency_flow_violation_cost', BigPenalty/2.), mutable=True)

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
                                            default=0.0, initialize=storage_attrs.get('min_discharge_rate', dict()))
    
    def maximum_power_output_validator_storage(m, v, s):
        return v >= value(m.MinimumPowerOutputStorage[s])
    
    model.MaximumPowerOutputStorage = Param(model.Storage, within=NonNegativeReals,
                                            validate=maximum_power_output_validator_storage, default=0.0,
                                            initialize=storage_attrs.get('max_discharge_rate', dict()))
    
    #Storage power input >0 when charging
    
    model.MinimumPowerInputStorage = Param(model.Storage, within=NonNegativeReals,
                                            default=0.0, initialize=storage_attrs.get('min_charge_rate', dict()))
    
    def maximum_power_input_validator_storage(m, v, s):
        return v >= value(m.MinimumPowerInputStorage[s])
    
    model.MaximumPowerInputStorage = Param(model.Storage, within=NonNegativeReals,
                                            validate=maximum_power_input_validator_storage, default=0.0,
                                            initialize=storage_attrs.get('max_charge_rate', dict()))
    
    ###############################################
    # storage ramp up/down rates. units are MW/h. #
    ###############################################
    
    # ramp rate limits when discharging
    model.NominalRampUpLimitStorageOutput    = Param(model.Storage, within=NonNegativeReals,
                                                        initialize=storage_attrs.get('ramp_up_output_60min', dict()))
    model.NominalRampDownLimitStorageOutput  = Param(model.Storage, within=NonNegativeReals,
                                                        initialize=storage_attrs.get('ramp_down_output_60min', dict()))
    
    # ramp rate limits when charging
    model.NominalRampUpLimitStorageInput     = Param(model.Storage, within=NonNegativeReals,
                                                        initialize=storage_attrs.get('ramp_up_input_60min', dict()))
    model.NominalRampDownLimitStorageInput   = Param(model.Storage, within=NonNegativeReals,
                                                        initialize=storage_attrs.get('ramp_down_input_60min', dict()))
    
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
                                        initialize=storage_attrs.get('energy_capacity', dict()))
    model.MinimumSocStorage = Param(model.Storage, within=PercentFraction, default=0.0,
                                        initialize=storage_attrs.get('minimum_state_of_charge', dict()))
    
    ################################################################################
    # round trip efficiency for each storage unit given as a fraction (i.e. [0,1]) #
    ################################################################################
    
    model.InputEfficiencyEnergy  = Param(model.Storage, within=PercentFraction, default=1.0,
                                            initialize=storage_attrs.get('charge_efficiency', dict()))
    model.OutputEfficiencyEnergy = Param(model.Storage, within=PercentFraction, default=1.0,
                                            initialize=storage_attrs.get('discharge_efficienty', dict()))
    model.RetentionRate          = Param(model.Storage, within=PercentFraction, default=1.0,
                                            initialize=storage_attrs.get('retention_rate_60min', dict())) ## assumed to be %/hr

    model.ChargeCost = Param(model.Storage, within=Reals, default=0.0, initialize=storage_attrs.get('charge_cost', dict()))
    model.DischargeCost = Param(model.Storage, within=Reals, default=0.0, initialize=storage_attrs.get('discharge_cost', dict()))

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
    def _end_point_soc(m, s):
        if s is None:
            return
        s_dict = storage[s]
        if 'end_state_of_charge' in s_dict:
            return s_dict['end_state_of_charge']
        if 'initial_state_of_charge' in s_dict:
            return s_dict['initial_state_of_charge']
        return 0.5
    
    model.EndPointSocStorage = Param(model.Storage, within=PercentFraction, default=0.5,
                                        initialize=_end_point_soc)
    
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
                                            initialize=storage_attrs.get('initial_discharge_rate', dict()))
    model.StoragePowerInputOnT0  = Param(model.Storage, within=NonNegativeReals,
                                            validate=t0_storage_power_input_validator,
                                            default=0.0,
                                            initialize=storage_attrs.get('initial_charge_rate', dict()))
    model.StorageSocOnT0         = Param(model.Storage, within=PercentFraction,
                                            default=0.5, initialize=storage_attrs.get('initial_state_of_charge', dict()))
    return model
