#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module contains several helper functions that are useful when
working with transmission models
"""

import egret.model_library.transmission.tx_calc as tx_calc
import math
from math import radians
import logging
import copy
from collections import OrderedDict
from egret.data.data_utils import map_items, zip_items


logger = logging.getLogger(__name__)


def get_unique_bus_pairs(md):
    branch_attrs = md.attributes(element_type='branch')
    bus_pairs = zip_items(branch_attrs['from_bus'], branch_attrs['to_bus'])
    unique_bus_pairs = list(OrderedDict((val, None) for idx, val in bus_pairs.items()))
    return unique_bus_pairs


def _get_out_of_service_gens(md):
    out_of_service_gens = list()
    for gen_name, gen_dict in md.elements(element_type='generator'):
        if not gen_dict['in_service']:
            out_of_service_gens.append(gen_name)
            gen_dict['in_service'] = True

    return out_of_service_gens


def _get_out_of_service_branches(md):
    out_of_service_branches = list()
    for branch_name, branch_dict in md.elements(element_type='branch'):
        if not branch_dict['in_service']:
            out_of_service_branches.append(branch_name)
            branch_dict['in_service'] = True

    return out_of_service_branches


def dicts_of_vr_vj(buses):
    """
    Create dictionaries of vr and vj values from the bus vm and va values
    """
    # TODO: Change api to be vr_vj_dicts_from_vm_va(bus_vm, bus_va)
    vr = dict()
    vj = dict()
    for bus_name, bus in buses.items():
        vr[bus_name] = tx_calc.calculate_vr_from_vm_va_degrees(bus['vm'], bus['va'])
        vj[bus_name] = tx_calc.calculate_vj_from_vm_va_degrees(bus['vm'], bus['va'])

    return vr, vj


def dict_of_bus_loads(buses, loads):
    """
    Create dictionaries of the p and q bus load values from the
    load elements
    """
    # loop over all the load elements and sum the loads at each of the buses
    # TODO: Make this dictionary so that it returns None when no load
    bus_p_loads = {k: 0 for k in buses.keys()}
    bus_q_loads = {k: 0 for k in buses.keys()}

    for load_name, load in loads.items():
        bus_name = load['bus']
        ## NOTE: for DC models we may not have q_load defined
        ##       making this the same for p_load too..?
        if 'p_load' in load:
            bus_p_loads[bus_name] += load['p_load']
        if 'q_load' in load:
            bus_q_loads[bus_name] += load['q_load']

    return bus_p_loads, bus_q_loads


def dict_of_bus_fixed_shunts(buses, shunts):
    """
    Create dictionaries of the p and q bus shunt values from the
    shunt elements
    """
    # loop over all the load elements and sum the loads at each of the buses
    # TODO: Make this dictionary so that it returns None when no shunt
    bus_bs_fixed_shunts = {k: 0 for k in buses.keys()}
    bus_gs_fixed_shunts = {k: 0 for k in buses.keys()}

    for shunt_name, shunt in shunts.items():
        if shunt['shunt_type'] == 'fixed':
            bus_name = shunt['bus']
            if shunt['bs'] != 0.0:
                bus_bs_fixed_shunts[bus_name] += shunt['bs']
            if shunt['gs'] != 0.0:
                bus_gs_fixed_shunts[bus_name] += shunt['gs']

    return bus_bs_fixed_shunts, bus_gs_fixed_shunts

def dict_of_branch_currents(branches, buses):
    """
    Create a dictionary of the branch currents
    (with subkeys ifr, ifj, itr, itj)
    """
    branch_currents = dict()
    branch_currents['ifr'] = dict()
    branch_currents['ifj'] = dict()
    branch_currents['itr'] = dict()
    branch_currents['itj'] = dict()
    for branch_name, branch in branches.items():
        from_bus = buses[branch['from_bus']]
        to_bus = buses[branch['to_bus']]
        ifr = 0
        ifj = 0
        itr = 0
        itj = 0
        if branch['in_service'] \
            and from_bus['vm'] is not None and from_bus['va'] is not None \
            and to_bus['vm'] is not None and to_bus['va'] is not None:
            # we have all the information we need
            y_matrix = tx_calc.calculate_y_matrix_from_branch(branch)
            vfr = tx_calc.calculate_vr_from_vm_va_degrees(from_bus['vm'], from_bus['va'])
            vfj = tx_calc.calculate_vj_from_vm_va_degrees(from_bus['vm'], from_bus['va'])
            vtr = tx_calc.calculate_vr_from_vm_va_degrees(to_bus['vm'], to_bus['va'])
            vtj = tx_calc.calculate_vj_from_vm_va_degrees(to_bus['vm'], to_bus['va'])
            ifr = tx_calc.calculate_ifr(vfr, vfj, vtr, vtj, y_matrix)
            ifj = tx_calc.calculate_ifj(vfr, vfj, vtr, vtj, y_matrix)
            itr = tx_calc.calculate_itr(vfr, vfj, vtr, vtj, y_matrix)
            itj = tx_calc.calculate_itj(vfr, vfj, vtr, vtj, y_matrix)
        branch_currents['ifr'][branch_name] = ifr
        branch_currents['ifj'][branch_name] = ifj
        branch_currents['itr'][branch_name] = itr
        branch_currents['itj'][branch_name] = itj
    return branch_currents


def dict_of_branch_powers(branches, buses):
    """
    Create a dictionary of the branch powers
    (with subkeys pf, qf, pt, qt)
    """
    branch_powers = dict()
    branch_powers['pf'] = dict()
    branch_powers['qf'] = dict()
    branch_powers['pt'] = dict()
    branch_powers['qt'] = dict()
    for branch_name, branch in branches.items():
        from_bus = buses[branch['from_bus']]
        to_bus = buses[branch['to_bus']]
        pf = 0
        qf = 0
        pt = 0
        qt = 0
        if branch['in_service'] \
            and from_bus['vm'] is not None and from_bus['va'] is not None \
            and to_bus['vm'] is not None and to_bus['va'] is not None:
            # we have all the information we need
            y_matrix = tx_calc.calculate_y_matrix_from_branch(branch)
            vfr = tx_calc.calculate_vr_from_vm_va_degrees(from_bus['vm'], from_bus['va'])
            vfj = tx_calc.calculate_vj_from_vm_va_degrees(from_bus['vm'], from_bus['va'])
            vtr = tx_calc.calculate_vr_from_vm_va_degrees(to_bus['vm'], to_bus['va'])
            vtj = tx_calc.calculate_vj_from_vm_va_degrees(to_bus['vm'], to_bus['va'])
            ifr = tx_calc.calculate_ifr(vfr, vfj, vtr, vtj, y_matrix)
            ifj = tx_calc.calculate_ifj(vfr, vfj, vtr, vtj, y_matrix)
            itr = tx_calc.calculate_itr(vfr, vfj, vtr, vtj, y_matrix)
            itj = tx_calc.calculate_itj(vfr, vfj, vtr, vtj, y_matrix)
            pf = tx_calc.calculate_p(ifr, ifj, vfr, vfj)
            qf = tx_calc.calculate_q(ifr, ifj, vfr, vfj)
            pt = tx_calc.calculate_p(itr, itj, vtr, vtj)
            qt = tx_calc.calculate_q(itr, itj, vtr, vtj)
        branch_powers['pf'][branch_name] = pf
        branch_powers['qf'][branch_name] = qf
        branch_powers['pt'][branch_name] = pt
        branch_powers['qt'][branch_name] = qt
    return branch_powers


def inlet_outlet_branches_by_bus(branches, buses):
    """
    Return dictionaries of the inlet and outlet branches
    to each bus
    """
    inlet_branches_by_bus = {k: list() for k in buses.keys()}
    outlet_branches_by_bus ={k: list() for k in buses.keys()}

    for branch_name, branch in branches.items():
        inlet_branches_by_bus[branch['to_bus']].append(branch_name)
        outlet_branches_by_bus[branch['from_bus']].append(branch_name)

    return inlet_branches_by_bus, outlet_branches_by_bus


def gens_by_bus(buses, gens):
    """
    Return a dictionary of the generators attached to each bus
    """
    gens_by_bus = {k: list() for k in buses.keys()}
    for gen_name, gen in gens.items():
        gens_by_bus[gen['bus']].append(gen_name)

    return gens_by_bus

def over_gen_limit(load, gens, gen_maxs):
    '''
    Calculates the maximum amount of over-generation
    given a load and set of generators with
    associated maximum outputs
    '''
    max_over_gen = 0.
    if load < 0.:
        max_over_gen += -load
    for g in gens:
        g_max = gen_maxs[g]
        if g_max > 0:
            max_over_gen += g_max

    return max_over_gen

def load_shed_limit(load, gens, gen_mins):
    '''
    Calculates the maximum amount of load shedding
    given a load and set of generators with
    associated minimum outputs
    '''
    max_load_shed = 0.
    if load > 0.:
        max_load_shed += load
    for g in gens:
        g_min = gen_mins[g]
        if g_min < 0:
            max_load_shed += -g_min

    return max_load_shed

## attributes which are scaled for power flow models
## tuple of supported ancillary services, as named in model data
ancillary_services = (
    'reserve',
    'spinning_reserve',
    'non_spinning_reserve',
    'regulation_up',
    'regulation_down',
    'supplemental_reserve',
    'flexible_ramp_up',
    'flexible_ramp_down',
)

## tuple of penalty prices (up/down reserves have the same penalty price)
penalty_prices = (
    'regulation_penalty_price',
    'spinning_reserve_penalty_price',
    'non_spinning_reserve_penalty_price',
    'supplemental_reserve_penalty_price',
    'flexible_ramp_penalty_price',
)

# construct this from the above items
ancillary_service_stack = [
        *(name+'_requirement' for name in ancillary_services),
        *(name+'_shortfall' for name in ancillary_services),
        *(name+'_price' for name in ancillary_services),
        *penalty_prices,
]


## TODO?: break apart by data that needed to be scaled down (capacity limits, power),
## vs. scaled up (costs, prices, etc)
scaled_attributes = {
                         ('element_type','generator', None): [
                                                          'p_min',
                                                          'p_max',
                                                          'p_min_agc',
                                                          'p_max_agc',
                                                          'q_min',
                                                          'q_max',
                                                          'startup_capacity',
                                                          'shutdown_capacity',
                                                          'startup_curve',
                                                          'shutdown_curve',
                                                          'ramp_up_60min',
                                                          'ramp_down_60min',
                                                          'initial_p_output',
                                                          'initial_q_output',
                                                          'pc1',
                                                          'pc2',
                                                          'qc1_min',
                                                          'qc1_max',
                                                          'qc2_min',
                                                          'qc2_max',
                                                          'ramp_agc',
                                                          'ramp_10',
                                                          'ramp_30',
                                                          'ramp_q',
                                                          'pg',
                                                          'qg',
                                                          'reserve_supplied',
                                                          'headroom',
                                                          'regulation_up_supplied',
                                                          'regulation_down_supplied',
                                                          'flexible_ramp_up_supplied',
                                                          'flexible_ramp_down_supplied',
                                                          'spinning_reserve_supplied',
                                                          'non_spinning_reserve_supplied',
                                                          'supplemental_reserve_supplied',
                                                          'p_cost',
                                                          'p_fuel',
                                                          'q_cost',
                                                          'agc_marginal_cost',
                                                          'spinning_cost',
                                                          'non_spinning_cost',
                                                          'supplemental_cost',
                                                          'spinning_capacity',
                                                          'non_spinning_capacity',
                                                          'supplemental_spinning_capacity',
                                                          'supplemental_non_spinning_capacity',
                                                       ],
                       ('element_type','storage', None): [
                                                        'energy_capacity',
                                                        'max_discharge_rate',
                                                        'min_discharge_rate',
                                                        'max_charge_rate',
                                                        'min_charge_rate',
                                                        'ramp_up_input_60min',
                                                        'ramp_down_input_60min',
                                                        'ramp_up_output_60min',
                                                        'ramp_down_output_60min',
                                                        'p_discharge',
                                                        'p_charge',
                                                        'charge_cost',
                                                        'discharge_cost',
                                                   ],
                       ('element_type','load', None) : [
                                                      'p_load',
                                                      'q_load',
                                                      'p_load_shed',
                                                      'q_load_shed',
                                                      'p_price',
                                                      'q_price',
                                                     ],
                       ('element_type','branch', None) : [
                                                       'rating_long_term',
                                                        'rating_short_term',
                                                        'rating_emergency',
                                                        'pf',
                                                        'qf',
                                                        'pt',
                                                        'qt',
                                                        'violation_penalty',
                                                        'pf_violation',
                                                        ],
                       ('element_type','dc_branch', None) : [
                                                        'rating_long_term',
                                                        'rating_short_term',
                                                        'rating_emergency',
                                                        'pf',
                                                        'qf',
                                                        'pt',
                                                        'qt',
                                                        'violation_penalty',
                                                        'pf_violation',
                                                        ],
                       ('element_type', 'shunt', None) : [
                                                      'bs',
                                                      'gs',
                                                      'bs_min',
                                                      'bs_max',
                                                      'gs_min',
                                                      'gs_max',
                                                     ],
                       ('element_type', 'area', None) : [
                                                      ] + \
                                                  ancillary_service_stack,
                       ('element_type', 'zone', None) : [
                                                      ] + \
                                                  ancillary_service_stack,
                       ('element_type', 'interface', None) : [ 
                                                         'minimum_limit',
                                                         'maximum_limit',
                                                         'pf',
                                                         'qf',
                                                         'pt',
                                                         'qt',
                                                         'violation_penalty',
                                                         'pf_violation',
                                                       ],
                       ('element_type', 'bus', None) : [ 
                                                    'p_balance_violation',
                                                    'q_balance_violation',
                                                    'lmp',
                                                    'q_lmp',
                                                    'qlmp',
                                                    'pl',
                                                    'ql',
                                                 ],
                       ('element_type', 'security_constraint', None) : [
                                                                       ],
                       ('element_type', 'security_constraint', 'pg') : [ 'lower_bound',
                                                                         'upper_bound',
                                                                         'violation_penalty',
                                                                         'pf',
                                                                         'pf_violation',
                                                                       ],
                       ('element_type', 'contingency', None) : [ 
                                                                 'violation_penalty',
                                                                 'pf',
                                                                 'pf_violation',
                                                               ],
                       ('element_type', 'fuel_supply', None) : [
                                                               ],
                       ('system_attributes', None, None ) : [
                                                        'load_mismatch_cost',
                                                        'q_load_mismatch_cost',
                                                        'transmission_flow_violation_cost',
                                                        'contingency_flow_violation_cost',
                                                        'interface_flow_violation_cost',
                                                        'reserve_shortfall_cost',
                                                     ] + \
                                                     ancillary_service_stack,
                   }

def element_types():
    ''' Get an iterable that yields each valid egret element type as a string
    '''
    return (key[1] for key in scaled_attributes.keys()
            if key[0] == 'element_type' and key[2] is None)

def scale_ModelData_to_pu(model_data, inplace=False):
    return _convert_modeldata_pu(model_data, _divide_by_baseMVA, inplace)


def unscale_ModelData_to_pu(model_data, inplace=False):
    return _convert_modeldata_pu(model_data, _multiply_by_baseMVA, inplace)


def _multiply_by_baseMVA(element, attr_name, attr, baseMVA, attributes):
    _scale_by_baseMVA(_mul, _div, element, attr_name, attr, baseMVA, attributes)
def _divide_by_baseMVA(element, attr_name, attr, baseMVA, attributes):
    _scale_by_baseMVA(_div, _mul, element, attr_name, attr, baseMVA, attributes)

def _mul(a,b):
    return a*b
def _div(a,b):
    return a/b

def _get_op(normal_op, inverse_op, attr_name):
    if ('cost' in attr_name) or ('price' in attr_name) or ('lmp' in attr_name) or ('penalty' in attr_name):
        return inverse_op 
    return normal_op

def _no_op(a,b):
    return a

def _recurse_deeper_dict(normal_op, inverse_op, element, attr_name, value, baseMVA, attributes):
    if 'data_type' in value:
        _scale_by_baseMVA(normal_op, inverse_op, element, attr_name, value, baseMVA, attributes)
    else: # recurse deeper
        for k,v in value.items():
            _scale_by_baseMVA(normal_op, inverse_op, value, k, v, baseMVA, attributes)

def _scale_by_baseMVA(normal_op, inverse_op, element, attr_name, attr, baseMVA, attributes):
    if attr is None:
        return
    if isinstance(attr, dict):
        if 'data_type' in attr:
            if attr['data_type'] == 'time_series':
                if attr_name in attributes:
                    op = _get_op(normal_op, inverse_op, attr_name)
                else:
                    op = _no_op
                values_list = attr['values']
                for time, value in enumerate(values_list):
                    if isinstance(value, dict):
                        _recurse_deeper_dict(normal_op, inverse_op, element, attr_name, value, baseMVA, attributes)
                    elif isinstance(value, list):
                        values_list[time] = [ op(v, baseMVA) for v in value ]
                    elif isinstance(value, tuple):
                        values_list[time] = tuple( op(v, baseMVA) for v in value )
                    else:
                        values_list[time] = op( value , baseMVA )
            elif attr['data_type'] == 'cost_curve':
                if attr['cost_curve_type'] == 'polynomial':
                    values_dict = attr['values']
                    if 'data_type' in values_dict:
                        _recurse_deeper_dict(normal_op, inverse_op, element, attr_name, values_dict, baseMVA, attributes)
                    else:
                        attr['values'] = { int(power): coeff*(inverse_op(1.,baseMVA)**int(power)) \
                                        for (power, coeff) in values_dict.items() }
                elif attr['cost_curve_type'] == 'piecewise':
                    values = attr['values']
                    if isinstance(values, list):
                        attr['values'] = [ (normal_op(point,baseMVA), cost) \
                                        for (point, cost) in values ]
                    elif isinstance(values, tuple):
                        attr['values'] = tuple( (normal_op(point,baseMVA), cost) \
                                            for (point, cost) in values )
                    elif isinstance(values, dict):
                        _recurse_deeper_dict(normal_op, inverse_op, element, attr_name, values, baseMVA, attributes)
                    else:
                        raise RuntimeError("Unexpected case converting piecewise cost curve")
            elif attr['data_type'] == 'fuel_curve':
                values = attr['values']
                if isinstance(values, list):
                    attr['values'] = [ (normal_op(point,baseMVA), fuel) \
                                    for (point, fuel) in values ]
                elif isinstance(values, tuple):
                    attr['values'] = tuple( (normal_op(point,baseMVA), fuel) \
                                        for (point, fuel) in values )
                elif isinstance(values, dict):
                    _recurse_deeper_dict(normal_op, inverse_op, element, attr_name, values, baseMVA, attributes)
                else:
                    raise RuntimeError("Unexpected case converting piecewise fuel curve")
            else: # potentially recurse deeper on the "values"
                if attr_name in attributes:
                    op = _get_op(normal_op, inverse_op, attr_name)
                else:
                    op = _no_op
                values = attr['values']
                if isinstance(values, dict):
                    _recurse_deeper_dict(normal_op, inverse_op, element, attr_name, values, baseMVA, attributes)
                elif isinstance(values, list):
                    attr['values'] = [ op(v, baseMVA) for v in values ]
                elif isinstance(value, tuple):
                    attr['values'] = tuple( op(v, baseMVA) for v in values )
                else:
                    attr['values'] = op( values , baseMVA )
        else: # recurse deeper AND we've already checked for data_type
            for k,v in attr.items():
                _scale_by_baseMVA(normal_op, inverse_op, attr, k, v, baseMVA, attributes)
    elif attr_name in attributes:
        op = _get_op(normal_op, inverse_op, attr_name)
        if isinstance(attr, list):
            element[attr_name] = [ op(a, baseMVA) for a in attr ]
        elif isinstance(attr, tuple):
            element[attr_name] = tuple( op(a, baseMVA) for a in attr )
        else:
            element[attr_name] = op( attr , baseMVA )
    else:
        return


## NOTE: ideally this would be done in the definitions of
##       these constraints. Futher, it is not obvious that
##       the baseMVA provides the best scaling
## NOTE: specifying inplace returns None
def _convert_modeldata_pu(model_data, transform_func, inplace):

    if inplace:
        md = model_data
    else:
        md = model_data.clone()
    baseMVA = float(md.data['system']['baseMVA'])

    for (attr_type, element_type, element_subtype), attributes in scaled_attributes.items():

        if attr_type == 'system_attributes':
            system_dict = md.data['system']
            assert element_type is None
            assert element_subtype is None
            for name, sys_attr in system_dict.items():
                transform_func(system_dict, name, sys_attr, baseMVA, attributes)
        
        elif attr_type == 'element_type':
            if element_type not in md.data['elements']:
                continue
            element_dict = md.data['elements'][element_type]
            if element_subtype is None:
                for name, element in element_dict.items():
                    for attr_name, attr in element.items():
                        transform_func(element, attr_name, attr, baseMVA, attributes)
            else: ## allow for different actions depending on the subtype
                for name, element in element_dict.items():
                    element_subtype_key = element_type+'_type'
                    if element_subtype == element[element_subtype_key]:
                        for attr_name, attr in element.items():
                            transform_func(element, attr_name, attr, baseMVA, attributes)

    if inplace:
        return
    else:
        return md

def radians_from_degrees_dict(d):
    return {k:radians(d[k]) for k in d}

def _insert_first_point(p_min, values, pop=True):
    first_slope = (values[1][1] - values[0][1]) / (values[1][0] - values[0][0])
    first_intercept = values[0][1] - first_slope * values[0][0]
    first_cost = first_slope * p_min + first_intercept
    if pop:
        values.pop(0)
    values.insert(0, (p_min, first_cost))

def _insert_last_point(p_max, values, pop=True):
    last_slope = (values[-1][1] - values[-2][1]) / (values[-1][0] - values[-2][0])
    last_intercept = values[-1][1] - last_slope * values[-1][0]
    last_cost = last_slope * p_max + last_intercept
    if pop:
        values.pop()
    values.append((p_max, last_cost))

def _remove_duplicates(values):
    duplicate_incides = []
    for idx, ((o1, c1), (o2, c2)) in enumerate(zip(values, values[1:])):
        if math.isclose(o1,o2) and math.isclose(c1, c2):
            duplicate_incides.append(idx)
    for idx in reversed(duplicate_incides):
        del values[idx]

def validate_and_clean_cost_curve(curve, curve_type, p_min, p_max, gen_name, t=None):
    """
    Parameters
    ----------
    curve: dict
    curve_type: str
    p_min: float
    p_max: float
    gen_name: str
    t: None or int

    Returns
    -------
    cleaned_values: list or dict
        If the curve is piecewise, this function will return a list. If the curve is
        a polynomial, this function will return a dict.
    """
    # validate that what we have is a cost_curve
    if curve['data_type'] != curve_type:
        raise ValueError(f"cost curve must be of data_type {curve_type}.")

    # get the values, check against something empty
    values = curve['values']
    if len(values) == 0:
        if t is None:
            logger.warning(f"WARNING: Generator {gen_name} has no cost information associated with it")
        else:
            logger.warning(f"WARNING: Generator {gen_name} has no cost information associated with it at time {t}")
        return copy.copy(values)


    # if we have a piecewise cost curve, ensure its convexity past p_min
    # if no curve_type+'_type' is specified, we assume piecewise (for backwards
    # compatibility with no 'fuel_curve_type')
    if curve_type + '_type' not in curve or \
            curve[curve_type + '_type'] == 'piecewise':

        # handle this case specially
        if len(values) == 1:
            o1, c1 = values[0]
            # allow and resolve some FP error
            if math.isclose(p_min, p_max) and (math.isclose(p_min, o1) or math.isclose(p_max, o1)):
                return [(p_min,c1)]
            else:
                at_time_t = "" if (t is None) else f"at time {t}"
                raise ValueError(f"Generator {gen_name} {at_time_t} has only a single point on its "
                                  "piecewise cost curve which is not covered by p_min or p_max.")

        _remove_duplicates(values)

        # print a warning (once) if we're extending the cost curve
        if (p_min < values[0][0] and not math.isclose(p_min, values[0][0])) or \
                (p_max > values[-1][0] and not math.isclose(p_max, values[-1][0])):
            msg = f"WARNING: Extending piecewise linear cost curve beyond p_min and/or p_max for generator {gen_name}"
            if validate_and_clean_cost_curve._printed_warning:
                logger.debug(msg)
            else:
                logger.warning(msg+" (and perhaps others)")
                validate_and_clean_cost_curve._printed_warning = True

            # we should have copied the user's data at this point, and this
            # will allow the user to see in the output *how* we modified their
            # provided cost curve, in case that's useful
            if (p_min < values[0][0] and not math.isclose(p_min, values[0][0])):
                _insert_first_point(p_min, values, pop=False)
            if (p_max > values[-1][0] and not math.isclose(p_max, values[-1][0])):
                _insert_last_point(p_max, values, pop=False)

        # now that we've extended the cost curve, we need to catch the
        # case that p_min == p_max and we're on one of the points
        # of the cost curve. The logic below fails in this case.
        if math.isclose(p_min, p_max):
            for o, c in values:
                if math.isclose(p_min, o):
                    return [(p_min,c)]
                if math.isclose(p_max, o):
                    return [(p_max,c)]

        cleaned_values = list()
        last_slope = None
        for (o1, c1), (o2, c2) in zip(values, values[1:]):
            if o2 <= p_min or math.isclose(p_min, o2):
                continue

            if p_max <= o1 or math.isclose(p_max, o1):
                continue

            if len(cleaned_values) == 0:
                cleaned_values.append((o1, c1))

            if math.isclose(o2, o1):
                if math.isclose(c2, c1):
                    continue
                raise ValueError(f"Found piecewise {curve_type} for generator {gen_name} at time {t} " +
                                 "with nearly infinite slope.")

            cleaned_values.append((o2, c2))

            # else p_min < o2
            if last_slope is None:
                last_slope = (c2 - c1) / (o2 - o1)
                continue
            this_slope = (c2 - c1) / (o2 - o1)
            if this_slope < last_slope and not math.isclose(this_slope, last_slope):
                raise ValueError(f"Piecewise {curve_type} must be convex above p_min. " +
                                 f"Found non-convex piecewise {curve_type} for generator {gen_name} at time {t}")

            # combine pieces if slope is the same
            if math.isclose(this_slope, last_slope):
                cleaned_values.pop(-2)
            else: # update last slope
                last_slope = this_slope

        # match first and last point with p_min and p_max
        _insert_first_point(p_min, cleaned_values, pop=True)
        _insert_last_point(p_max, cleaned_values, pop=True)

        # If p_min is p_max, we'll have two identical
        # (or close to identical) points. In this case
        # we'll just take the last one.
        if math.isclose(p_min, p_max):
            cleaned_values.pop(0)

    # if we have a quadratic cost curve, ensure its convexity
    elif curve[curve_type + '_type'] == 'polynomial':
        if set(values.keys()) <= {0, 1, 2}:
            if 2 in values and values[2] < 0:
                raise ValueError(f"Polynomial {curve_type}s must be convex. " +
                                 f"Found non-convex {curve_type} for generator {gen_name} at time {t}.")
        else:
            raise ValueError(f"Polynomial {curve_type}s must be quadratic. " +
                             f"Found non-quatric {curve_type} for generator {gen_name} at time {t}.")
        cleaned_values = copy.copy(values)
    else:
        raise Exception(f"Unexpected {curve_type}_type")

    return cleaned_values

validate_and_clean_cost_curve._printed_warning = False
