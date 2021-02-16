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


def dicts_of_vr_vj(buses):
    """
    Create dictionaries of vr and vj values from the bus vm and va values
    """
    # TODO: Change api to be vr_vj_dicts_from_vm_va(bus_vm, bus_va)
    vr = dict()
    vj = dict()
    for bus_name, bus in buses.items():
        vr[bus_name] = tx_calc.calculate_vr_from_vm_va(bus['vm'], bus['va'])
        vj[bus_name] = tx_calc.calculate_vj_from_vm_va(bus['vm'], bus['va'])

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
            vfr = tx_calc.calculate_vr_from_vm_va(from_bus['vm'], from_bus['va'])
            vfj = tx_calc.calculate_vj_from_vm_va(from_bus['vm'], from_bus['va'])
            vtr = tx_calc.calculate_vr_from_vm_va(to_bus['vm'], to_bus['va'])
            vtj = tx_calc.calculate_vj_from_vm_va(to_bus['vm'], to_bus['va'])
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
            vfr = tx_calc.calculate_vr_from_vm_va(from_bus['vm'], from_bus['va'])
            vfj = tx_calc.calculate_vj_from_vm_va(from_bus['vm'], from_bus['va'])
            vtr = tx_calc.calculate_vr_from_vm_va(to_bus['vm'], to_bus['va'])
            vtj = tx_calc.calculate_vj_from_vm_va(to_bus['vm'], to_bus['va'])
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
ancillary_service_stack = [
                            'reserve_requirement',
                            'spinning_reserve_requirement',
                            'non_spinning_reserve_requirement',
                            'regulation_up_requirement',
                            'regulation_down_requirement',
                            'flexible_ramp_up_requirement',
                            'flexible_ramp_down_requirement',
                            'supplemental_reserve_requirement',
                            'reserve_shortfall',
                            'spinning_reserve_shortfall',
                            'non_spinning_reserve_shortfall',
                            'regulation_up_shortfall',
                            'regulation_down_shortfall',
                            'flexible_ramp_up_shortfall',
                            'flexible_ramp_down_shortfall',
                            'supplemental_shortfall',
                            'reserve_price',
                            'spinning_reserve_price',
                            'non_spinning_reserve_price',
                            'regulation_up_price',
                            'regulation_down_price',
                            'flexible_ramp_up_price',
                            'flexible_ramp_down_price',
                            'supplemental_price',
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
                                                          'rg',
                                                          'headroom',
                                                          'reg_up_supplied',
                                                          'reg_down_supplied',
                                                          'spin_supplied',
                                                          'flex_up_supplied',
                                                          'flex_down_supplied',
                                                          'non_spinning_supplied',
                                                          'supplemental_supplied',
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
                       ('system_attributes', None, None ) : [
                                                        'load_mismatch_cost',
                                                        'q_load_mismatch_cost',
                                                        'reserve_shortfall_cost',
                                                     ] + \
                                                     ancillary_service_stack,
                   }


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

def _scale_by_baseMVA(normal_op, inverse_op, element, attr_name, attr, baseMVA, attributes):
    if attr is None:
        return
    if isinstance(attr, dict):
        if 'data_type' in attr and attr['data_type'] == 'time_series':
            op = _get_op(normal_op, inverse_op, attr_name)
            values_list = attr['values']
            for time, value in enumerate(values_list):
                if isinstance(value, dict):
                    _scale_by_baseMVA(normal_op, inverse_op, element, attr_name, value, baseMVA, attributes)
                else:
                    values_list[time] = op( value , baseMVA )
        elif 'data_type' in attr and attr['data_type'] == 'cost_curve':
            if attr['cost_curve_type'] == 'polynomial':
                values_dict = attr['values']
                new_values = { int(power): coeff*(inverse_op(1.,baseMVA)**int(power)) \
                                for (power, coeff) in values_dict.items() }
                attr['values'] = new_values
            elif attr['cost_curve_type'] == 'piecewise':
                values_list_of_tuples = attr['values']
                new_values = [ ( normal_op(point,baseMVA), cost) \
                                for (point, cost) in values_list_of_tuples ]
                attr['values'] = new_values
        elif 'data_type' in attr and attr['data_type'] == 'fuel_curve':
            values_list_of_tuples = attr['values']
            new_values = [ ( normal_op(point,baseMVA), fuel) \
                            for (point, fuel) in values_list_of_tuples ]
            attr['values'] = new_values
        else: # recurse deeper
            for k,v in attr.items():
                _scale_by_baseMVA(normal_op, inverse_op, attr, k, v, baseMVA, attributes)
    elif attr_name in attributes:
        op = _get_op(normal_op, inverse_op, attr_name)
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
