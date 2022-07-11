#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module provides functions that create the modules for typical ACOPF formulations.

#TODO: document this with examples
"""
import pyomo.environ as pe
import operator as op
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen
from egret.model_library.defn import FlowType, CoordinateType
from egret.data.data_utils import map_items, zip_items
from math import pi, radians, degrees, cos, sin
from collections import OrderedDict
from egret.data.networkx_utils import get_networkx_graph
import networkx
from pyomo.common.collections.orderedset import OrderedSet
from pyomo.contrib.fbbt import interval

def _include_feasibility_slack(model, bus_names, bus_p_loads, bus_q_loads,
                               gens_by_bus, gen_attrs,
                               p_marginal_slack_penalty, q_marginal_slack_penalty):
    
    import egret.model_library.decl as decl

    p_over_gen_bounds = {k: (0, tx_utils.over_gen_limit(bus_p_loads[k], gens_by_bus[k], gen_attrs['p_max'])) for k in bus_names}
    decl.declare_var('p_over_generation', model=model, index_set=bus_names,
                     initialize=0., bounds=p_over_gen_bounds
                     )

    p_load_shed_bounds  = {k: (0, tx_utils.load_shed_limit(bus_p_loads[k], gens_by_bus[k], gen_attrs['p_min'])) for k in bus_names}
    decl.declare_var('p_load_shed', model=model, index_set=bus_names,
                     initialize=0., bounds=p_load_shed_bounds
                     )

    q_over_gen_bounds = {k: (0, tx_utils.over_gen_limit(bus_q_loads[k], gens_by_bus[k], gen_attrs['q_max'])) for k in bus_names}
    decl.declare_var('q_over_generation', model=model, index_set=bus_names,
                     initialize=0., bounds=q_over_gen_bounds
                     )

    q_load_shed_bounds  = {k: (0, tx_utils.load_shed_limit(bus_q_loads[k], gens_by_bus[k], gen_attrs['q_min'])) for k in bus_names}
    decl.declare_var('q_load_shed', model=model, index_set=bus_names,
                     initialize=0., bounds=q_load_shed_bounds
                     )
    p_rhs_kwargs = {'include_feasibility_load_shed':'p_load_shed', 'include_feasibility_over_generation':'p_over_generation'}
    q_rhs_kwargs = {'include_feasibility_load_shed':'q_load_shed', 'include_feasibility_over_generation':'q_over_generation'}

    penalty_expr = sum(p_marginal_slack_penalty * (model.p_over_generation[bus_name] + model.p_load_shed[bus_name])
                     + q_marginal_slack_penalty * (model.q_over_generation[bus_name] + model.q_load_shed[bus_name])
                    for bus_name in bus_names)
    return p_rhs_kwargs, q_rhs_kwargs, penalty_expr


def _validate_and_extract_slack_penalties(model_data):
    assert('load_mismatch_cost' in model_data.data['system'])
    assert('q_load_mismatch_cost' in model_data.data['system'])
    return model_data.data['system']['load_mismatch_cost'], model_data.data['system']['q_load_mismatch_cost']


def _create_base_power_ac_model(model_data, include_feasibility_slack=False, pw_cost_model='delta', keep_vars_for_out_of_service_elements=False):
    if keep_vars_for_out_of_service_elements:
        out_of_service_gens = tx_utils._get_out_of_service_gens(model_data)
        out_of_service_branches = tx_utils._get_out_of_service_branches(model_data)
    else:
        out_of_service_gens = list()
        out_of_service_branches = list()
                
    md = model_data.clone_in_service()
    tx_utils.scale_ModelData_to_pu(md, inplace=True)

    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    loads = dict(md.elements(element_type='load'))
    shunts = dict(md.elements(element_type='shunt'))

    gen_attrs = md.attributes(element_type='generator')
    bus_attrs = md.attributes(element_type='bus')
    branch_attrs = md.attributes(element_type='branch')

    inlet_branches_by_bus, outlet_branches_by_bus = \
        tx_utils.inlet_outlet_branches_by_bus(branches, buses)
    gens_by_bus = tx_utils.gens_by_bus(buses, gens)

    unique_bus_pairs = tx_utils.get_unique_bus_pairs(md)

    model = pe.ConcreteModel()

    ### declare (and fix) the loads at the buses
    bus_p_loads, bus_q_loads = tx_utils.dict_of_bus_loads(buses, loads)

    libbus.declare_var_pl(model, bus_attrs['names'], initialize=bus_p_loads)
    libbus.declare_var_ql(model, bus_attrs['names'], initialize=bus_q_loads)
    model.pl.fix()
    model.ql.fix()

    ### declare the fixed shunts at the buses
    bus_bs_fixed_shunts, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    libbus.declare_var_vmsq(model=model,
                            index_set=bus_attrs['names'],
                            initialize={k: v**2 for k, v in bus_attrs['vm'].items()},
                            bounds=zip_items({k: v**2 for k, v in bus_attrs['v_min'].items()},
                                             {k: v**2 for k, v in bus_attrs['v_max'].items()}))

    branch_w_index = {(v['from_bus'], v['to_bus']): v for v in branches.values()}

    def _c_bounds_rule(m, from_bus, to_bus):
        bdat = branch_w_index[(from_bus, to_bus)]

        theta_bounds = (bdat['angle_diff_min']*pi/180, bdat['angle_diff_max']*pi/180)
        vf_bounds = (bus_attrs['v_min'][from_bus], bus_attrs['v_max'][from_bus])
        vt_bounds = (bus_attrs['v_min'][to_bus], bus_attrs['v_max'][to_bus])
        c_bounds = interval.cos(*theta_bounds)
        c_bounds = interval.mul(*vf_bounds, *c_bounds)
        c_bounds = interval.mul(*vt_bounds, *c_bounds)
        return c_bounds
    libbranch.declare_var_c(model=model,
                            index_set=unique_bus_pairs,
                            initialize=1,
                            bounds=_c_bounds_rule)

    def _s_bounds_rule(m, from_bus, to_bus):
        bdat = branch_w_index[(from_bus, to_bus)]

        theta_bounds = (bdat['angle_diff_min']*pi/180, bdat['angle_diff_max']*pi/180)
        vf_bounds = (bus_attrs['v_min'][from_bus], bus_attrs['v_max'][from_bus])
        vt_bounds = (bus_attrs['v_min'][to_bus], bus_attrs['v_max'][to_bus])
        s_bounds = interval.sin(*theta_bounds)
        s_bounds = interval.mul(*vf_bounds, *s_bounds)
        s_bounds = interval.mul(*vt_bounds, *s_bounds)
        return s_bounds
    libbranch.declare_var_s(model=model, index_set=unique_bus_pairs, initialize=0,
                            bounds=_s_bounds_rule)

    ### include the feasibility slack for the bus balances
    p_rhs_kwargs = {}
    q_rhs_kwargs = {}
    if include_feasibility_slack:
        p_marginal_slack_penalty, q_marginal_slack_penalty = _validate_and_extract_slack_penalties(md)
        p_rhs_kwargs, q_rhs_kwargs, penalty_expr = _include_feasibility_slack(model, bus_attrs['names'],
                                                                              bus_p_loads, bus_q_loads,
                                                                              gens_by_bus, gen_attrs,
                                                                              p_marginal_slack_penalty,
                                                                              q_marginal_slack_penalty)

    ### declare the generator real and reactive power
    pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model, gen_attrs['names'], initialize=pg_init,
                          bounds=zip_items(gen_attrs['p_min'], gen_attrs['p_max'])
                          )

    qg_init = {k: (gen_attrs['q_min'][k] + gen_attrs['q_max'][k]) / 2.0 for k in gen_attrs['qg']}
    libgen.declare_var_qg(model, gen_attrs['names'], initialize=qg_init,
                          bounds=zip_items(gen_attrs['q_min'], gen_attrs['q_max'])
                          )

    ### declare the current flows in the branches
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(radians(bus_attrs['va'][k])) for k in bus_attrs['vm']}
    vj_init = {k: bus_attrs['vm'][k] * pe.sin(radians(bus_attrs['va'][k])) for k in bus_attrs['vm']}
    s_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
    s_lbub = dict()
    for k in branches.keys():
        if s_max[k] is None:
            s_lbub[k] = (None, None)
        else:
            s_lbub[k] = (-s_max[k],s_max[k])
    pf_bounds = s_lbub
    pt_bounds = s_lbub
    qf_bounds = s_lbub
    qt_bounds = s_lbub
    pf_init = dict()
    pt_init = dict()
    qf_init = dict()
    qt_init = dict()
    for branch_name, branch in branches.items():
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']
        y_matrix = tx_calc.calculate_y_matrix_from_branch(branch)
        ifr_init = tx_calc.calculate_ifr(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        ifj_init = tx_calc.calculate_ifj(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        itr_init = tx_calc.calculate_itr(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        itj_init = tx_calc.calculate_itj(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        pf_init[branch_name] = tx_calc.calculate_p(ifr_init, ifj_init, vr_init[from_bus], vj_init[from_bus])
        pt_init[branch_name] = tx_calc.calculate_p(itr_init, itj_init, vr_init[to_bus], vj_init[to_bus])
        qf_init[branch_name] = tx_calc.calculate_q(ifr_init, ifj_init, vr_init[from_bus], vj_init[from_bus])
        qt_init[branch_name] = tx_calc.calculate_q(itr_init, itj_init, vr_init[to_bus], vj_init[to_bus])

    libbranch.declare_var_pf(model=model,
                             index_set=branch_attrs['names'],
                             initialize=pf_init,
                             bounds=pf_bounds
                             )
    libbranch.declare_var_pt(model=model,
                             index_set=branch_attrs['names'],
                             initialize=pt_init,
                             bounds=pt_bounds
                             )
    libbranch.declare_var_qf(model=model,
                             index_set=branch_attrs['names'],
                             initialize=qf_init,
                             bounds=qf_bounds
                             )
    libbranch.declare_var_qt(model=model,
                             index_set=branch_attrs['names'],
                             initialize=qt_init,
                             bounds=qt_bounds
                             )

    ### declare the branch power flow constraints
    libbranch.declare_eq_branch_power(model=model,
                                      index_set=branch_attrs['names'],
                                      branches=branches
                                      )

    ### declare the pq balances
    libbus.declare_eq_p_balance(model=model,
                                index_set=bus_attrs['names'],
                                bus_p_loads=bus_p_loads,
                                gens_by_bus=gens_by_bus,
                                bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                inlet_branches_by_bus=inlet_branches_by_bus,
                                outlet_branches_by_bus=outlet_branches_by_bus,
                                **p_rhs_kwargs
                                )

    libbus.declare_eq_q_balance(model=model,
                                index_set=bus_attrs['names'],
                                bus_q_loads=bus_q_loads,
                                gens_by_bus=gens_by_bus,
                                bus_bs_fixed_shunts=bus_bs_fixed_shunts,
                                inlet_branches_by_bus=inlet_branches_by_bus,
                                outlet_branches_by_bus=outlet_branches_by_bus,
                                **q_rhs_kwargs
                                )

    ### declare the thermal limits
    libbranch.declare_ineq_s_branch_thermal_limit(model=model,
                                                  index_set=branch_attrs['names'],
                                                  branches=branches,
                                                  s_thermal_limits=s_max,
                                                  flow_type=FlowType.POWER
                                                  )

    # declare angle difference limits on interconnected buses
    libbranch.declare_ineq_angle_diff_branch_lbub_c_s(model=model,
                                                      index_set=branch_attrs['names'],
                                                      branches=branches
                                                      )

    # declare the generator cost objective
    p_costs = gen_attrs['p_cost']
    pw_pg_cost_gens = list(libgen.pw_gen_generator(gen_attrs['names'], costs=p_costs))
    if len(pw_pg_cost_gens) > 0:
        if pw_cost_model == 'delta':
            libgen.declare_var_delta_pg(model=model, index_set=pw_pg_cost_gens, p_costs=p_costs)
            libgen.declare_pg_delta_pg_con(model=model, index_set=pw_pg_cost_gens, p_costs=p_costs)
        else:
            libgen.declare_var_pg_cost(model=model, index_set=pw_pg_cost_gens, p_costs=p_costs)
            libgen.declare_piecewise_pg_cost_cons(model=model, index_set=pw_pg_cost_gens, p_costs=p_costs)
    libgen.declare_expression_pg_operating_cost(model=model, index_set=gen_attrs['names'], p_costs=p_costs, pw_formulation=pw_cost_model)
    obj_expr = sum(model.pg_operating_cost[gen_name] for gen_name in model.pg_operating_cost)
    q_costs = gen_attrs.get('q_cost', None)
    if q_costs is not None:
        pw_qg_cost_gens = list(libgen.pw_gen_generator(gen_attrs['names'], costs=q_costs))
        if len(pw_qg_cost_gens) > 0:
            if pw_cost_model == 'delta':
                libgen.declare_var_delta_qg(model=model, index_set=pw_qg_cost_gens, q_costs=q_costs)
                libgen.declare_qg_delta_qg_con(model=model, index_set=pw_qg_cost_gens, q_costs=q_costs)
            else:
                libgen.declare_var_qg_cost(model=model, index_set=pw_qg_cost_gens, q_costs=q_costs)
                libgen.declare_piecewise_qg_cost_cons(model=model, index_set=pw_qg_cost_gens, q_costs=q_costs)
        libgen.declare_expression_qg_operating_cost(model=model, index_set=gen_attrs['names'], q_costs=q_costs, pw_formulation=pw_cost_model)
        obj_expr += sum(model.qg_operating_cost[gen_name] for gen_name in model.qg_operating_cost)

    if include_feasibility_slack:
        obj_expr += penalty_expr

    model.obj = pe.Objective(expr=obj_expr)

    out_of_service_gens_set = set(out_of_service_gens)

    for gen_name, e in model.pg_operating_cost.items():
        if gen_name in out_of_service_gens_set:
            e.expr = 0

    if len(pw_pg_cost_gens) > 0:
        if pw_cost_model == 'delta':
            for gen_name, ndx in model.delta_pg_set:
                if gen_name in out_of_service_gens_set:
                    model.delta_pg[gen_name, ndx].set_value(0, skip_validation=True)
                    model.delta_pg[gen_name, ndx].fix()
                    model.pg_delta_pg_con[gen_name].deactivate()
        else:
            for gen_name, ndx in model.pg_piecewise_cost_set:
                if gen_name in out_of_service_gens_set:
                    model.pg_cost[gen_name].set_value(0, skip_validation=True)
                    model.pg_cost[gen_name].fix()
                    model.pg_piecewise_cost_cons[gen_name, ndx].deactivate()

    if q_costs is not None:
        for gen_name, e in model.qg_operating_cost.items():
            if gen_name in out_of_service_gens_set:
                e.expr = 0
        if len(pw_qg_cost_gens) > 0:
            if pw_cost_model == 'delta':
                for gen_name, ndx in model.delta_qg_set:
                    if gen_name in out_of_service_gens_set:
                        model.delta_qg[gen_name, ndx].set_value(0, skip_validation=True)
                        model.delta_qg[gen_name, ndx].fix()
                        model.qg_delta_qg_con[gen_name].deactivate()
            else:
                for gen_name, ndx in model.qg_piecewise_cost_set:
                    if gen_name in out_of_service_gens_set:
                        model.qg_cost[gen_name].set_value(0, skip_validation=True)
                        model.qg_cost[gen_name].fix()
                        model.qg_piecewise_cost_cons[gen_name, ndx].deactivate()

    for gen_name in out_of_service_gens:
        model.pg[gen_name].set_value(0, skip_validation=True)
        model.qg[gen_name].set_value(0, skip_validation=True)
        model.pg[gen_name].fix()
        model.qg[gen_name].fix()
        model_data.data['elements']['generator'][gen_name]['in_service'] = False
        md.data['elements']['generator'][gen_name]['in_service'] = False
    for branch_name in out_of_service_branches:
        model.pf[branch_name].set_value(0, skip_validation=True)
        model.pt[branch_name].set_value(0, skip_validation=True)
        model.qf[branch_name].set_value(0, skip_validation=True)
        model.qt[branch_name].set_value(0, skip_validation=True)
        model.pf[branch_name].fix()
        model.pt[branch_name].fix()
        model.qf[branch_name].fix()
        model.qt[branch_name].fix()
        model.eq_pf_branch[branch_name].deactivate()
        model.eq_pt_branch[branch_name].deactivate()
        model.eq_qf_branch[branch_name].deactivate()
        model.eq_qt_branch[branch_name].deactivate()
        model.ineq_sf_branch_thermal_limit[branch_name].deactivate()
        model.ineq_st_branch_thermal_limit[branch_name].deactivate()
        model.ineq_angle_diff_branch_lb[branch_name].deactivate()
        model.ineq_angle_diff_branch_ub[branch_name].deactivate()
        model_data.data['elements']['branch'][branch_name]['in_service'] = False
        md.data['elements']['branch'][branch_name]['in_service'] = False

    unique_bus_pairs_set = set(unique_bus_pairs)
    ubp = set(tx_utils.get_unique_bus_pairs(md))
    ubp_diff = unique_bus_pairs_set - ubp
    for from_bus, to_bus in ubp_diff:
        model.c[from_bus, to_bus].setlb(None)
        model.c[from_bus, to_bus].setub(None)
        model.s[from_bus, to_bus].setlb(None)
        model.s[from_bus, to_bus].setub(None)

    return model, md


def create_atan_acopf_model(model_data, include_feasibility_slack=False, pw_cost_model='delta', keep_vars_for_out_of_service_elements=False):
    model, md = _create_base_power_ac_model(model_data, include_feasibility_slack=include_feasibility_slack,
                                            pw_cost_model=pw_cost_model, keep_vars_for_out_of_service_elements=keep_vars_for_out_of_service_elements)

    unique_bus_pairs = tx_utils.get_unique_bus_pairs(md)
    for fb, tb in unique_bus_pairs:
        assert (tb, fb) not in unique_bus_pairs

    graph = get_networkx_graph(md)
    ref_bus = md.data['system']['reference_bus']
    cycle_basis = networkx.algorithms.cycle_basis(graph, root=ref_bus)

    cycle_basis_bus_pairs = OrderedSet()
    for cycle in cycle_basis:
        for ndx in range(len(cycle) - 1):
            b1 = cycle[ndx]
            b2 = cycle[ndx + 1]
            assert (b1, b2) in unique_bus_pairs or (b2, b1) in unique_bus_pairs
            if (b1, b2) in unique_bus_pairs:
                cycle_basis_bus_pairs.add((b1, b2))
            else:
                cycle_basis_bus_pairs.add((b2, b1))
        b1 = cycle[-1]
        b2 = cycle[0]
        assert (b1, b2) in unique_bus_pairs or (b2, b1) in unique_bus_pairs
        if (b1, b2) in unique_bus_pairs:
            cycle_basis_bus_pairs.add((b1, b2))
        else:
            cycle_basis_bus_pairs.add((b2, b1))

    branches = dict(md.elements(element_type='branch'))
    branch_w_index = {(v['from_bus'], v['to_bus']): v for v in branches.values()}
    def _dva_bounds_rule(m, from_bus, to_bus):
        bdat = branch_w_index[(from_bus, to_bus)]
        return max(-pi/2, bdat['angle_diff_min'] * pi / 180), min(pi/2, bdat['angle_diff_max'] * pi / 180)
    libbranch.declare_var_dva(model=model,
                              index_set=list(cycle_basis_bus_pairs),
                              initialize=0,
                              bounds=_dva_bounds_rule)
    libbranch.declare_eq_dva_arctan(model=model, index_set=list(cycle_basis_bus_pairs))
    libbranch.declare_eq_dva_cycle_sum(model=model, cycle_basis=cycle_basis, valid_bus_pairs=cycle_basis_bus_pairs)
    libbranch.declare_ineq_soc(model=model, index_set=list(unique_bus_pairs), use_outer_approximation=False)
    libbranch.declare_ineq_soc_ub(model=model, index_set=list(unique_bus_pairs))

    return model, md


def create_psv_acopf_model(model_data, include_feasibility_slack=False, pw_cost_model='delta', keep_vars_for_out_of_service_elements=False):
    model, md = _create_base_power_ac_model(model_data, include_feasibility_slack=include_feasibility_slack,
                                            pw_cost_model=pw_cost_model, keep_vars_for_out_of_service_elements=keep_vars_for_out_of_service_elements)
    bus_attrs = md.attributes(element_type='bus')
    unique_bus_pairs = tx_utils.get_unique_bus_pairs(md)

    # declare the polar voltages
    libbranch.declare_var_dva(model=model,
                              index_set=unique_bus_pairs,
                              initialize=0,
                              bounds=(-pi/2, pi/2))
    libbus.declare_var_vm(model,
                          bus_attrs['names'],
                          initialize=bus_attrs['vm'],
                          bounds=zip_items(bus_attrs['v_min'], bus_attrs['v_max']))

    va_bounds = {k: (-pi, pi) for k in bus_attrs['va']}
    libbus.declare_var_va(model,
                          bus_attrs['names'],
                          initialize=tx_utils.radians_from_degrees_dict(bus_attrs['va']),
                          bounds=va_bounds)

    # fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    ref_angle = md.data['system']['reference_bus_angle']
    model.va[ref_bus].fix(radians(ref_angle))

    # relate c, s, and vmsq to vm and va
    libbranch.declare_eq_delta_va(model=model,
                                  index_set=unique_bus_pairs)
    libbus.declare_eq_vmsq(model=model,
                           index_set=bus_attrs['names'],
                           coordinate_type=CoordinateType.POLAR)
    libbranch.declare_eq_c(model=model,
                           index_set=unique_bus_pairs,
                           coordinate_type=CoordinateType.POLAR)
    libbranch.declare_eq_s(model=model,
                           index_set=unique_bus_pairs,
                           coordinate_type=CoordinateType.POLAR)

    return model, md


def create_rsv_acopf_model(model_data, include_feasibility_slack=False, pw_cost_model='delta', keep_vars_for_out_of_service_elements=False):
    model, md = _create_base_power_ac_model(model_data, include_feasibility_slack=include_feasibility_slack,
                                            pw_cost_model=pw_cost_model, keep_vars_for_out_of_service_elements=keep_vars_for_out_of_service_elements)
    bus_attrs = md.attributes(element_type='bus')
    unique_bus_pairs = tx_utils.get_unique_bus_pairs(md)

    # declare the rectangular voltages
    neg_v_max = map_items(op.neg, bus_attrs['v_max'])
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(radians(bus_attrs['va'][k])) for k in bus_attrs['vm']}
    libbus.declare_var_vr(model, bus_attrs['names'], initialize=vr_init,
                          bounds=zip_items(neg_v_max, bus_attrs['v_max'])
                          )

    vj_init = {k: bus_attrs['vm'][k] * pe.sin(radians(bus_attrs['va'][k])) for k in bus_attrs['vm']}
    libbus.declare_var_vj(model, bus_attrs['names'], initialize=vj_init,
                          bounds=zip_items(neg_v_max, bus_attrs['v_max'])
                          )

    # fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    ref_angle = md.data['system']['reference_bus_angle']
    if ref_angle != 0.0:
        libbus.declare_eq_ref_bus_nonzero(model, ref_angle, ref_bus)
    else:
        model.vj[ref_bus].fix(0.0)
        model.vr[ref_bus].setlb(bus_attrs['v_min'][ref_bus])

    # relate c, s, and vmsq to vm and va
    libbus.declare_eq_vmsq(model=model,
                           index_set=bus_attrs['names'],
                           coordinate_type=CoordinateType.RECTANGULAR)
    libbranch.declare_eq_c(model=model,
                           index_set=unique_bus_pairs,
                           coordinate_type=CoordinateType.RECTANGULAR)
    libbranch.declare_eq_s(model=model,
                           index_set=unique_bus_pairs,
                           coordinate_type=CoordinateType.RECTANGULAR)

    return model, md


def create_riv_acopf_model(model_data, include_feasibility_slack=False, pw_cost_model='delta'):
    md = model_data.clone_in_service()
    tx_utils.scale_ModelData_to_pu(md, inplace = True)

    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    loads = dict(md.elements(element_type='load'))
    shunts = dict(md.elements(element_type='shunt'))

    gen_attrs = md.attributes(element_type='generator')
    bus_attrs = md.attributes(element_type='bus')
    branch_attrs = md.attributes(element_type='branch')
    load_attrs = md.attributes(element_type='load')
    shunt_attrs = md.attributes(element_type='shunt')

    inlet_branches_by_bus, outlet_branches_by_bus = \
        tx_utils.inlet_outlet_branches_by_bus(branches, buses)
    gens_by_bus = tx_utils.gens_by_bus(buses, gens)

    model = pe.ConcreteModel()

    ### declare (and fix) the loads at the buses
    bus_p_loads, bus_q_loads = tx_utils.dict_of_bus_loads(buses, loads)

    libbus.declare_var_pl(model, bus_attrs['names'], initialize=bus_p_loads)
    libbus.declare_var_ql(model, bus_attrs['names'], initialize=bus_q_loads)
    model.pl.fix()
    model.ql.fix()

    ### declare the fixed shunts at the buses
    bus_bs_fixed_shunts, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    ### declare the rectangular voltages
    neg_v_max = map_items(op.neg, bus_attrs['v_max'])
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(radians(bus_attrs['va'][k])) for k in bus_attrs['vm']}
    libbus.declare_var_vr(model, bus_attrs['names'], initialize=vr_init,
                          bounds=zip_items(neg_v_max, bus_attrs['v_max'])
                          )

    vj_init = {k: bus_attrs['vm'][k] * pe.sin(radians(bus_attrs['va'][k])) for k in bus_attrs['vm']}
    libbus.declare_var_vj(model, bus_attrs['names'], initialize=vj_init,
                          bounds=zip_items(neg_v_max, bus_attrs['v_max'])
                          )

    ### include the feasibility slack for the bus balances
    p_rhs_kwargs = {}
    q_rhs_kwargs = {}
    if include_feasibility_slack:
        p_marginal_slack_penalty, q_marginal_slack_penalty = _validate_and_extract_slack_penalties(md)
        p_rhs_kwargs, q_rhs_kwargs, penalty_expr = _include_feasibility_slack(model, bus_attrs['names'],
                                                                              bus_p_loads, bus_q_loads,
                                                                              gens_by_bus, gen_attrs,
                                                                              p_marginal_slack_penalty,
                                                                              q_marginal_slack_penalty)

    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    ref_angle = md.data['system']['reference_bus_angle']
    if ref_angle != 0.0:
        libbus.declare_eq_ref_bus_nonzero(model, ref_angle, ref_bus)
    else:
        model.vj[ref_bus].fix(0.0)
        model.vr[ref_bus].setlb(0.0)

    ### declare the generator real and reactive power
    pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model, gen_attrs['names'], initialize=pg_init,
                          bounds=zip_items(gen_attrs['p_min'], gen_attrs['p_max'])
                          )

    qg_init = {k: (gen_attrs['q_min'][k] + gen_attrs['q_max'][k]) / 2.0 for k in gen_attrs['qg']}
    libgen.declare_var_qg(model, gen_attrs['names'], initialize=qg_init,
                          bounds=zip_items(gen_attrs['q_min'], gen_attrs['q_max'])
                          )

    ### declare the current flows in the branches
    branch_currents = tx_utils.dict_of_branch_currents(branches, buses)
    s_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
    if_bounds = dict()
    it_bounds = dict()
    ifr_init = dict()
    ifj_init = dict()
    itr_init = dict()
    itj_init = dict()
    for branch_name, branch in branches.items():
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']
        y_matrix = tx_calc.calculate_y_matrix_from_branch(branch)
        ifr_init[branch_name] = tx_calc.calculate_ifr(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                                      vj_init[to_bus], y_matrix)
        ifj_init[branch_name] = tx_calc.calculate_ifj(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                                      vj_init[to_bus], y_matrix)
        itr_init[branch_name] = tx_calc.calculate_itr(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                                      vj_init[to_bus], y_matrix)
        itj_init[branch_name] = tx_calc.calculate_itj(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                                      vj_init[to_bus], y_matrix)
        if s_max[branch_name] is None:
            if_bounds[branch_name] = (None, None)
            it_bounds[branch_name] = (None, None)
        else:
            if_max = s_max[branch_name] / buses[branches[branch_name]['from_bus']]['v_min']
            it_max = s_max[branch_name] / buses[branches[branch_name]['to_bus']]['v_min']
            if_bounds[branch_name] = (-if_max, if_max)
            it_bounds[branch_name] = (-it_max, it_max)

    libbranch.declare_var_ifr(model=model,
                              index_set=branch_attrs['names'],
                              initialize=ifr_init,
                              bounds=if_bounds
                              )
    libbranch.declare_var_ifj(model=model,
                              index_set=branch_attrs['names'],
                              initialize=ifj_init,
                              bounds=if_bounds
                              )
    libbranch.declare_var_itr(model=model,
                              index_set=branch_attrs['names'],
                              initialize=itr_init,
                              bounds=it_bounds
                              )
    libbranch.declare_var_itj(model=model,
                              index_set=branch_attrs['names'],
                              initialize=itj_init,
                              bounds=it_bounds
                              )

    ir_init = dict()
    ij_init = dict()
    for bus_name, bus in buses.items():
        ir_expr = sum([ifr_init[branch_name] for branch_name in outlet_branches_by_bus[bus_name]])
        ir_expr += sum([itr_init[branch_name] for branch_name in inlet_branches_by_bus[bus_name]])
        ij_expr = sum([ifj_init[branch_name] for branch_name in outlet_branches_by_bus[bus_name]])
        ij_expr += sum([itj_init[branch_name] for branch_name in inlet_branches_by_bus[bus_name]])

        if bus_gs_fixed_shunts[bus_name] != 0.0:
            ir_expr += bus_gs_fixed_shunts[bus_name] * vr_init[bus_name]
            ij_expr += bus_gs_fixed_shunts[bus_name] * vj_init[bus_name]
        if bus_bs_fixed_shunts[bus_name] != 0.0:
            ir_expr += bus_bs_fixed_shunts[bus_name] * vj_init[bus_name]
            ij_expr += bus_bs_fixed_shunts[bus_name] * vr_init[bus_name]

        ir_init[bus_name] = ir_expr
        ij_init[bus_name] = ij_expr

    # TODO: Implement better bounds (?) for these aggregated variables -- note, these are unbounded in old Egret
    libbus.declare_var_ir_aggregation_at_bus(model=model,
                                             index_set=bus_attrs['names'],
                                             initialize=ir_init,
                                             bounds=(None,None)
                                             )
    libbus.declare_var_ij_aggregation_at_bus(model=model,
                                             index_set=bus_attrs['names'],
                                             initialize=ij_init,
                                             bounds=(None,None)
                                             )

    ### declare the branch current flow constraints
    libbranch.declare_eq_branch_current(model=model,
                                        index_set=branch_attrs['names'],
                                        branches=branches
                                        )

    ### declare the ir/ij_aggregation constraints
    libbus.declare_eq_i_aggregation_at_bus(model=model,
                                           index_set=bus_attrs['names'],
                                           bus_bs_fixed_shunts=bus_bs_fixed_shunts,
                                           bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                           inlet_branches_by_bus=inlet_branches_by_bus,
                                           outlet_branches_by_bus=outlet_branches_by_bus
                                           )

    ### declare the pq balances
    libbus.declare_eq_p_balance_with_i_aggregation(model=model,
                                                   index_set=bus_attrs['names'],
                                                   bus_p_loads=bus_p_loads,
                                                   gens_by_bus=gens_by_bus,
                                                   **p_rhs_kwargs
                                                   )

    libbus.declare_eq_q_balance_with_i_aggregation(model=model,
                                                   index_set=bus_attrs['names'],
                                                   bus_q_loads=bus_q_loads,
                                                   gens_by_bus=gens_by_bus,
                                                   **q_rhs_kwargs
                                                   )

    ### declare the thermal limits
    libbranch.declare_ineq_s_branch_thermal_limit(model=model,
                                                  index_set=branch_attrs['names'],
                                                  branches=branches,
                                                  s_thermal_limits=s_max,
                                                  flow_type=FlowType.CURRENT
                                                  )

    ### declare the voltage min and max inequalities
    libbus.declare_ineq_vm_bus_lbub(model=model,
                                    index_set=bus_attrs['names'],
                                    buses=buses,
                                    coordinate_type=CoordinateType.RECTANGULAR
                                    )

    ### declare angle difference limits on interconnected buses
    libbranch.declare_ineq_angle_diff_branch_lbub(model=model,
                                                  index_set=branch_attrs['names'],
                                                  branches=branches,
                                                  coordinate_type=CoordinateType.RECTANGULAR)

    ### declare the generator cost objective
    p_costs = gen_attrs['p_cost']
    pw_pg_cost_gens = list(libgen.pw_gen_generator(gen_attrs['names'], costs=p_costs))
    if len(pw_pg_cost_gens) > 0:
        if pw_cost_model == 'delta':
            libgen.declare_var_delta_pg(model=model, index_set=pw_pg_cost_gens, p_costs=p_costs)
            libgen.declare_pg_delta_pg_con(model=model, index_set=pw_pg_cost_gens, p_costs=p_costs)
        else:
            libgen.declare_var_pg_cost(model=model, index_set=pw_pg_cost_gens, p_costs=p_costs)
            libgen.declare_piecewise_pg_cost_cons(model=model, index_set=pw_pg_cost_gens, p_costs=p_costs)
    libgen.declare_expression_pg_operating_cost(model=model, index_set=gen_attrs['names'], p_costs=p_costs, pw_formulation=pw_cost_model)
    obj_expr = sum(model.pg_operating_cost[gen_name] for gen_name in model.pg_operating_cost)
    q_costs = gen_attrs.get('q_cost', None)
    if q_costs is not None:
        pw_qg_cost_gens = list(libgen.pw_gen_generator(gen_attrs['names'], costs=q_costs))
        if len(pw_qg_cost_gens) > 0:
            if pw_cost_model == 'delta':
                libgen.declare_var_delta_qg(model=model, index_set=pw_qg_cost_gens, q_costs=q_costs)
                libgen.declare_qg_delta_qg_con(model=model, index_set=pw_qg_cost_gens, q_costs=q_costs)
            else:
                libgen.declare_var_qg_cost(model=model, index_set=pw_qg_cost_gens, q_costs=q_costs)
                libgen.declare_piecewise_qg_cost_cons(model=model, index_set=pw_qg_cost_gens, q_costs=q_costs)
        libgen.declare_expression_qg_operating_cost(model=model, index_set=gen_attrs['names'], q_costs=q_costs, pw_formulation=pw_cost_model)
        obj_expr += sum(model.qg_operating_cost[gen_name] for gen_name in model.qg_operating_cost)

    if include_feasibility_slack:
        obj_expr += penalty_expr

    model.obj = pe.Objective(expr=obj_expr)

    return model, md


def solve_acopf(model_data,
                solver,
                timelimit = None,
                solver_tee = True,
                symbolic_solver_labels = False,
                options = None,
                acopf_model_generator = create_psv_acopf_model,
                return_model = False,
                return_results = False,
                **kwargs):
    '''
    Create and solve a new acopf model

    Parameters
    ----------
    model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded.
    solver : str or pyomo.opt.base.solvers.OptSolver
        Either a string specifying a pyomo solver name, or an instantiated pyomo solver
    timelimit : float (optional)
        Time limit for dcopf run. Default of None results in no time
        limit being set.
    solver_tee : bool (optional)
        Display solver log. Default is True.
    symbolic_solver_labels : bool (optional)
        Use symbolic solver labels. Useful for debugging; default is False.
    options : dict (optional)
        Other options to pass into the solver. Default is dict().
    acopf_model_generator : function (optional)
        Function for generating the acopf model. Default is
        egret.models.acopf.create_psv_acopf_model
    return_model : bool (optional)
        If True, returns the pyomo model object
    return_results : bool (optional)
        If True, returns the pyomo results object
    kwargs : dictionary (optional)
        Additional arguments for building model
    '''

    import pyomo.environ as pe
    from pyomo.environ import value
    from egret.common.solver_interface import _solve_model
    from egret.model_library.transmission.tx_utils import \
        scale_ModelData_to_pu, unscale_ModelData_to_pu

    m, md = acopf_model_generator(model_data, **kwargs)

    m.dual = pe.Suffix(direction=pe.Suffix.IMPORT)

    m, results = _solve_model(m,solver,timelimit=timelimit,solver_tee=solver_tee,
                              symbolic_solver_labels=symbolic_solver_labels,solver_options=options)

    # save results data to ModelData object
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))

    md.data['system']['total_cost'] = value(m.obj)

    for g,g_dict in gens.items():
        g_dict['pg'] = value(m.pg[g])
        g_dict['qg'] = value(m.qg[g])

    for b,b_dict in buses.items():
        b_dict['lmp'] = value(m.dual[m.eq_p_balance[b]])
        b_dict['qlmp'] = value(m.dual[m.eq_q_balance[b]])
        b_dict['pl'] = value(m.pl[b])
        if hasattr(m, 'vj'):
            b_dict['vm'] = tx_calc.calculate_vm_from_vj_vr(value(m.vj[b]), value(m.vr[b]))
            b_dict['va'] = tx_calc.calculate_va_degrees_from_vj_vr(value(m.vj[b]), value(m.vr[b]))
        else:
            b_dict['vm'] = value(m.vm[b])
            b_dict['va'] = degrees(value(m.va[b]))
        if hasattr(m, 'p_load_shed'):
            b_dict['p_balance_violation'] = value(m.p_load_shed[b]) - value(m.p_over_generation[b])
        if hasattr(m, 'q_load_shed'):
            b_dict['q_balance_violation'] = value(m.q_load_shed[b]) - value(m.q_over_generation[b])

    for k, k_dict in branches.items():
        if hasattr(m,'pf'):
            k_dict['pf'] = value(m.pf[k])
            k_dict['pt'] = value(m.pt[k])
            k_dict['qf'] = value(m.qf[k])
            k_dict['qt'] = value(m.qt[k])
        if hasattr(m,'irf'):
            b = k_dict['from_bus']
            k_dict['pf'] = value(tx_calc.calculate_p(value(m.ifr[k]), value(m.ifj[k]), value(m.vr[b]), value(m.vj[b])))
            k_dict['qf'] = value(tx_calc.calculate_q(value(m.ifr[k]), value(m.ifj[k]), value(m.vr[b]), value(m.vj[b])))
            b = k_dict['to_bus']
            k_dict['pt'] = value(tx_calc.calculate_p(value(m.itr[k]), value(m.itj[k]), value(m.vr[b]), value(m.vj[b])))
            k_dict['qt'] = value(tx_calc.calculate_q(value(m.itr[k]), value(m.itj[k]), value(m.vr[b]), value(m.vj[b])))

    unscale_ModelData_to_pu(md, inplace=True)

    if return_model and return_results:
        return md, m, results
    elif return_model:
        return md, m
    elif return_results:
        return md, results
    return md


if __name__ == '__main__':
    import os
    from egret.parsers.matpower_parser import create_ModelData

    path = os.path.dirname(__file__)
    filename = 'pglib_opf_case14_ieee.m'
    matpower_file = os.path.join(path, '../../download/pglib-opf/', filename)
    model_data = create_ModelData(matpower_file)
    kwargs = {'include_feasibility_slack':False}
    md,m,results = solve_acopf(model_data, "ipopt",acopf_model_generator=create_psv_acopf_model,return_model=True, return_results=True,**kwargs)
    md,m,results = solve_acopf(model_data, "ipopt",acopf_model_generator=create_rsv_acopf_model,return_model=True, return_results=True,**kwargs)
    md,m,results = solve_acopf(model_data, "ipopt",acopf_model_generator=create_riv_acopf_model,return_model=True, return_results=True,**kwargs)


