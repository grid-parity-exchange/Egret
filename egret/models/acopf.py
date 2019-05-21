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
from egret.data.model_data import map_items, zip_items
from math import pi


def create_psv_acopf_model(model_data):
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

    ### declare the polar voltages
    libbus.declare_var_vm(model, bus_attrs['names'], initialize=bus_attrs['vm'],
                          bounds=zip_items(bus_attrs['v_min'], bus_attrs['v_max'])
                          )

    va_bounds = {k: (-pi, pi) for k in bus_attrs['va']}
    libbus.declare_var_va(model, bus_attrs['names'], initialize=bus_attrs['va'],
                          bounds=va_bounds
                          )

    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    model.va[ref_bus].fix(0.0)

    ref_angle = md.data['system']['reference_bus_angle']
    if ref_angle != 0.0:
        raise ValueError('The RIV ACOPF formulation currently only supports'
                         ' a reference bus angle of 0 degrees, but an angle'
                         ' of {} degrees was found.'.format(ref_angle))

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
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    vj_init = {k: bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]) for k in bus_attrs['vm']}
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
                                      branches=branches,
                                      branch_attrs=branch_attrs,
                                      coordinate_type=CoordinateType.POLAR
                                      )

    ### declare the pq balances
    libbus.declare_eq_p_balance(model=model,
                                index_set=bus_attrs['names'],
                                bus_p_loads=bus_p_loads,
                                gens_by_bus=gens_by_bus,
                                bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                inlet_branches_by_bus=inlet_branches_by_bus,
                                outlet_branches_by_bus=outlet_branches_by_bus,
                                coordinate_type=CoordinateType.POLAR
                                )

    libbus.declare_eq_q_balance(model=model,
                                index_set=bus_attrs['names'],
                                bus_q_loads=bus_q_loads,
                                gens_by_bus=gens_by_bus,
                                bus_bs_fixed_shunts=bus_bs_fixed_shunts,
                                inlet_branches_by_bus=inlet_branches_by_bus,
                                outlet_branches_by_bus=outlet_branches_by_bus,
                                coordinate_type=CoordinateType.POLAR
                                )

    ### declare the thermal limits
    libbranch.declare_ineq_s_branch_thermal_limit(model=model,
                                                  index_set=branch_attrs['names'],
                                                  branches=branches,
                                                  s_thermal_limits=s_max,
                                                  flow_type=FlowType.POWER
                                                  )

    ### declare the voltage min and max inequalities
    libbus.declare_ineq_vm_bus_lbub(model=model,
                                    index_set=bus_attrs['names'],
                                    buses=buses,
                                    coordinate_type=CoordinateType.POLAR
                                    )

    ### declare angle difference limits on interconnected buses
    libbranch.declare_ineq_angle_diff_branch_lbub(model=model,
                                                  index_set=branch_attrs['names'],
                                                  branches=branches,
                                                  coordinate_type=CoordinateType.POLAR
                                                  )

    ### declare the generator cost objective
    libgen.declare_expression_pgqg_operating_cost(model=model,
                                                  index_set=gen_attrs['names'],
                                                  p_costs=gen_attrs['p_cost'],
                                                  q_costs=gen_attrs.get('q_cost', None)
                                                  )

    obj_expr = sum(model.pg_operating_cost[gen_name] for gen_name in model.pg_operating_cost)
    if hasattr(model, 'qg_operating_cost'):
        obj_expr += sum(model.qg_operating_cost[gen_name] for gen_name in model.qg_operating_cost)

    model.obj = pe.Objective(expr=obj_expr)

    return model


def create_rsv_acopf_model(model_data):
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
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    libbus.declare_var_vr(model, bus_attrs['names'], initialize=vr_init,
                          bounds=zip_items(neg_v_max, bus_attrs['v_max'])
                          )

    vj_init = {k: bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    libbus.declare_var_vj(model, bus_attrs['names'], initialize=vj_init,
                          bounds=zip_items(neg_v_max, bus_attrs['v_max'])
                          )

    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    model.vj[ref_bus].fix(0.0)
    model.vr[ref_bus].setlb(0.0)

    ref_angle = md.data['system']['reference_bus_angle']
    if ref_angle != 0.0:
        raise ValueError('The RSV ACOPF formulation currently only supports'
                         ' a reference bus angle of 0 degrees, but an angle'
                         ' of {} degrees was found.'.format(ref_angle))

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
                                      branches=branches,
                                      branch_attrs=branch_attrs,
                                      coordinate_type=CoordinateType.RECTANGULAR
                                      )

    ### declare the pq balances
    libbus.declare_eq_p_balance(model=model,
                                index_set=bus_attrs['names'],
                                bus_p_loads=bus_p_loads,
                                gens_by_bus=gens_by_bus,
                                bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                inlet_branches_by_bus=inlet_branches_by_bus,
                                outlet_branches_by_bus=outlet_branches_by_bus,
                                coordinate_type=CoordinateType.RECTANGULAR
                                )

    libbus.declare_eq_q_balance(model=model,
                                index_set=bus_attrs['names'],
                                bus_q_loads=bus_q_loads,
                                gens_by_bus=gens_by_bus,
                                bus_bs_fixed_shunts=bus_bs_fixed_shunts,
                                inlet_branches_by_bus=inlet_branches_by_bus,
                                outlet_branches_by_bus=outlet_branches_by_bus,
                                coordinate_type=CoordinateType.RECTANGULAR
                                )

    ### declare the thermal limits
    libbranch.declare_ineq_s_branch_thermal_limit(model=model,
                                                  index_set=branch_attrs['names'],
                                                  branches=branches,
                                                  s_thermal_limits=s_max,
                                                  flow_type=FlowType.POWER
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
                                                  coordinate_type=CoordinateType.RECTANGULAR
                                                  )

    ### declare the generator cost objective
    libgen.declare_expression_pgqg_operating_cost(model=model,
                                                  index_set=gen_attrs['names'],
                                                  p_costs=gen_attrs['p_cost'],
                                                  q_costs=gen_attrs.get('q_cost', None)
                                                  )

    obj_expr = sum(model.pg_operating_cost[gen_name] for gen_name in model.pg_operating_cost)
    if hasattr(model, 'qg_operating_cost'):
        obj_expr += sum(model.qg_operating_cost[gen_name] for gen_name in model.qg_operating_cost)

    model.obj = pe.Objective(expr=obj_expr)

    return model


def create_riv_acopf_model(model_data):
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
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    libbus.declare_var_vr(model, bus_attrs['names'], initialize=vr_init,
                          bounds=zip_items(neg_v_max, bus_attrs['v_max'])
                          )

    vj_init = {k: bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    libbus.declare_var_vj(model, bus_attrs['names'], initialize=vj_init,
                          bounds=zip_items(neg_v_max, bus_attrs['v_max'])
                          )

    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    model.vj[ref_bus].fix(0.0)
    model.vr[ref_bus].setlb(0.0)

    ref_angle = md.data['system']['reference_bus_angle']
    if ref_angle != 0.0:
        raise ValueError('The RIV ACOPF formulation currently only supports'
                         ' a reference bus angle of 0 degrees, but an angle'
                         ' of {} degrees was found.'.format(ref_angle))

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
                                                   )

    libbus.declare_eq_q_balance_with_i_aggregation(model=model,
                                                   index_set=bus_attrs['names'],
                                                   bus_q_loads=bus_q_loads,
                                                   gens_by_bus=gens_by_bus,
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
                                                  coordinate_type=CoordinateType.RECTANGULAR
                                                  )

    ### declare the generator cost objective
    libgen.declare_expression_pgqg_operating_cost(model=model,
                                                  index_set=gen_attrs['names'],
                                                  p_costs=gen_attrs['p_cost'],
                                                  q_costs=gen_attrs.get('q_cost', None)
                                                  )

    obj_expr = sum(model.pg_operating_cost[gen_name] for gen_name in model.pg_operating_cost)
    if hasattr(model, 'qg_operating_cost'):
        obj_expr += sum(model.qg_operating_cost[gen_name] for gen_name in model.qg_operating_cost)

    model.obj = pe.Objective(expr=obj_expr)

    return model


def solve_acopf(model_data,
                solver,
                timelimit = None,
                solver_tee = True,
                symbolic_solver_labels = False,
                options = None,
                acopf_model_generator = create_psv_acopf_model,
                return_model = False):
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

    '''

    import pyomo.environ as pe
    from pyomo.environ import value
    from egret.common.solver_interface import _solve_model
    from egret.model_library.transmission.tx_utils import \
        scale_ModelData_to_pu, unscale_ModelData_to_pu

    m = acopf_model_generator(model_data)

    m.dual = pe.Suffix(direction=pe.Suffix.IMPORT)

    m, results = _solve_model(m,solver,timelimit=timelimit,solver_tee=solver_tee,
                              symbolic_solver_labels=symbolic_solver_labels,options=options)

    md = model_data

    # save results data to ModelData object
    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    storage = dict(md.elements(element_type='storage'))
    zones = dict(md.elements(element_type='zone'))
    areas = dict(md.elements(element_type='area'))

    md.data['system']['total_cost'] = value(m.obj)

    for g,g_dict in gens.items():
        g_dict['pg'] = value(m.pg[g])
        g_dict['qg'] = value(m.qg[g])


    for b,b_dict in buses.items():
        b_dict['lmp'] = value(m.dual[m.eq_p_balance[b]])
        b_dict['qlmp'] = value(m.dual[m.eq_q_balance[b]])
        b_dict['vm'] = value(m.vm[b])
        b_dict['va'] = value(m.va[b])
        b_dict['pl'] = value(m.pl[b])

    for k, k_dict in branches.items():
        k_dict['pf'] = value(m.pf[k])
        k_dict['pt'] = value(m.pt[k])
        k_dict['qf'] = value(m.qf[k])
        k_dict['qt'] = value(m.qt[k])

    unscale_ModelData_to_pu(md, inplace=True)


    if return_model:
        return md, m
    return md

if __name__ == '__main__':
    import os
    from egret.parsers.matpower_parser import create_ModelData

    path = os.path.dirname(__file__)
    filename = 'pglib_opf_case3_lmbd.m'
    matpower_file = os.path.join(path, '../../download/pglib-opf/', filename)
    md = create_ModelData(matpower_file)
    solve_acopf(md, "ipopt")