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


def create_psv_acpf_model(model_data):
    md = tx_utils.scale_ModelData_to_pu(model_data)

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
    s_lbub = {k: (-s_max[k],s_max[k]) for k in branches.keys()}
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

    ### declare the voltage min and max inequalities
    libbus.declare_ineq_vm_bus_lbub(model=model,
                                    index_set=bus_attrs['names'],
                                    buses=buses,
                                    coordinate_type=CoordinateType.POLAR
                                    )

    model.obj = pe.Objective(expr=1.0)

    # fix pg at all PV buses (except the refbus)
    found_refbus = False

    for bus_name in bus_attrs['names']:
        for gen_name in gens_by_bus[bus_name]:
            model.pg[gen_name].fix(0.0)
            model.qg[gen_name].fix(0.0)

    for g, gen in grid.gens():
        if not gen.in_service:
            m.pg[g].fix(0.0)
            m.qg[g].fix(0.0)
            continue
        b = gen.attached_bus_name
        bus = grid.get_bus(b)
        if b == grid.reference_bus_name:
            if found_refbus:
                m.pg[g].fix(gen.pg)
            else:
                found_refbus = True
        else:
            m.pg[g].fix(gen.pg)

    # fix vm at all PV buses
    for g, gen in grid.gens():
        if gen.in_service:
            b = gen.attached_bus_name
            m.vm[b].fix(gen.Vg)

    return model
    model.pl.fix()
    model.ql.fix()

if __name__ == '__main__':
    import os

    path = os.path.dirname(__file__)

    matpower_file = os.path.join(path, '../../download/pglib-opf/pglib_opf_case14_ieee.m')

    case = 'pglib_opf_case14_ieee'

    assert os.path.isfile(matpower_file)

    from egret.parsers.matpower_parser import create_ModelData

    md = create_ModelData(matpower_file)

    model = create_riv_acopf_model(md)

    solver = pe.SolverFactory('ipopt')

    results = solver.solve(model, tee=True)


    model = create_rsv_acopf_model(md)

    solver = pe.SolverFactory('ipopt')

    results = solver.solve(model, tee=True)



    model = create_psv_acopf_model(md)

    solver = pe.SolverFactory('ipopt')

    results = solver.solve(model, tee=True)

