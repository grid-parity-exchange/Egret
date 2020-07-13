
#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

import pyomo.environ as pe
import pao.bilevel as bi
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen
import egret.model_library.extensions.subproblem_bilevel_nk as subcons
from egret.model_library.defn import CoordinateType, ApproximationType
from math import pi, radians


def create_gdp_subproblem(model, model_data, include_angle_diff_limits=False):
    md = model_data
    tx_utils.scale_ModelData_to_pu(md, inplace = True)

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

    model.subproblem = bi.SubModel(fixed=(model.u, model.v, model.w))

    ### declare (and fix) the loads at the buses
    bus_p_loads, _ = tx_utils.dict_of_bus_loads(buses, loads)
    buses_with_loads = list(k for k in bus_p_loads.keys() if bus_p_loads[k] != 0.)

    libbus.declare_var_pl(model.subproblem, bus_attrs['names'], initialize=bus_p_loads)
    model.subproblem.pl.fix()

    ### declare the fixed shunts at the buses
    _, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)

    ### declare the polar voltages
    va_bounds = {k: (-pi, pi) for k in bus_attrs['va']}
    libbus.declare_var_va(model.subproblem, bus_attrs['names'], initialize=bus_attrs['va'],
                          bounds=va_bounds
                          )

    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    ref_angle = md.data['system']['reference_bus_angle']
    model.subproblem.va[ref_bus].fix(radians(ref_angle))

    ### declare the generator real power
    pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model.subproblem, gen_attrs['names'], initialize=pg_init,
                          bounds=zip_items(gen_attrs['p_min'], gen_attrs['p_max'])
                          )

    ### declare the current flows in the branches
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    vj_init = {k: bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    p_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
    p_lbub = {k: (-p_max[k],p_max[k]) for k in branches.keys()}
    pf_bounds = p_lbub
    pf_init = dict()
    for branch_name, branch in branches.items():
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']
        y_matrix = tx_calc.calculate_y_matrix_from_branch(branch)
        ifr_init = tx_calc.calculate_ifr(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        ifj_init = tx_calc.calculate_ifj(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        pf_init[branch_name] = tx_calc.calculate_p(ifr_init, ifj_init, vr_init[from_bus], vj_init[from_bus])

    libbranch.declare_var_pf(model=model.subproblem,
                             index_set=branch_attrs['names'],
                             initialize=pf_init,
                             bounds=pf_bounds
                             )

    # need to include variable references on subproblem to variables, which exist on the master block
    bi.components.varref(model.subproblem)

    ### declare the branch power flow disjuncts (LHS is status quo, RHS is compromised)
    libbranch.declare_eq_branch_power_btheta_approx(model=model.subproblem,
                                                    index_set=branch_attrs['names'],
                                                    branches=branches
                                                    )
    subcons.declare_eq_branch_power_off(model=model.subproblem,
                                        index_set=branch_attrs['names'],
                                        branches=branches
                                        )
    subcons.disjunctify(model=model.subproblem,
                        indicator_name='pf_branch_indicator',
                        disjunct_name='pf_branch_disjunct',
                        LHS_disjunct_set=model.subproblem.eq_pf_branch,
                        RHS_disjunct_set=model.subproblem.eq_pf_branch_off
                        )

    ### declare the load shed disjuncts (LHS is status quo, RHS is compromised)
    subcons.declare_ineq_load_shed_ub(model=model.subproblem,
                                      index_set=buses_with_loads)
    subcons.declare_ineq_load_shed_lb(model=model.subproblem,
                                      index_set=buses_with_loads)
    subcons.declare_ineq_load_shed_lb_off(model=model.subproblem,
                                      index_set=buses_with_loads)
    subcons.disjunctify(model=model.subproblem,
                        indicator_name='load_shed_indicator',
                        disjunct_name='load_shed_disjunct',
                        LHS_disjunct_set=model.subproblem.ineq_load_shed_lb,
                        RHS_disjunct_set=model.subproblem.ineq_load_shed_lb_off
                        )

    ### declare the generator disjuncts (LHS is status quo, RHS is compromised)
    subcons.declare_ineq_gen_on(model=model.subproblem,
                             index_set=gen_attrs['names'],
                             gens=gens)
    subcons.declare_ineq_gen_off(model=model.subproblem,
                                 index_set=gen_attrs['names'],
                                 gens=gens)
    subcons.disjunctify(model=model.subproblem,
                        indicator_name='gen_indicator',
                        disjunct_name='gen_disjunct',
                        LHS_disjunct_set=model.subproblem.ineq_gen,
                        RHS_disjunct_set=model.subproblem.ineq_gen_off
                        )

    ### declare the p balance
    rhs_kwargs = {'include_feasibility_slack_neg':'load_shed'}
    libbus.declare_eq_p_balance_dc_approx(model=model.subproblem,
                                          index_set=bus_attrs['names'],
                                          bus_p_loads=bus_p_loads,
                                          gens_by_bus=gens_by_bus,
                                          bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                          inlet_branches_by_bus=inlet_branches_by_bus,
                                          outlet_branches_by_bus=outlet_branches_by_bus,
                                          approximation_type=ApproximationType.BTHETA,
                                          **rhs_kwargs
                                          )

    ### declare the real power flow limits
    libbranch.declare_ineq_p_branch_thermal_lbub(model=model.subproblem,
                                                 index_set=branch_attrs['names'],
                                                 branches=branches,
                                                 p_thermal_limits=p_max,
                                                 approximation_type=ApproximationType.BTHETA
                                                 )

    ### declare angle difference limits on interconnected buses
    if include_angle_diff_limits:
        libbranch.declare_ineq_angle_diff_branch_lbub(model=model.subproblem,
                                                      index_set=branch_attrs['names'],
                                                      branches=branches,
                                                      coordinate_type=CoordinateType.POLAR
                                                      )

    model.subproblem.obj = pe.Objective(expr=sum(model.load_shed[l] for l in buses_with_loads), sense=pe.minimize)

    return model, md
