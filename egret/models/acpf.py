#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module provides functions that create the modules for typical ACPF formulations.
#TODO: document this with examples
"""
import pyomo.environ as pe
import egret.model_library.decl as decl
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen
from egret.data.data_utils import zip_items
from egret.model_library.defn import CoordinateType
from math import pi, radians
from collections import OrderedDict

def _create_base_acpf_model(model_data):
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
    inlet_branches_by_bus, outlet_branches_by_bus = \
        tx_utils.inlet_outlet_branches_by_bus(branches, buses)
    gens_by_bus = tx_utils.gens_by_bus(buses, gens)

    bus_pairs = zip_items(branch_attrs['from_bus'], branch_attrs['to_bus'])
    unique_bus_pairs = list(OrderedDict((val, None) for idx, val in bus_pairs.items()))

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
                            initialize={k: v**2 for k, v in bus_attrs['vm'].items()}
                            )
    libbranch.declare_var_c(model=model, index_set=unique_bus_pairs)
    libbranch.declare_var_s(model=model, index_set=unique_bus_pairs)

    ### declare the generator real and reactive power
    #pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model, gen_attrs['names'], initialize=gen_attrs['pg'])

    #qg_init = {k: (gen_attrs['q_min'][k] + gen_attrs['q_max'][k]) / 2.0 for k in gen_attrs['qg']}
    libgen.declare_var_qg(model, gen_attrs['names'], initialize=gen_attrs['qg'])

    ### declare the current flows in the branches
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    vj_init = {k: bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]) for k in bus_attrs['vm']}
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
                             initialize=pf_init
                             )
    libbranch.declare_var_pt(model=model,
                             index_set=branch_attrs['names'],
                             initialize=pt_init
                             )
    libbranch.declare_var_qf(model=model,
                             index_set=branch_attrs['names'],
                             initialize=qf_init
                             )
    libbranch.declare_var_qt(model=model,
                             index_set=branch_attrs['names'],
                             initialize=qt_init
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
                                outlet_branches_by_bus=outlet_branches_by_bus
                                )

    libbus.declare_eq_q_balance(model=model,
                                index_set=bus_attrs['names'],
                                bus_q_loads=bus_q_loads,
                                gens_by_bus=gens_by_bus,
                                bus_bs_fixed_shunts=bus_bs_fixed_shunts,
                                inlet_branches_by_bus=inlet_branches_by_bus,
                                outlet_branches_by_bus=outlet_branches_by_bus
                                )

    # if there are multiple generators at the same bus, we will
    # have unwanted degrees of freedom in qg
    # therefore, we add a constraint making them equal
    # if the reference bus has multiple generators, we will also
    # have unwanted degrees of freedom in pg
    #ref_bus = md.data['system']['reference_bus']
    qg_equality_tuples = list()
    #pg_equality_tuples = list()
    for b, genlist in gens_by_bus.items():
        if len(genlist) > 1:
            # we have more than one generator at this bus
            for i in range(1,len(genlist)):
                qg_equality_tuples.append((genlist[0], genlist[i]))
    #        if b == ref_bus:
    #            pg_equality_tuples.append((genlist[0], genlist[i]))

    def _qg_equalities(m,i,j):
        return m.qg[i] == m.qg[j]
    model.qg_equalities = pe.Constraint(qg_equality_tuples, rule=_qg_equalities)

    #def _pg_equalities(m,i,j):
    #    return m.pg[i] == m.pg[j]
    #model.pg_equalities = pe.Constraint(pg_equality_tuples, rule=_pg_equalities)

    model.obj = pe.Objective(expr=0.0)
    return model, md

def create_psv_acpf_model(model_data):
    model, md = _create_base_acpf_model(model_data)

    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    bus_attrs = md.attributes(element_type='bus')
    branch_attrs = md.attributes(element_type='branch')
    gens_by_bus = tx_utils.gens_by_bus(buses, gens)
    buses_with_gens = _buses_with_gens(gens)
    bus_pairs = zip_items(branch_attrs['from_bus'], branch_attrs['to_bus'])
    unique_bus_pairs = list(OrderedDict((val, None) for idx, val in bus_pairs.items()).keys())

    # declare the polar voltages
    libbranch.declare_var_dva(model=model,
                              index_set=unique_bus_pairs,
                              initialize=0
                              )

    libbus.declare_var_vm(model,
                          bus_attrs['names'],
                          initialize=bus_attrs['vm']
                          )

    libbus.declare_var_va(model,
                          bus_attrs['names'],
                          initialize=bus_attrs['va']
                          )

    ### In a system with N buses and G generators, there are then 2(N-1)-(G-1) unknowns.
    ### fix the reference bus
    ref_bus = md.data['system']['reference_bus']
    ref_angle = md.data['system']['reference_bus_angle']
    model.va[ref_bus].fix(radians(ref_angle))
    model.vm[ref_bus].fixed = True

    # if there is more than one generator at the reference
    # bus, then we fix the pg for all but one
    for i,g in enumerate(gens_by_bus[ref_bus]):
        if i > 0:
            model.pg[g].fixed = True

    for bus_name in bus_attrs['names']:
        if bus_name != ref_bus and bus_name in buses_with_gens:
            model.vm[bus_name].fixed = True
            for gen_name in gens_by_bus[bus_name]:
                model.pg[gen_name].fixed = True

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


def solve_acpf(model_data,
                solver,
                timelimit = None,
                solver_tee = True,
                symbolic_solver_labels = False,
                options = None,
                acpf_model_generator = create_psv_acpf_model,
                return_model = False,
                return_results = False,
                **kwargs):
    '''
    Create and solve a new acpf model
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
    acpf_model_generator : function (optional)
        Function for generating the acpf model. Default is
        egret.models.acpf.create_psv_acpf_model
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

    m, md = acpf_model_generator(model_data, **kwargs)

    m.dual = pe.Suffix(direction=pe.Suffix.IMPORT)

    m, results = _solve_model(m,solver,timelimit=timelimit,solver_tee=solver_tee,
                              symbolic_solver_labels=symbolic_solver_labels,solver_options=options)

    if results is not None:
        # save results data to ModelData object
        gens = dict(md.elements(element_type='generator'))
        buses = dict(md.elements(element_type='bus'))
        branches = dict(md.elements(element_type='branch'))

        md.data['system']['total_cost'] = value(m.obj)

        for g,g_dict in gens.items():
            g_dict['pg'] = value(m.pg[g])
            g_dict['qg'] = value(m.qg[g])

        for b,b_dict in buses.items():
            b_dict['pl'] = value(m.pl[b])
            if hasattr(m,'p_slack_pos'):
                b_dict['p_slack_pos'] = value(m.p_slack_pos[b])
            if hasattr(m, 'p_slack_neg'):
                b_dict['p_slack_neg'] = value(m.p_slack_neg[b])
            if hasattr(m, 'q_slack_pos'):
                b_dict['q_slack_pos'] = value(m.q_slack_pos[b])
            if hasattr(m, 'q_slack_neg'):
                b_dict['q_slack_neg'] = value(m.q_slack_neg[b])
            if hasattr(m, 'vj'):
                b_dict['vm'] = tx_calc.calculate_vm_from_vj_vr(value(m.vj[b]), value(m.vr[b]))
                b_dict['va'] = tx_calc.calculate_va_from_vj_vr(value(m.vj[b]), value(m.vr[b]))
            else:
                b_dict['vm'] = value(m.vm[b])
                b_dict['va'] = value(m.va[b])

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

# Todo: This should be moved to tx_utils following initial commit
def _buses_with_gens(gens):
    """
    Return a list of buses with generators
    """
    buses_with_gens = list()
    for gen_name, gen in gens.items():
        if not gen['bus'] in buses_with_gens:
            buses_with_gens.append(gen['bus'])

    return buses_with_gens
