#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module contains several constraints that are useful when
working with the attacker-defender bilevel model
"""

import egret.model_library.decl as decl
import pyomo.environ as pe


def declare_budget(model, k, relays):
    """
    Create the budget constraint for attacks
    """
    m = model

    m.budget = pe.Constraint(expr=sum(m.delta[r] for r in relays) <= k)

def declare_physical_budget(model, k):
    """
    Create the budget constraint for physical attacks
    """

    m = model
    m.budget = pe.Constraint(expr=sum(5*m.delta_gen[g] for g in m.delta_gen.index_set()) +\
        sum(1*m.delta_branch[k] for k in m.delta_branch.index_set()) +\
        sum(3*m.delta_load[b] for b in m.delta_load.index_set()) +\
        sum(1*m.delta_bus[b] for b in m.delta_bus.index_set()) == k)

def declare_load_compromised(model, index_set):

    m = model
    con_set = decl.declare_set("_con_load_compromised", model=model, index_set=index_set)

    m.load_compromised = pe.Constraint(con_set)

    for (r,l) in con_set:
        m.load_compromised[(r,l)] = m.u[l] <= 1 - m.delta[r]


def declare_load_uncompromised(model, index_set, load_relays):

    m = model
    con_set = decl.declare_set("_con_load_uncompromised", model=model, index_set=index_set)

    m.load_uncompromised = pe.Constraint(con_set)

    for l in con_set:
        m.load_uncompromised[l] = sum((1 - m.delta[r]) for r in load_relays[l]) - len(load_relays[l]) + 1 <= m.u[l]


def declare_branch_compromised(model, index_set):

    m = model
    con_set = decl.declare_set("_con_branch_compromised", model=model, index_set=index_set)

    m.branch_compromised = pe.Constraint(con_set)

    for (r,b) in con_set:
        m.branch_compromised[(r,b)] = m.w[b] <= 1 - m.delta[r]


def declare_branch_uncompromised(model, index_set, branch_relays):

    m = model
    con_set = decl.declare_set("_con_branch_uncompromised", model=model, index_set=index_set)

    m.branch_uncompromised = pe.Constraint(con_set)

    for b in con_set:
        m.branch_uncompromised[b] = sum((1 - m.delta[r]) for r in branch_relays[b]) - len(branch_relays[b]) + 1 <= m.w[b]


def declare_gen_compromised(model, index_set):

    m = model
    con_set = decl.declare_set("_con_gen_compromised", model=model, index_set=index_set)

    m.gen_compromised = pe.Constraint(con_set)

    for (r,g) in con_set:
        m.gen_compromised[(r,g)] = m.v[g] <= 1 - m.delta[r]


def declare_gen_uncompromised(model, index_set, gen_relays):

    m = model
    con_set = decl.declare_set("_con_gen_uncompromised", model=model, index_set=index_set)

    m.gen_uncompromised = pe.Constraint(con_set)

    for g in con_set:
        m.gen_uncompromised[g] = sum((1 - m.delta[r]) for r in gen_relays[g]) - len(gen_relays[g]) + 1 <= m.v[g]

def declare_physical_load_compromised(model, index_set):

    m = model
    con_set = decl.declare_set("_con_physical_load_compromised", model=model, index_set=index_set)

    m.physical_load_compromised = pe.Constraint(con_set)
    m.physical_load_compromised_by_bus = pe.Constraint(con_set)
    
    for l in con_set:
        m.physical_load_compromised[l] = m.u[l] <= 1 - m.delta_load[l]
        m.physical_load_compromised_by_bus[l] = m.u[l] <= 1 - m.delta_bus[l]

def declare_physical_load_uncompromised(model, index_set):
    
    m = model
    con_set = decl.declare_set("_con_physical_load_uncompromised", model=model, index_set=index_set)

    m.physical_load_uncompromised = pe.Constraint(con_set)

    for l in con_set:
        m.physical_load_uncompromised[l] = m.u[l] >= 1 - m.delta_bus[l] - m.delta_load[l]

def declare_physical_gen_uncompromised(model, index_set, bus_by_gen):
    
    m = model
    con_set = decl.declare_set("_con_physical_gen_compromised", model=model, index_set=index_set)

    m.physical_gen_compromised = pe.Constraint(con_set)
    m.physical_gen_compromised_by_bus = pe.Constraint(con_set)

    for g in con_set:
        m.physical_gen_compromised[g] = m.v[g] <= 1 - m.delta_gen[g]
        m.physical_gen_compromised_by_bus[g] = m.v[g] <= 1 - m.delta_bus[bus_by_gen[g]]

def declare_physical_gen_compromised(model, index_set, bus_by_gen):

    m = model
    con_set = decl.declare_set("_con_physical_gen_uncompromised", model=model, index_set=index_set)

    m.physical_gen_uncompromised = pe.Constraint(con_set)

    for g in con_set:
        m.physical_gen_uncompromised[g] = m.v[g] >= 1 - m.delta_gen[g] - m.delta_bus[bus_by_gen[g]] 

def declare_physical_line_compromised(model, index_set, buses_by_branch):

    m = model
    con_set = decl.declare_set("_con_physical_line_compromised", model=model, index_set=index_set)

    m.physical_line_compromised = pe.Constraint(con_set)
    m.physical_line_compromised_by_from_bus = pe.Constraint(con_set)
    m.physical_line_compromised_by_to_bus = pe.Constraint(con_set)

    for k in con_set:
        m.physical_line_compromised[k] = m.w[k] <= 1 - m.delta_branch[k]
        m.physical_line_compromised_by_from_bus[k] = m.w[k] <= 1 - m.delta_bus[buses_by_branch[k][0]]
        m.physical_line_compromised_by_to_bus[k] = m.w[k] <= 1 - m.delta_bus[buses_by_branch[k][1]]

def declare_physical_line_uncompromised(model, index_set, buses_by_branch):

    m = model
    con_set = decl.declare_set("_con_physical_line_uncompromised", model=model, index_set=index_set)

    m.physical_line_uncompromised = pe.Constraint(con_set)

    for k in con_set:
        m.physical_line_uncompromised[k] = m.w[k] >= 1 - m.delta_branch[k] - m.delta_bus[buses_by_branch[k][0]] - m.delta_bus[buses_by_branch[k][1]]