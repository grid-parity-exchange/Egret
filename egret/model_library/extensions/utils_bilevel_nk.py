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
working with the attacker-defender bilevel model
"""

def dict_of_relay_branches(relays, branches):
    """
    Create dictionaries of the branch keys from the
    relay elements
    """
    relay_branches = {k: list() for k in relays.keys()}

    for branch_name, branch in branches.items():
        branch_relay_mappings = branch['relay']
        for r in branch_relay_mappings:
            relay_branches[r].append(branch_name)

    return relay_branches


def dict_of_branch_relays(relays, branches):
    """
    Create dictionaries of the relay keys from the
    branch elements
    """
    branch_relays = {k: list() for k in branches.keys()}

    for relay_name, relay in relays.items():
        branch_relay_mappings = relay['branch']
        if branch_relay_mappings:
            for b in branch_relay_mappings:
                branch_relays[b].append(relay_name)

    return branch_relays


def relay_branch_tuple(relay_branches):
    """
    Create list of (relay, branch) tuples
    """
    relay_branch_tuple = list()

    for r, branches in relay_branches.items():
        for b in branches:
            relay_branch_tuple.append((r,b))
    return relay_branch_tuple


def dict_of_relay_gens(relays, gens):
    """
    Create dictionaries of the gen keys from the
    relay elements
    """
    relay_gens = {k: list() for k in relays.keys()}

    for gen_name, gen in gens.items():
        gen_relay_mappings = gen['relay']
        for r in gen_relay_mappings:
            relay_gens[r].append(gen_name)

    return relay_gens


def dict_of_gen_relays(relays, gens):
    """
    Create dictionaries of the relay keys from the
    gen elements
    """
    gen_relays = {k: list() for k in gens.keys()}

    for relay_name, relay in relays.items():
        gen_relay_mappings = relay['gen']
        if gen_relay_mappings:
            for b in gen_relay_mappings:
                gen_relays[b].append(relay_name)

    return gen_relays


def relay_gen_tuple(relay_gens):
    """
    Create list of (relay, gen) tuples
    """
    relay_gen_tuple = list()

    for r, gens in relay_gens.items():
        for g in gens:
            relay_gen_tuple.append((r,g))
    return relay_gen_tuple


def dict_of_relay_loads(relays, loads, buses_with_loads):
    """
    Create dictionaries of the load keys from the
    relay elements
    """
    relay_loads = {k: list() for k in relays.keys()}

    for load_name, load in loads.items():
        if load['bus'] in buses_with_loads:
            load_relay_mappings = load['relay']
            for r in load_relay_mappings:
                relay_loads[r].append(load['bus'])

    return relay_loads


def dict_of_load_relays(relays, buses_with_loads):
    """
    Create dictionaries of the relay keys from the
    load elements
    """
    load_relays = {k: list() for k in buses_with_loads}

    for relay_name, relay in relays.items():
        load_relay_mappings = relay['load']
        if load_relay_mappings:
            for b in load_relay_mappings:
                load_relays[b].append(relay_name)

    return load_relays


def relay_load_tuple(relay_loads):
    """
    Create list of (relay, load) tuples
    """
    relay_load_tuple = list()

    for r, loads in relay_loads.items():
        for l in loads:
            relay_load_tuple.append((r,l))
    return relay_load_tuple
