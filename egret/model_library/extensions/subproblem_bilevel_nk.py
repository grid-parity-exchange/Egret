#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module contains several disjunctive subproblem constraints that
are useful when working with the attacker-defender bilevel model
"""

import egret.model_library.decl as decl
import pyomo.environ as pe
import pao.bilevel as bi
import pyomo.gdp as gdp
import egret.model_library.transmission.branch as libbranch


def declare_gdp_branch_power_btheta_approx(model, index_set, branches):
    """
    Create the budget constraint for attacks
    """
    def _callee_function(disjunct, idx):
        branch = branches[idx]
        bi.varref(disjunct)

        ### declare the branch power flow approximation constraints
        libbranch.declare_eq_branch_power_btheta_approx(disjunct, idx, branch)

    model.gdp_branch_power_btheta_approx = gdp.Disjunct(index_set, rule=_callee_function)




