#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## functions for adding the basic status varibles
from pyomo.environ import *
import math

from .uc_utils import add_model_attr 
component_name = 'status_vars'

def _is_relaxed(model):
    if hasattr(model, 'relax_binaries') and model.relax_binaries:
        return True
    else:
        return False

def _add_unit_on_vars(model, relaxed=False):
    # indicator variables for each generator, at each time period.
    if relaxed:
        model.UnitOn = Var(model.ThermalGenerators, model.TimePeriods, within=UnitInterval) 
    else:
        model.UnitOn = Var(model.ThermalGenerators, model.TimePeriods, within=Binary) 

def _add_unit_start_vars(model, relaxed=False):
    # unit start
    if relaxed:
        model.UnitStart=Var(model.ThermalGenerators,model.TimePeriods, within=UnitInterval)
    else:
        model.UnitStart=Var(model.ThermalGenerators,model.TimePeriods, within=Binary)

def _add_unit_stop_vars(model, relaxed=False):

    if relaxed:
        model.UnitStop=Var(model.ThermalGenerators,model.TimePeriods, within=UnitInterval)

    else:
        model.UnitStop=Var(model.ThermalGenerators,model.TimePeriods, within=Binary)


@add_model_attr(component_name, requires = {'data_loader': None} )
def CA_1bin_vars(model):
    '''
    This adds only a binary variable for unit-on, as in

    Carrion, M. and Arroyo, J. (2006) A Computationally Efficient Mixed-Integer
    Liner Formulation for the Thermal Unit Commitment Problem. IEEE Transactions
    on Power Systems, Vol. 21, No. 3, Aug 2006.
    '''
    if _is_relaxed(model):
        _add_unit_on_vars(model, True)
    else:
        _add_unit_on_vars(model)

@add_model_attr(component_name, requires = {'data_loader': None} )
def garver_3bin_vars(model):
    '''
    This add the common 3-binary variables per generator per time period.
    One for start, one for stop, and one for on, as originally proposed in

    L. L. Garver. Power generation scheduling by integer programming-development
    of theory. Power Apparatus and Systems, Part III. Transactions of the
    American Institute of Electrical Engineers, 81(3): 730–734, April 1962. ISSN
    0097-2460.

    '''

    if _is_relaxed(model):
        _add_unit_on_vars(model, True)
        _add_unit_start_vars(model, True)
        _add_unit_stop_vars(model, True)
    else:
        _add_unit_on_vars(model)
        _add_unit_start_vars(model)
        _add_unit_stop_vars(model)

    return

@add_model_attr(component_name, requires = {'data_loader': None} )
def garver_2bin_vars(model):
    '''
    This adds the unit start and unit on variables, and causes the
    unit stop variable to be projected out.
    '''

    if _is_relaxed(model):
        _add_unit_on_vars(model, True)
        _add_unit_start_vars(model, True)
    else:    
        _add_unit_on_vars(model)
        _add_unit_start_vars(model)
    
    # unit stop
    def unit_stop_expr_rule(m, g, t):
        if t == value(m.InitialTime):
            return m.UnitOnT0[g] - m.UnitOn[g,t] + m.UnitStart[g,t]
        return m.UnitOn[g,t-1] - m.UnitOn[g,t] + m.UnitStart[g,t] 
    model.UnitStop=Expression(model.ThermalGenerators,model.TimePeriods, rule=unit_stop_expr_rule)

    return


@add_model_attr(component_name, requires = {'data_loader': None} )
def garver_3bin_relaxed_stop_vars(model):
    '''
    This adds the 3-binary variables, but relaxes the integrality on the stop
    variable, like the "MILP-3R" formulation from

    Carrion, M. and Arroyo, J. (2006) A Computationally Efficient Mixed-Integer
    Liner Formulation for the Thermal Unit Commitment Problem. IEEE Transactions
    on Power Systems, Vol. 21, No. 3, Aug 2006.
    '''

    if _is_relaxed(model):
        _add_unit_on_vars(model, True)
        _add_unit_start_vars(model, True)
    else:    
        _add_unit_on_vars(model)
        _add_unit_start_vars(model)

    _add_unit_stop_vars(model, True)

    return


@add_model_attr(component_name, requires = {'data_loader': None} )
def ALS_state_transition_vars(model):
    '''
    These are the state-transition variables proposed in
    
    Atakan, Semih, Guglielmo Lulli, and Suvrajeet Sen. "A State Transition MIP
    Formulation for the Unit Commitment Problem." IEEE Transactions on Power
    Systems 33.1 (2018): 736-748.
    '''

    if _is_relaxed(model):
        model.UnitStayOn = Var(model.ThermalGenerators, model.TimePeriods, within=UnitInterval)
        _add_unit_start_vars(model, True)
        _add_unit_stop_vars(model, True)
    else:
        model.UnitStayOn = Var(model.ThermalGenerators, model.TimePeriods, within=Binary)
        _add_unit_start_vars(model)
        _add_unit_stop_vars(model)


    def unit_on_expr_rule(m, g, t):
        return m.UnitStayOn[g,t] + m.UnitStart[g,t]
    model.UnitOn = Expression(model.ThermalGenerators, model.TimePeriods, rule=unit_on_expr_rule)
