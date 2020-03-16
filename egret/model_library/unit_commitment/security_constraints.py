#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

## for adding generic constraints to the unit commitment model

from pyomo.environ import *
from .uc_utils import add_model_attr, uc_time_helper
from egret.data.data_utils import _recurse_into_time_index


@add_model_attr('security_constraints', requires = {'data_loader' : None,
                                                    'status_vars' : None,
                                                    'power_vars' : None,
                                                    'storage_service' : None,
                                                    'power_balance' : None,
                                                    })
def security_constraint_model(model):
    '''
    Add generic security constraints to the unit commitment model of 
    the form lb <= \sum_{g \in Gens} \\alpha_g (real power generation) \
                 + \sum_{s \in Storage} \\beta_s (net storage injection) <= ub, 
    with optional slack
    '''
    md = model.model_data
    time_keys = md.data['system']['time_keys']

    if not hasattr(model, 'TransmissionBlock'):
        model.TransmissionBlock = Block(model.TimePeriods, concrete=True)

    model.SecurityConstraintViolationCost = Expression(model.TimePeriods)

    thermal_generators = set(model.ThermalGenerators)
    renewable_generators = set(model.AllNondispatchableGenerators)

    pg_security_constraints = md.attributes(element_type='security_constraint', security_constraint_type='pg')
    pg_defined_keys = { 'generator', 'storage'}

    for i,t in enumerate(model.TimePeriods):
        b = model.TransmissionBlock[t]

        b_pg_sec_constr = _recurse_into_time_index(pg_security_constraints, i)
        needed_keys = ['violation_penalty', 'lower_bound', 'upper_bound']
        for nk in needed_keys:
            if nk not in b_pg_sec_constr:
                b_pg_sec_constr[nk] = dict()
            b_pg_sc_nk = b_pg_sec_constr[nk]
            for s in b_pg_sec_constr['names']:
                if s in b_pg_sc_nk and b_pg_sc_nk[s] is None:
                    del b_pg_sec_constr[nk][s]

        violation_penalty, lower_bound, upper_bound = (b_pg_sec_constr[nk] for nk in needed_keys)

        b.pgSecuritySet = Set(initialize=b_pg_sec_constr['names'])

        b.pgRelaxedSecuritySet = Set(within=b.pgSecuritySet, initialize=violation_penalty.keys())
        b.pgSecurityPenalty = Param(b.pgRelaxedSecuritySet, initialize=violation_penalty)

        ## check the lower_bound and upper_bound
        for s in b_pg_sec_constr['names']:
            if s in lower_bound:
                continue
            if s in upper_bound:
                continue
            logger.warning('Security constraint {} found to have no lower_bound and no upper_bound at time {}'.format(s, time_keys[i]))

        b.pgSecuritySlackPos = Var(b.pgRelaxedSecuritySet, within=NonNegativeReals)
        b.pgSecuritySlackNeg = Var(b.pgRelaxedSecuritySet, within=NonNegativeReals)

        b.pgSecurityExpression = Expression(b.pgSecuritySet)
        b.pgSecurityConstraint = Constraint(b.pgSecuritySet)

        for s, coefs in b_pg_sec_constr['coefs'].items():
            assert coefs['data_type'] == 'coefs'
            coefs = coefs['values']

            if not set(coefs.keys()) <= pg_defined_keys:
                raise Exception('Unrecognized or unsupported element type in coefs for security constraint {} at time {}'.format(s, time_keys[i]))
            for k in pg_defined_keys:
                if k not in coefs:
                    coefs[k] = dict()
            
            generator_coefs = coefs['generator']
            thermal_gen_coefs = dict()
            renewable_gen_coefs = dict()
            for g, val in generator_coefs.items():
                if g in thermal_generators:
                    thermal_gen_coefs[g] = val
                elif g in renewable_generators:
                    renewable_gen_coefs[g] = val
                else:
                    raise Exception('Unrecognized generator {} for security constraint {} at time {}'.format(g,s,time_keys[i]))

            ## TODO: Make LinearExpression?
            b.pgSecurityExpression[s] = sum( val*model.PowerGenerated[g,t] for g,val in thermal_gen_coefs.items() ) \
                                            + sum( val*model.NondispatchablePowerUsed[g,t] for g,val in renewable_gen_coefs.items() ) \
                                            + sum( val*(model.PowerOutputStorage[s,t]-model.PowerInputStorage[s,t]) for s,val in coefs['storage'].items())

            lb = lower_bound.get(s)
            ub = upper_bound.get(s)
            if lb is None and ub is None:
                b.pgSecurityConstraint[s] = Constraint.Feasible
                continue

            constr_expr = b.pgSecurityExpression[s]
            if s in b.pgRelaxedSecuritySet:
                constr_expr += b.pgSecuritySlackNeg[s]
                constr_expr -= b.pgSecuritySlackPos[s]

            b.pgSecurityConstraint[s] = (lb, constr_expr, ub)

        model.SecurityConstraintViolationCost[t] = sum(model.TimePeriodLengthHours*b.pgSecurityPenalty[s]*( \
                                                         b.pgSecuritySlackNeg[s] + b.pgSecuritySlackPos[s] )
                                                         for s in b.pgRelaxedSecuritySet)
