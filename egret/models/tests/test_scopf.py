#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
scopf tester
'''

import os
import math
import unittest

import pyomo.environ as pyo

from egret.models.scopf import solve_scopf
from egret.models.dcopf import create_btheta_dcopf_model

from egret.data.model_data import ModelData

current_dir = os.path.dirname(os.path.abspath(__file__))

case_300_fn = os.path.join(current_dir, 'transmission_test_instances', 'pglib-opf-master', 'pglib_opf_case300_ieee.m')

class TestSCOPF(unittest.TestCase):
    show_output = True

    def test_scopf_model(self):
        md = ModelData.read(case_300_fn)
        
        contingencies = ['390', '150', '100', '240', '382']
        
        contingencies = { cn : {'branch_contingency': cn } for cn in contingencies }
        
        md.data['elements']['contingency'] = contingencies
        
        mdos, pyos = solve_scopf(md, 'cbc', return_model=True)

        # spot checks on mdos
        assert math.isclose( mdos.data['system']['total_cost'], 517915.03835644305 )
        assert math.isclose( mdos.data['elements']['contingency']['390']['monitored_branches']['137']['pf'], -815. )
        assert math.isclose( mdos.data['elements']['bus']['664']['lmp'], 39.69733343779804 )
        
        PTDF = pyos._PTDF
        PFV, _, VA = PTDF.calculate_PFV(pyos)
        
        ## formulate dcopf without this branch
        mdbo = md
        #del mdbo.data['elements']['branch'][branch_out]
        pyobo, _ = create_btheta_dcopf_model(mdbo)
        
        ## set generation levels equal to security-constraint solution
        for g in pyos.pg:
            pyobo.pg[g].value = pyo.value(pyos.pg[g])
            pyobo.pg[g].fix()
        
        ## turn off lower/upper bounds
        pyobo.ineq_pf_branch_thermal_lb.deactivate()
        pyobo.ineq_pf_branch_thermal_ub.deactivate()
        pyobo.pf.setlb(None)
        pyobo.pf.setub(None)
        pyobo.obj.expr = 0.
        
        s = pyo.SolverFactory('cbc')
        
        for cn, cont in contingencies.items():
            branch_out = cont['branch_contingency']
        
            pf_branch_con = pyobo.eq_pf_branch[branch_out]
            pf_branch_con.deactivate()
            
            pf_branch_var = pyobo.pf[branch_out]
            pf_branch_var.value = 0.
            pf_branch_var.fix()
            
            s.solve(pyobo, tee=False)
            
            ## check and see if the flows are the same
            for (c_all, bn) in pyos.ineq_pf_contingency_branch_thermal_bounds:
                if c_all == cn:
                    print(f"checking model flows {(cn, bn)}")
                    if not math.isclose( pyo.value(pyos.pfc[cn,bn]), pyo.value(pyobo.pf[bn]), rel_tol=1e-6, abs_tol=1e-10):
                        print(f"pyos.pfc[cn,bn] : {pyo.value(pyos.pfc[cn,bn])}, pyobo.pf[bn] : {pyo.value(pyobo.pf[bn])}")
                        assert False
        
            PFV_delta = PTDF.calculate_PFV_delta(cn, PFV, VA)
            PFC = PFV + PFV_delta
            print(f"checking measured flows {cn}")
            for i, bn in enumerate(PTDF.branches_keys):
                if not math.isclose( PFC[i], pyo.value(pyobo.pf[bn]),rel_tol=1e-6, abs_tol=1e-10 ):
                    print(f"PFC[i] : {PFC[i]}, pyobo.pf[bn] : {pyo.value(pyobo.pf[bn])}")
                    assert False
        
            #print(f"FLOWS EQUAL for contingency {cn}")
        
            ## reset
            pf_branch_con.activate()
        
            pf_branch_var.unfix()
