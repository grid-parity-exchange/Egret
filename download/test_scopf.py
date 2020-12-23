import sys
import math

import pyomo.environ as pyo
import egret.model_library.transmission.tx_calc as tx_calc

from egret.models.scopf import solve_scopf
from egret.models.dcopf import create_btheta_dcopf_model

from egret.data.model_data import ModelData

tee=True
md = ModelData.read(sys.argv[1])

#branch_out = sys.argv[2]
mapping_bus_to_idx = { b : idx for idx,b in enumerate(md.data['elements']['bus']) }
branches = {b: bd for b,bd in md.data['elements']['branch'].items() if bd['in_service']}
graph = tx_calc.construct_connection_graph(branches, mapping_bus_to_idx)
branches_not_disconnecting = tx_calc.get_N_minus_1_branches(graph, branches, mapping_bus_to_idx)

contingencies = { bn : {'branch_contingency': bn } for bn in branches_not_disconnecting }
#contingencies = {}
print(f"Number of contingencies: {len(contingencies)}")

md.data['elements']['contingency'] = contingencies

mdos, pyos = solve_scopf(md, 'xpress_persistent', solver_tee=tee, return_model=True)

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

s = pyo.SolverFactory('xpress_direct')
#s.solve(pyobo, tee=True)
#s.set_instance(pyobo)

for cn, cont in contingencies.items():
    branch_out = cont['branch_contingency']

    pf_branch_con = pyobo.eq_pf_branch[branch_out]
    pf_branch_con.deactivate()
    #s.remove_constraint(pf_branch_con)
    
    pf_branch_var = pyobo.pf[branch_out]
    pf_branch_var.value = 0.
    pf_branch_var.fix()
    #s.update_var(pf_branch_var)
    
    s.solve(pyobo, tee=False)
    
    ## check and see if the flows are the same
    for (c_all, bn) in pyos.ineq_pf_contingency_branch_thermal_bounds:
        if c_all == cn:
            print(f"checking {(cn, bn)}")
            assert math.isclose( pyo.value(pyos.pfc[cn,bn]), pyo.value(pyobo.pf[bn]) )
    
    #print(f"FLOWS EQUAL for contingency {cn}")

    ## reset
    pf_branch_con.activate()
    #s.add_constraint(pf_branch_con)

    pf_branch_var.unfix()
    #s.update_var(pf_branch_var)

print("CONTINGENCY FLOWS VERIFIED")
