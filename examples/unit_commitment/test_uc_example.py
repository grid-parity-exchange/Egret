

## Example of solving a unit commitment problem in egret format
## from the test instance library
import os

from egret.data.model_data import ModelData
from egret.models.unit_commitment import solve_unit_commitment

this_module_path = os.path.dirname(os.path.abspath(__file__))
## Create an Egret "ModelData" object, which is just a lightweight
## wrapper around a python dictionary, from an Egret json test instance
print('Creating and solving tiny_uc_tc')
md = ModelData.read(os.path.join(this_module_path, '..','..','egret','models',
                                 'tests','uc_test_instances','tiny_uc_tc.json'))

## solve the unit commitment instance using solver cbc
md_sol = solve_unit_commitment(md, 'cbc', mipgap=0.01, timelimit=300, solver_tee=True)
print('Solved!')

## print the objective value to the screen
print('Objective value:', md_sol.data['system']['total_cost'])

## write the solution to an Egret *.json file
md_sol.write(os.path.join(this_module_path, 'tiny_uc_tc_solution.json'))
print('Wrote solution to tiny_uc_tc_solution.json')
