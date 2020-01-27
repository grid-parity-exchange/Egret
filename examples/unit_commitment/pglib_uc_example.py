

## Example of solving a unit commitment problem from pglib-uc
import os

from egret.data.model_data import ModelData
from egret.models.unit_commitment import solve_unit_commitment


this_module_path = os.path.dirname(os.path.abspath(__file__))
if not os.path.isdir(os.path.join(this_module_path, '..', '..', 'download', 'pglib-uc-master')):
    print('Downloading pglib-uc...')
    from egret.thirdparty.get_pglib_uc import get_pglib_uc
    get_pglib_uc(download_dir=os.path.join(this_module_path,'..','..', 'download'))
    print('Downloaded pglib-uc!')


## Create an Egret 'ModelData' object, which is just a lightweight
## wrapper around a python dictionary, from a pglib-uc instance
print('Creating and solving rts_gmlc/2020-01-27.json ...')
md = ModelData.read(os.path.join(this_module_path, '..', '..', 'download', 
                                 'pglib-uc-master', 'rts_gmlc', '2020-01-27.json'),
                    file_type = 'pglib-uc')

## solve the unit commitment instance using solver cbc -- could use 'gurobi', 'cplex',
## or any valid Pyomo solver name, provided its available
md_sol = solve_unit_commitment(md, 'cbc', mipgap=0.01, timelimit=300, solver_tee=True)
print('Solved!')

## print the objective value to the screen
print('Objective value:', md_sol.data['system']['total_cost'])

## write the solution to an Egret *.json file
md_sol.write(os.path.join(this_module_path, 'rts_2020-01-27_solution.json'))
print('Wrote solution to rts_2020-01-27_solution.json')
