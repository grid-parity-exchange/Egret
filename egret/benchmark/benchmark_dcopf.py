import sys

from egret.models.dcopf import (solve_dcopf,
        create_btheta_dcopf_model,
        create_ptdf_dcopf_model,
        )
from egret.data.model_data import ModelData

solver = "xpress_persistent"
if len(sys.argv) < 3:
    print("Useage: python benchmark_ptdf.py matpower_file dcopf_type")
    print("""dcopf_type is one of "btheta", "fullptdf", "virtualptdf" """)
    sys.exit(0)

matpower_file = sys.argv[1]
dcopf_type = sys.argv[2].lower().replace("-","").replace("_","")

md = ModelData.read(matpower_file)

if dcopf_type == "btheta":
    mdo = solve_dcopf(md, solver, dcopf_model_generator=create_btheta_dcopf_model) 
elif "ptdf" in dcopf_type:
    kwargs = {}
    if dcopf_type == "fullptdf":
        kwargs["virtual_ptdf"] = False
    elif dcopf_type == "virtualptdf":
        kwargs["virtual_ptdf"] = True
    else:
        raise RuntimeError(f"Unrecognized dcopf_type {dcopf_type}")
    mdo = solve_dcopf(md, solver, dcopf_model_generator=create_ptdf_dcopf_model, **kwargs)
else:
    raise RuntimeError(f"Unrecognized dcopf_type {dcopf_type}")
