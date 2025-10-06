import os
import glob

import pytest

from egret.data.model_data import ModelData
from egret.models.scopf import solve_scopf
from egret.model_library.transmission.tx_calc import construct_connection_graph, get_N_minus_1_branches

# TODO: automatic solver detection
solver = "xpress_persistent"

GENERATE_BASELINE = False
ABS_TOL = 1e-6
NUMBER_PIECEWISE_POINTS = 4
LOAD_MISMATCH_COST = 1e+6

_this_directory = os.path.dirname(os.path.abspath(__file__))
_pglib_opf_dir = os.path.join(_this_directory, "../egret/thirdparty/pglib-opf-master")

_baseline_dir = os.path.join(_this_directory, "linearized_pglib_scopf_baseline")

# sort the matpower files by size we we're more likely to hit failure on small models
_all_matpower_files = sorted(glob.glob(_pglib_opf_dir + "/*.m"), key=os.path.getsize)

_ptdf_options = {
    "rel_ptdf_tol" : 0.0,
    "abs_ptdf_tol" : 0.0,
    "abs_flow_tol" : 0.0,
    "rel_flow_tol" : 0.0,
}

def _get_matpower_file_name(matpower_file):
    return os.path.basename(matpower_file)


class TestLinearizedSCOPF:

    @pytest.mark.parametrize("matpower_file", _all_matpower_files, ids=_get_matpower_file_name)
    def test_linearized_pglib_scopf(self, matpower_file):
    
        baseline_file = os.path.join(_baseline_dir, os.path.splitext(_get_matpower_file_name(matpower_file))[0]+".json")
        if GENERATE_BASELINE and os.path.exists(baseline_file):
            return
    
        try:
            md = ModelData.read(matpower_file)
        except ValueError as ve:
            # some PGLIB cases have disconnected buses, which Egret currently does not support
            if str(ve) == "Encountered an unsupported bus type: 4 when parsing MATPOWER input file":
                return pytest.skip("Egret is not compatible with this .m file")
            else:
                raise ve
    
        # only want in-service items
        md = md.clone_in_service()
    
        _linearize_objective(md)
    
        _add_nondisconnecting_contingencies(md)
    
        _add_load_mismatch_cost(md)
    
        try:
            mdo = solve_scopf(md, solver, include_feasibility_slack=True, ptdf_options=_ptdf_options)
        except ZeroDivisionError:
            # some PGLIB cases have 0 reactance ??
            return pytest.skip("Problem with input data; likely a line with 0 reactance")
    
        if GENERATE_BASELINE:
            mdo.write(baseline_file)
        else:
            baseline = ModelData.read(baseline_file)
            # TODO: primal degenercy can be an issue for this test
            # results_gens = mdo.data["elements"]["generator"]
            # for g, gd in baseline.elements("generator"):
            #     assert results_gens[g]["pg"] == pytest.approx( gd["pg"], abs=ABS_TOL )
            results_buses = mdo.data["elements"]["bus"]
            for b, bd in baseline.elements("bus"):
                assert results_buses[b]["lmp"] == pytest.approx( bd["lmp"], abs=ABS_TOL )


def _add_nondisconnecting_contingencies(md):
    mapping_bus_to_idx = { k : i for i,k in enumerate(md.data['elements']['bus'].keys())}
    graph = construct_connection_graph(md.data['elements']['branch'], mapping_bus_to_idx)
    contingency_list = get_N_minus_1_branches(graph, md.data['elements']['branch'], mapping_bus_to_idx)
    contingency_dict = { cn : {'branch_contingency':cn} for cn in contingency_list}

    md.data['elements']['contingency'] = contingency_dict


def _add_load_mismatch_cost(md):
    md.data["system"]["load_mismatch_cost"] = LOAD_MISMATCH_COST


def _linearize_objective(md):

    for _,gd in md.elements("generator"):
        p_cost = gd["p_cost"]
        if p_cost["cost_curve_type"] == "piecewise":
            continue
        assert p_cost["cost_curve_type"] == "polynomial"
        p_min = gd["p_min"]
        p_max = gd["p_max"]
        step_size = (p_max - p_min) / (NUMBER_PIECEWISE_POINTS-1)
        pnts = [p_min + step_size * i for i in range(NUMBER_PIECEWISE_POINTS)]
        pnts[-1] = p_max

        assert len(pnts) == NUMBER_PIECEWISE_POINTS

        vals = [ sum(coef*(pnt**exp) for exp, coef in p_cost["values"].items()) for pnt in pnts ]  

        p_cost["cost_curve_type"] = "piecewise"
        p_cost["values"] = [ (pnt, val) for pnt, val in zip(pnts, vals) ]
