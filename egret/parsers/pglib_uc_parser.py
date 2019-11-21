#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module provides supporting functions for parsing (and writing) input files
from the pglib-uc instances

.. todo::
    documentation and examples
"""

import os.path
import logging
import warnings
import json
import egret.data.model_data as md
import numpy as np

logger = logging.getLogger('egret.parsers.pglib_uc_parser')

def create_ModelData(pglib_uc_filename):
    """
    Parse a pglib-uc input file into a ModelData object containing
    the model_data dictionary

    Parameters
    ----------
    pglib_uc_filename : str
        Path and filename of the pglib-uc input file you wish to load

    Returns
    -------
        ModelData
    """
    data = create_model_data_dict(pglib_uc_filename)
    return md.ModelData(data)


def create_model_data_dict(pglib_uc_filename):
    """
    Parse a pglib-uc input file into a model_data dictionary

    Parameters
    ----------
    pglib_uc_filename : str
        Path and filename of the pglib-uc input file you wish to load

    Returns
    -------
        dict : Returns a dictionary in the format required for the ModelData
               object.
    """

    # create the model data object
    pglib_uc_dict = json.load(open(pglib_uc_filename, 'r'))

    model_data = md.ModelData.empty_model_data_dict()
    system = model_data["system"]
    elements = model_data["elements"]

    time_periods = list(range(1, pglib_uc_dict["time_periods"]+1))

    system["reserve_requirement"] = { "data_type" : "time_series", "values" : pglib_uc_dict["reserves"] }
    system["time_keys"] = time_periods
    system["baseMVA"] = 1.
    system["time_period_length_minutes"] = 60
    system["load_mismatch_cost"] = 10000
    system["reserve_shortfall_cost"] = 1000
    system["reference_bus"] = "copperplate"
    system["reference_bus_angle"] = 0.

    elements["bus"] = {"copperplate": dict()}
    
    elements["load"] = {"demand": { "bus" : "copperplate",
                                        "in_service" : True,
                                        "p_load" : {"data_type" : "time_series", "values" : pglib_uc_dict["demand"] }
                                      }
                       }

    elements["branch"] = dict()
    elements["zone"] = dict()

    generators = dict()

    name_mapping = {
                    "p_min" : "power_output_minimum",
                    "p_max" : "power_output_maximum",
                    "ramp_up_60min" : "ramp_up_limit",
                    "ramp_down_60min" : "ramp_down_limit",
                    "startup_capacity" : "ramp_startup_limit",
                    "shutdown_capacity" : "ramp_shutdown_limit",
                    "min_up_time" : "time_up_minimum",
                    "min_down_time" : "time_down_minimum",
                    "initial_p_output" : "power_output_t0",
                   }

    for g,g_dict in pglib_uc_dict["thermal_generators"].items():
        md_g_dict = {"generator_type" : "thermal", "bus" : "copperplate", "fuel" : "G", "in_service" : True}
        for mdn, pgucn in name_mapping.items():
            md_g_dict[mdn] = g_dict[pgucn]
        unit_on_t0 = g_dict["unit_on_t0"]
        md_g_dict["initial_status"] = unit_on_t0*g_dict["time_up_t0"] - (1-unit_on_t0)*g_dict["time_down_t0"]
        md_g_dict["fixed_commitment"] = 1 if bool(g_dict["must_run"]) else None
        md_g_dict["startup_cost"] = list( (s["lag"], s["cost"]) for s in g_dict["startup"] )
        md_g_dict["p_cost"] = { "data_type":"cost_curve", "cost_curve_type":"piecewise", 
                                "values": list( (c["mw"], c["cost"]) for c in g_dict["piecewise_production"]) }

        ## avoid name conflicts
        generators[g_dict["name"]+"_T"] = md_g_dict

    for g,g_dict in pglib_uc_dict["renewable_generators"].items():
        md_g_dict = {"generator_type" : "renewable", "bus" : "copperplate", "fuel" : "W", "in_service" : True}
        md_g_dict["p_min"] = {"data_type" : "time_series", "values" : g_dict["power_output_minimum"] }
        md_g_dict["p_max"] = {"data_type" : "time_series", "values" : g_dict["power_output_maximum"] }

        ## avoid name conflicts
        generators[g_dict["name"]+"_R"] = md_g_dict

    elements["generator"] = generators

    return model_data
