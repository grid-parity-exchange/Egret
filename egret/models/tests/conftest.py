#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
pytest configuration options for test_unit_commitment.py,
per the pytest examples
'''
import pytest

def pytest_addoption(parser):
    parser.addoption("--runmip", action="store_true", default=False,
                     help="If enabled, this solves the MIP for each unit "
                          "commitment instance. Either the Gurobi or CPLEX "
                          " solver is required for this test."
                    )

def pytest_collection_modifyitems(config, items):
    if not config.getoption("--runmip"):
        skip_mip = pytest.mark.skip(reason="need --runmip option to run")
        for item in items:
            if "mip" in item.keywords:
                item.add_marker(skip_mip)
