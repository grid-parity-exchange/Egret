#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module contains several helper functions that are useful when
modifying the data dictionary
"""
import egret.model_library.transmission.tx_calc as tx_calc
from egret.model_library.defn import BasePointType
import numpy as np

def create_dicts_of_ptdf(md,base_point=BasePointType.FLATSTART):
    branches = dict(md.elements(element_type='branch'))
    buses = dict(md.elements(element_type='bus'))
    branch_attrs = md.attributes(element_type='branch')
    bus_attrs = md.attributes(element_type='bus')

    reference_bus = md.data['system']['reference_bus']
    ptdf = tx_calc.calculate_ptdf(branches,buses,branch_attrs['names'],bus_attrs['names'],reference_bus,base_point)

    _len_bus = len(bus_attrs['names'])
    _len_branch = len(branch_attrs['names'])
    _mapping_branch = {i: branch_attrs['names'][i] for i in list(range(0,_len_branch))}

    for idx,branch_name in _mapping_branch.items():
        branch = md.data['elements']['branch'][branch_name]
        _row_ptdf = {bus_attrs['names'][i]: ptdf[idx,i] for i in list(range(0,_len_bus))}
        branch['ptdf'] = _row_ptdf

