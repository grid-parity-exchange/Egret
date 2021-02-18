#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

from __future__ import annotations

reserve_name_map = {
    'Spin_Up': 'spinning_reserve_requirement',
    'Reg_Up': 'regulation_up_requirement',
    'Reg_Down': 'regulation_down_requirement',
    'Flex_Up': 'flexible_ramp_up_requirement',
    'Flex_Down': 'flexible_ramp_down_requirement'
}

def is_valid_reserve_name(name:str, model_dict:dict=None):
    if name in reserve_name_map:
        return True
    if name.find('_R') < 0:
        return False
    res, area = name.split('_R', 1)
    return (res in reserve_name_map) and \
           ((model_dict is None) or (area in model_dict['elements']['area']))
