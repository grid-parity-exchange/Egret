#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

from __future__ import annotations

from typing import NamedTuple, Optional, Sequence

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

class ScalarReserveValue(NamedTuple):
    ''' A reserve type, scope, and scalar value.

    If area_name is None, this is a global reserve value.
    '''
    reserve_type: str
    area_name: Optional[str]
    value: float

class ScalarReserveData():
    ''' Scalar reserve values that should only be applied to one type of model.
    '''
    def __init__(self, 
                 da_scalars: Sequence[ScalarReserveValue], 
                 rt_scalars: Sequence[ScalarReserveValue]):
        self._da_scalars = da_scalars
        self._rt_scalars = rt_scalars

    @property
    def da_scalars(self) -> Sequence[ScalarReserveValue]:
        ''' Scalar reserve values that only apply to DAY_AHEAD models
        '''
        return self._da_scalars

    @property
    def rt_scalars(self) -> Sequence[ScalarReserveValue]:
        ''' Scalar reserve values that only apply to REAL_TIME models
        '''
        return self._rt_scalars

    def get_simulation_reserves(self, simulation_type:str) -> Sequence[ScalarReserveValue]:
        ''' Get scalar reserve values that only apply the requested simulation type
        '''
        return self._rt_scalars if simulation_type == 'REAL_TIME' else self._da_scalars
