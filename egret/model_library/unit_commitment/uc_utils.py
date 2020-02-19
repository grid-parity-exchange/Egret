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
working with unit commitment models
"""

## some useful function decorators for building these dynamic models
from functools import wraps
import warnings

import logging
logger = logging.getLogger('egret.model_library.unit_commitment.uc_utils')

def add_model_attr(attr, requires = {}):
    def actual_decorator(func):
        @wraps(func)
        def wrapper(*args, **kwds):
            ## tag this function in the model with the appropriate attribute
            model = args[0]
            if hasattr(model, attr):
                msg = "Warning: adding %s! Model already has %s %s! You may only add one type of %s!"%(func.__name__, attr, getattr(model,attr), attr)
                logger.warning(msg)
                warnings.warn(msg)
            # this checks to see if the required components were already added
            for base_attr in requires:
                if not hasattr(model, base_attr):
                    msg = "Warning: adding %s! %s requires some %s to be added first!"%(func.__name__, func.__name__, base_attr)
                    logger.warning(msg)
                    warnings.warn(msg)
                ## None in this context means there is no specific requirement
                if requires[base_attr] is None:
                    continue
                if getattr(model, base_attr) not in requires[base_attr]:
                    msg = "Warning: adding %s! %s requires one of: "%(func.__name__, func.__name__) + ", ".join(requires[base_attr]) + ", to be added first."
                    logger.warning(msg)
                    warnings.warn(msg)
            setattr(model, attr, func.__name__)
            return func(*args, **kwds)
        return wrapper
    return actual_decorator

## provides a view on grid_data attributes that
## is handy for building pyomo params
## Assums the last key is time
def uc_time_helper(_data):
    ## if there is no data,
    ## we return None to the initializer
    if _data is None:
        return None
    def init_rule(m, *key):
        ## last key is time
        pm_t = key[-1]
        key = key[:-1]
        if len(key) == 0:
            return get_time_attr(_data, pm_t)
        if len(key) == 1:
            key = key[0]
        if key in _data:
            return get_time_attr(_data[key], pm_t)
        else:
            return None

    def get_time_attr(att, pm_t):
        if isinstance(att, dict):
            if 'data_type' in att and att['data_type'] == 'time_series':
                return att['values'][pm_t-1]
            else:
                raise Exception("Unexpected dictionary {}".format(att))
        else:
            return att

    return init_rule
