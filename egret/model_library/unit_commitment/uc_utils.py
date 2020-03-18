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
                if (not hasattr(model, base_attr)) or (getattr(model, base_attr) is None):
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
def uc_time_helper(model_time_periods):
    TimePeriods = list(model_time_periods)

    def dict_constructor(_data):
        return_dict = dict()
        ## if there is no data,
        ## we return dict() to the initializer
        if _data is None or _data == dict():
            return return_dict
        ## if the _data is a non-empty dictionary,
        ## then either this "thing" is itself indexed,
        ## or it is a time series for one thing
        if isinstance(_data,dict):
            ## in this case, this is a time series of one "thing"
            if 'data_type' in _data and _data['data_type'] == 'time_series':
                values = _data['values']
                for i,t in enumerate(TimePeriods):
                    return_dict[t] = values[i]
            else: ## it's a dictionary of things, which are potentially time indexed
                for key, att in _data.items():
                    if isinstance(att, dict):
                        if 'data_type' in att and att['data_type'] == 'time_series':
                            values = att['values']
                            for i,t in enumerate(TimePeriods):
                                return_dict[key,t] = values[i]
                        else: ## assume we know what to do with it, not copying
                            for t in TimePeriods:
                                return_dict[key,t] = att
                    else:
                        for t in TimePeriods:
                            return_dict[key,t] = att
        else:
            for t in TimePeriods:
                return_dict[t] = _data
        return return_dict

    return dict_constructor
