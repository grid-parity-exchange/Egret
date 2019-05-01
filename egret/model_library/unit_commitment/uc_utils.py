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

def add_model_attr(attr, requires = {}):
    def actual_decorator(func):
        @wraps(func)
        def wrapper(*args, **kwds):
            ## tag this function in the model with the appropriate attribute
            model = args[0]
            if hasattr(model, attr):
                raise Exception("Exception adding %s! Model already has %s %s! You may only add one type of %s!"%(func.__name__, attr, getattr(model,attr), attr)) 
            # this checks to see if the required components were already added
            for base_attr in requires:
                if not hasattr(model, base_attr):
                    raise Exception("Exception adding %s! %s requires some %s to be added first!"%(func.__name__, func.__name__, base_attr)) 
                ## None in this context means there is no specific requirement
                if requires[base_attr] is None:
                    continue
                if getattr(model, base_attr) not in requires[base_attr]:
                    raise Exception("Exception adding %s! %s requires one of: "%(func.__name__, func.__name__) + ", ".join(requires[base_attr]) + ", to be added first.")
            setattr(model, attr, func.__name__)
            return func(*args, **kwds)
        return wrapper
    return actual_decorator

## provides a view on grid_data attributes that
## is handy for building pyomo params
def build_uc_time_mapping(md_timeperiods):
    ## Assums the last key is time
    def uc_time_helper(_data):
        ## if there is no data,
        ## we return None to the initializer
        if _data is None:
            return None
        def init_rule(m, *key):
            if not isinstance(key, tuple):
                ## time indexed
                pm_t = key
                att = _data
                return self.get_time_attr(att, pm_t)
            else:
                ## last key is time
                pm_t = key[-1]
                key = key[:-1]
                if len(key) == 0:
                    att = _data
                elif len(key) == 1:
                    key = key[0]
                    att = _data[key]
                else:
                    att = _data[key]
                return get_time_attr(att, pm_t)

        def get_time_attr(att, pm_t):
            if isinstance(att, dict):
                if 'data_type' in att and att['data_type'] == 'time_series':
                    return att['values'][md_timeperiods[pm_t-1]]
                else:
                    raise Exception("Unexpected dictionary {}".format(att))
            else:
                return att

        return init_rule
    return uc_time_helper
