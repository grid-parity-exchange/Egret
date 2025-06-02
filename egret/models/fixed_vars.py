import pyomo.environ as pe
import weakref


"""
The functions below are motivated by the fact that pyomo now 
raises exceptions when a variable is fixed to a value outside. 
Suppose a generator with a strictly positive lower bound on 
power output is not in service. The power output needs to be fixed 
to zero, but that is outside the bounds. The functions below 
allow one to automatically remove and cache the variable bounds 
when fixing a variable and later restore those bounds when 
unfixing a variable.

Doing this is a bit tricky because variables are not hashable 
and we don't want a memory leak due to the cache. If variables 
were hashable, I would just use a WeakKeyDictionary. Instead,
I'm using the id of the variable as the key in the cache. As long 
as the object exists, its id is unique and constant. Therefore, 
it is a reliable key as long as it gets removed from the cache 
when the variable is garbage collected. We ensure this using 
weakref.finalize.
"""


_bounds_cache = {}


def _remove_var_from_cache(vid):
    # print(f'removing {vid}')
    del _bounds_cache[vid]


def fix_var_and_remove_bounds(v, val):
    # print(f'adding {id(v)}')
    if id(v) not in _bounds_cache:
        weakref.finalize(v, _remove_var_from_cache, vid=id(v))

    _bounds_cache[id(v)] = (v.domain, v._lb, v._ub)

    v.domain = pe.Reals
    v.setlb(None)
    v.setub(None)
    v.set_value(val, skip_validation=True)
    v.fix()


def unfix_var_and_restore_bounds(v):
    v.unfix()
    if id(v) in _bounds_cache:
        domain, lb, ub = _bounds_cache[id(v)]
        v.domain = domain
        v.setlb(lb)
        v.setub(ub)
