from .acopf import _create_base_power_ac_model, create_rsv_acopf_model, create_psv_acopf_model, create_atan_acopf_model
import egret.model_library.transmission.branch as libbranch
from egret.data.data_utils import map_items, zip_items
from collections import OrderedDict
import coramin
from coramin.relaxations.custom_block import declare_custom_block
from coramin.relaxations.relaxations_base import BaseRelaxationData, ComponentWeakRef
from pyomo.core.expr.numvalue import NumericConstant, is_constant
import pyomo.environ as pe
import math
from math import sqrt
import logging
from pyomo.common.collections.orderedset import OrderedSet


logger = logging.getLogger(__name__)


def _relaxation_helper(model, md, include_soc, use_linear_relaxation, use_fbbt=True):
    coramin.relaxations.relax(model,
                              in_place=True,
                              use_fbbt=use_fbbt,
                              fbbt_options={'deactivate_satisfied_constraints': True,
                                            'max_iter': 2})
    if not use_linear_relaxation:
        for b in coramin.relaxations.relaxation_data_objects(model, descend_into=True, active=True, sort=True):
            if not isinstance(b, coramin.relaxations.PWMcCormickRelaxationData):
                b.use_linear_relaxation = False
                b.rebuild()

    if include_soc:
        branch_attrs = md.attributes(element_type='branch')
        bus_pairs = zip_items(branch_attrs['from_bus'], branch_attrs['to_bus'])
        unique_bus_pairs = list(OrderedDict((val, None) for idx, val in bus_pairs.items()).keys())
        libbranch.declare_ineq_soc(model=model, index_set=unique_bus_pairs,
                                   use_outer_approximation=use_linear_relaxation)


def _check_bounds(v, lb, ub):
    if lb is None:
        lb = -math.inf
    if ub is None:
        ub = math.inf
    if ub < lb:
        raise ValueError('Variable lb is greater than ub')
    if not v.is_fixed() and (ub - lb) <= 1e-8:
        logger.warning('Variable is not fixed, but (ub - lb) is small; var: {0}; lb: {1}; ub: {2}'.format(str(v), lb, ub))


def _get_bounds(v):
    if v.is_fixed():
        lb = ub = v.value
    else:
        lb = v.lb
        ub = v.ub
    _check_bounds(v, lb, ub)
    return lb, ub


@declare_custom_block(name='SOCEdgeCuts')
class SOCEdgeCutsData(BaseRelaxationData):
    """
    Relaxation for

      0 >= vmsq_1 * vmsq_2 - c**2 - s**2

    based on

    Kocuk, B., Dey, S.S. & Sun, X.A. Matrix minor reformulation
    and SOCP-based spatial branch-and-cut method for the AC
    optimal power flow problem. Math. Prog. Comp. 10, 557–596
    (2018). https://doi.org/10.1007/s12532-018-0150-9
    """
    def __init__(self, component):
        BaseRelaxationData.__init__(self, component)
        self._cref = ComponentWeakRef(None)
        self._sref = ComponentWeakRef(None)
        self._vmsq_1_ref = ComponentWeakRef(None)
        self._vmsq_2_ref = ComponentWeakRef(None)
        self._aux_var = NumericConstant(0)

    @property
    def _c(self):
        return self._cref.get_component()

    @property
    def _s(self):
        return self._sref.get_component()

    @property
    def _vmsq_1(self):
        return self._vmsq_1_ref.get_component()

    @property
    def _vmsq_2(self):
        return self._vmsq_2_ref.get_component()

    def get_rhs_vars(self):
        return [self._c, self._s, self._vmsq_1, self._vmsq_2]

    def get_rhs_expr(self):
        return self._vmsq_1 * self._vmsq_2 - self._c**2 - self._s**2

    def vars_with_bounds_in_relaxation(self):
        return self.get_rhs_vars()

    def set_input(self, c, s, vmsq_1, vmsq_2, persistent_solvers=None):
        self._set_input(relaxation_side=coramin.utils.RelaxationSide.UNDER,
                        persistent_solvers=persistent_solvers,
                        use_linear_relaxation=True)
        self._cref.set_component(c)
        self._sref.set_component(s)
        self._vmsq_1_ref.set_component(vmsq_1)
        self._vmsq_2_ref.set_component(vmsq_2)

    def build(self, c, s, vmsq_1, vmsq_2, persistent_solvers=None):
        self.set_input(c=c, s=s, vmsq_1=vmsq_1, vmsq_2=vmsq_2, persistent_solvers=persistent_solvers)
        self.rebuild()

    def _build_relaxation(self):
        clb, cub = _get_bounds(self._c)
        slb, sub = _get_bounds(self._s)
        vmsq_1_lb, vmsq_1_ub = _get_bounds(self._vmsq_1)
        vmsq_2_lb, vmsq_2_ub = _get_bounds(self._vmsq_2)

        if None in {clb, cub, slb, sub, vmsq_1_lb, vmsq_1_ub, vmsq_2_lb, vmsq_2_ub}:
            return None

        if vmsq_1_lb == vmsq_1_ub and vmsq_2_lb == vmsq_2_ub:
            if clb == cub or slb == sub:
                rhs = [vmsq_1_lb * vmsq_2_lb]
            else:
                rhs = [math.sqrt(vmsq_1_lb * vmsq_2_lb)]
        elif vmsq_1_lb == vmsq_1_ub:
            if clb == cub or slb == sub:
                rhs = [vmsq_1_lb * self._vmsq_2]
            else:
                m = (math.sqrt(vmsq_2_ub) - math.sqrt(vmsq_2_lb)) / (vmsq_2_ub - vmsq_2_lb)
                b = math.sqrt(vmsq_2_ub) - m * vmsq_2_ub
                rhs = [math.sqrt(vmsq_1_lb) * (m * self._vmsq_2 + b)]
        elif vmsq_2_lb == vmsq_2_ub:
            if clb == cub or slb == sub:
                rhs = [vmsq_2_lb * self._vmsq_1]
            else:
                m = (math.sqrt(vmsq_1_ub) - math.sqrt(vmsq_1_lb)) / (vmsq_1_ub - vmsq_1_lb)
                b = math.sqrt(vmsq_1_ub) - m * vmsq_1_ub
                rhs = [math.sqrt(vmsq_2_lb) * (m * self._vmsq_1 + b)]
        else:
            if clb == cub or slb == sub:
                rhs = [vmsq_1_lb * self._vmsq_2 + self._vmsq_1 * vmsq_2_lb - vmsq_1_lb * vmsq_2_lb,
                       vmsq_1_ub * self._vmsq_2 + self._vmsq_1 * vmsq_2_ub - vmsq_1_ub * vmsq_2_ub]
            else:
                cm = [sqrt(vmsq_1_lb) / (sqrt(vmsq_2_lb) + sqrt(vmsq_2_ub)),
                      sqrt(vmsq_1_ub) / (sqrt(vmsq_2_lb) + sqrt(vmsq_2_ub))]
                bm = [sqrt(vmsq_2_lb) / (sqrt(vmsq_1_lb) + sqrt(vmsq_1_ub)),
                      sqrt(vmsq_2_ub) / (sqrt(vmsq_1_lb) + sqrt(vmsq_1_ub))]
                am = [sqrt(vmsq_1_lb * vmsq_2_lb) - bm[0] * vmsq_1_lb - cm[0] * vmsq_2_lb,
                      sqrt(vmsq_1_ub * vmsq_2_ub) - bm[1] * vmsq_1_ub - cm[1] * vmsq_2_ub]
                rhs = list()
                for i in [0, 1]:
                    rhs.append(am[i] + bm[i] * self._vmsq_1 + cm[i] * self._vmsq_2)

        if clb == cub and slb == sub:
            lhs = [clb**2 + slb**2]
        elif clb == cub:
            m = (sub**2 - slb**2) / (sub - slb)
            b = sub**2 - m * sub
            lhs = [clb**2 + m * self._s + b]
        elif slb == sub:
            m = (cub**2 - clb**2) / (cub - clb)
            b = cub**2 - m * cub
            lhs = [m * self._c + b + slb**2]
        else:
            if sqrt(cub ** 2 + sub ** 2) + sqrt(clb ** 2 + slb ** 2) - sqrt(cub ** 2 + slb ** 2) - sqrt(
                    clb ** 2 + sub ** 2) <= 0:
                cn = [(sqrt(clb ** 2 + sub ** 2) - sqrt(clb ** 2 + slb ** 2)) / (sub - slb),
                      (sqrt(cub ** 2 + sub ** 2) - sqrt(cub ** 2 + slb ** 2)) / (sub - slb)]
                bn = [(sqrt(cub ** 2 + slb ** 2) - sqrt(clb ** 2 + slb ** 2)) / (cub - clb),
                      (sqrt(cub ** 2 + sub ** 2) - sqrt(clb ** 2 + sub ** 2)) / (cub - clb)]
                an = [sqrt(clb ** 2 + slb ** 2) - bn[0] * clb - cn[0] * slb,
                      sqrt(cub ** 2 + sub ** 2) - bn[1] * cub - cn[1] * sub]
            else:
                cn = [(sqrt(clb ** 2 + sub ** 2) - sqrt(clb ** 2 + slb ** 2)) / (sub - slb),
                      (sqrt(cub ** 2 + sub ** 2) - sqrt(cub ** 2 + slb ** 2)) / (sub - slb)]
                bn = [(sqrt(cub ** 2 + sub ** 2) - sqrt(clb ** 2 + sub ** 2)) / (cub - clb),
                      (sqrt(cub ** 2 + slb ** 2) - sqrt(clb ** 2 + slb ** 2)) / (cub - clb)]
                an = [sqrt(clb ** 2 + slb ** 2) - bn[0] * clb - cn[0] * slb,
                      sqrt(cub ** 2 + sub ** 2) - bn[1] * cub - cn[1] * sub]
            lhs = list()
            for i in [0, 1]:
                lhs.append(an[i] + bn[i] * self._c + cn[i] * self._s)

        self.underestimators = pe.ConstraintList()
        for _lhs in lhs:
            for _rhs in rhs:
                if is_constant(_lhs) and is_constant(_rhs):
                    continue
                self.underestimators.add(_lhs >= _rhs)

    def is_rhs_convex(self):
        return False

    def is_rhs_concave(self):
        return False

    @property
    def use_linear_relaxation(self):
        return self._use_linear_relaxation

    @use_linear_relaxation.setter
    def use_linear_relaxation(self, val):
        if val is not True:
            raise ValueError('The SOCEdgeCuts class can only produce linear relaxations')
        self._use_linear_relaxation = True


def create_soc_relaxation(model_data,
                          use_linear_relaxation=True,
                          include_feasibility_slack=False,
                          use_fbbt=True,
                          keep_vars_for_out_of_service_elements=False):
    model, md = _create_base_power_ac_model(model_data, include_feasibility_slack=include_feasibility_slack,
                                            keep_vars_for_out_of_service_elements=keep_vars_for_out_of_service_elements)
    if use_linear_relaxation:
        _relaxation_helper(model=model,
                           md=md,
                           include_soc=True,
                           use_linear_relaxation=use_linear_relaxation,
                           use_fbbt=use_fbbt)
    else:
        branch_attrs = md.attributes(element_type='branch')
        bus_pairs = zip_items(branch_attrs['from_bus'], branch_attrs['to_bus'])
        unique_bus_pairs = list(OrderedDict((val, None) for idx, val in bus_pairs.items()).keys())
        libbranch.declare_ineq_soc(model=model, index_set=unique_bus_pairs,
                                   use_outer_approximation=use_linear_relaxation)
    return model, md


def create_atan_relaxation(model_data,
                           use_linear_relaxation=True,
                           include_feasibility_slack=False,
                           use_soc_edge_cuts=False,
                           use_fbbt=True,
                           keep_vars_for_out_of_service_elements=False):
    model, md = create_atan_acopf_model(model_data=model_data, include_feasibility_slack=include_feasibility_slack,
                                        keep_vars_for_out_of_service_elements=keep_vars_for_out_of_service_elements)
    del model.ineq_soc
    del model._con_ineq_soc
    if use_soc_edge_cuts:
        del model.ineq_soc_ub
        del model._con_ineq_soc_ub
    _relaxation_helper(model=model,
                       md=md,
                       include_soc=True,
                       use_linear_relaxation=use_linear_relaxation,
                       use_fbbt=use_fbbt)
    if use_soc_edge_cuts:
        branch_attrs = md.attributes(element_type='branch')
        bus_pairs = zip_items(branch_attrs['from_bus'], branch_attrs['to_bus'])
        unique_bus_pairs = OrderedSet(val for val in bus_pairs.values())
        model.ineq_soc_ub_set = pe.Set(initialize=list(unique_bus_pairs))
        model.ineq_soc_ub = SOCEdgeCuts(model.ineq_soc_ub_set)
        for from_bus, to_bus in model.ineq_soc_ub_set:
            model.ineq_soc_ub[from_bus, to_bus].build(c=model.c[from_bus, to_bus],
                                                      s=model.s[from_bus, to_bus],
                                                      vmsq_1=model.vmsq[from_bus],
                                                      vmsq_2=model.vmsq[to_bus])
    return model, md


def create_polar_acopf_relaxation(model_data,
                                  include_soc=True,
                                  use_linear_relaxation=True,
                                  include_feasibility_slack=False,
                                  use_fbbt=True,
                                  keep_vars_for_out_of_service_elements=False):
    model, md = create_psv_acopf_model(model_data, include_feasibility_slack=include_feasibility_slack,
                                       keep_vars_for_out_of_service_elements=keep_vars_for_out_of_service_elements)
    _relaxation_helper(model=model,
                       md=md,
                       include_soc=include_soc,
                       use_linear_relaxation=use_linear_relaxation,
                       use_fbbt=use_fbbt)
    return model, md


def create_rectangular_acopf_relaxation(model_data,
                                        include_soc=True,
                                        use_linear_relaxation=True,
                                        include_feasibility_slack=False,
                                        use_fbbt=True,
                                        keep_vars_for_out_of_service_elements=False):
    model, md = create_rsv_acopf_model(model_data, include_feasibility_slack=include_feasibility_slack,
                                       keep_vars_for_out_of_service_elements=keep_vars_for_out_of_service_elements)
    _relaxation_helper(model=model,
                       md=md,
                       include_soc=include_soc,
                       use_linear_relaxation=use_linear_relaxation,
                       use_fbbt=use_fbbt)
    return model, md
