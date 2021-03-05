from .acopf import _create_base_power_ac_model, create_rsv_acopf_model, create_psv_acopf_model, create_atan_acopf_model
import egret.model_library.transmission.branch as libbranch
from egret.data.data_utils import map_items, zip_items
from collections import OrderedDict
try:
    import coramin
    coramin_available = True
except ImportError:
    coramin_available = False


def _relaxation_helper(model, md, include_soc, use_linear_relaxation):
    if not coramin_available:
        raise ImportError('Cannot create relaxation unless coramin is available.')
    coramin.relaxations.relax(model,
                              in_place=True,
                              use_fbbt=True,
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


def create_soc_relaxation(model_data, use_linear_relaxation=True, include_feasibility_slack=False):
    model, md = _create_base_power_ac_model(model_data, include_feasibility_slack=include_feasibility_slack)
    if use_linear_relaxation:
        _relaxation_helper(model=model, md=md, include_soc=True, use_linear_relaxation=use_linear_relaxation)
    else:
        branch_attrs = md.attributes(element_type='branch')
        bus_pairs = zip_items(branch_attrs['from_bus'], branch_attrs['to_bus'])
        unique_bus_pairs = list(OrderedDict((val, None) for idx, val in bus_pairs.items()).keys())
        libbranch.declare_ineq_soc(model=model, index_set=unique_bus_pairs,
                                   use_outer_approximation=use_linear_relaxation)
    return model, md


def create_atan_relaxation(model_data, use_linear_relaxation=True, include_feasibility_slack=False):
    model, md = create_atan_acopf_model(model_data=model_data, include_feasibility_slack=include_feasibility_slack)
    del model.ineq_soc
    del model._con_ineq_soc
    _relaxation_helper(model=model, md=md, include_soc=True, use_linear_relaxation=use_linear_relaxation)
    return model, md


def create_polar_acopf_relaxation(model_data, include_soc=True, use_linear_relaxation=True, include_feasibility_slack=False):
    model, md = create_psv_acopf_model(model_data, include_feasibility_slack=include_feasibility_slack)
    _relaxation_helper(model=model, md=md, include_soc=include_soc, use_linear_relaxation=use_linear_relaxation)
    return model, md


def create_rectangular_acopf_relaxation(model_data, include_soc=True, use_linear_relaxation=True, include_feasibility_slack=False):
    model, md = create_rsv_acopf_model(model_data, include_feasibility_slack=include_feasibility_slack)
    _relaxation_helper(model=model, md=md, include_soc=include_soc, use_linear_relaxation=use_linear_relaxation)
    return model, md
