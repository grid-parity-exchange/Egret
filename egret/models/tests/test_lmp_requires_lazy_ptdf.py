#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

'''
LMP extraction from a PTDF network model requires the lazy formulation.

`ptdf_utils.VirtualPTDFMatrix.calculate_LMP` is written against the lazy
constraint generation model shape.  It reads

  * `self.branches_keys_masked` / `self.B_dA_masked` and the other "masked"
    attributes, which are built by `_get_filtered_lines()` -- called only when
    `ptdf_options['lazy']` is True, and
  * `mb.ineq_pf_branch_thermal_bounds` / `mb.ineq_pf_interface_bounds`, which
    are declared only on the lazy path (the non-lazy path declares
    `ineq_pf_branch_thermal_lb` / `_ub` instead).

Both callers invoke it without checking `lazy`:

  * `dcopf.solve_dcopf` for `create_ptdf_dcopf_model`, and
  * `unit_commitment._save_uc_results` when `relaxed=True`.

So `lazy=False` used to fail with an unhelpful
`AttributeError: 'VirtualPTDFMatrix' object has no attribute
'branches_keys_masked'`.  These tests pin the clear error, and pin that the
lazy path still works.

They also cover the two solver-capability calls on the lazy path, which used to
be made unconditionally and made GLPK/CBC/HiGHS unusable:

  * `warmstart=True` was always passed to `solver.solve()`.  Pyomo passes
    unrecognized keywords through to the file writer, which rejects the key
    itself, so the failure was `ValueError: key 'warmstart' not defined for
    ConfigDict '' and implicit (undefined) keys are not allowed` raised from
    `pyomo/repn/plugins/lp_writer.py`.
  * `solver.reset()` was called even though the adjacent comment says "If the
    solver implements it, reset the solver".
'''
import os

import pytest

from egret.data.model_data import ModelData
from egret.models import unit_commitment as uc
from egret.models.dcopf import solve_dcopf, create_ptdf_dcopf_model
from egret.models.unit_commitment import solve_unit_commitment

current_dir = os.path.dirname(os.path.abspath(__file__))

DCOPF_CASE = os.path.join(current_dir, 'transmission_test_instances',
                          'test_instances', 'hvdc_test_case3.json')
UC_CASE = os.path.join(current_dir, 'uc_test_instances', 'tiny_uc_tc.json')


def _slack_type():
    return uc.SlackType.TRANSMISSION_LIMITS


def _available_solver():
    '''
    First available solver, following the convention used by
    test_unit_commitment.py so these tests run wherever the suite runs.
    '''
    from pyomo.opt import SolverFactory

    for name in ('xpress_persistent', 'gurobi_persistent', 'cplex_persistent',
                 'gurobi', 'cplex', 'cbc', 'glpk', 'highs'):
        try:
            if SolverFactory(name).available(exception_flag=False):
                return name
        except Exception:
            continue
    return None


def _linearize_costs(md):
    '''
    Drop the >= 2nd order terms from polynomial cost curves.

    hvdc_test_case3 has quadratic curves, which is why the upstream DCOPF tests
    run it with ipopt.  Removing the higher order terms lets LP/MILP solvers
    (GLPK, HiGHS, CPLEX) solve it, so these tests do not need a nonlinear
    solver.  Note that a zero coefficient still produces a nonlinear
    expression, so the term has to be removed rather than zeroed.
    '''
    for _, gen in md.elements(element_type='generator'):
        p_cost = gen.get('p_cost')
        if p_cost and p_cost.get('cost_curve_type') == 'polynomial':
            values = dict(p_cost['values'])
            if any(int(order) >= 2 for order in values):
                p_cost['values'] = {order: coef for order, coef in values.items()
                                    if int(order) < 2}


# ---------------------------------------------------------------------------
# The guards
# ---------------------------------------------------------------------------

def test_dcopf_ptdf_lazy_false_raises_a_clear_error():
    '''
    create_ptdf_dcopf_model accepts ptdf_options={'lazy': False}, but
    solve_dcopf cannot extract LMPs from that model.  It must say so instead of
    raising AttributeError.  The guard fires right after the model is built,
    so no solver is needed here.
    '''
    md = ModelData(DCOPF_CASE)
    with pytest.raises(RuntimeError, match='requires lazy PTDF'):
        solve_dcopf(md, 'glpk',
                    dcopf_model_generator=create_ptdf_dcopf_model,
                    solver_tee=False,
                    ptdf_options={'lazy': False})


def test_dcopf_ptdf_lazy_true_still_produces_lmps():
    '''
    Control: the default (lazy=True) path must keep working and still produce
    LMPs -- the guard must not block the supported configuration.
    '''
    solver = _available_solver()
    if solver is None:
        pytest.skip('no LP/MIP solver available')

    md = ModelData(DCOPF_CASE)
    _linearize_costs(md)
    md = solve_dcopf(md, solver,
                     dcopf_model_generator=create_ptdf_dcopf_model,
                     solver_tee=False,
                     ptdf_options={'lazy': True})

    buses = md.data['elements']['bus']
    assert len(buses) > 0
    for bus_name, bus in buses.items():
        assert bus.get('lmp') is not None, \
            'bus {} has no lmp on the supported lazy path'.format(bus_name)


def test_uc_relaxed_lazy_false_raises_a_clear_error():
    '''
    solve_unit_commitment(relaxed=True) asks for LMPs, which needs lazy=True.
    '''
    solver = _available_solver()
    if solver is None:
        pytest.skip('no LP/MIP solver available')

    md = ModelData(UC_CASE)
    with pytest.raises(RuntimeError, match='requires lazy PTDF'):
        solve_unit_commitment(md, solver, relaxed=True,
                              slack_type=_slack_type(),
                              network_constraints='ptdf_power_flow',
                              ptdf_options={'lazy': False})


def test_uc_relaxed_lazy_true_still_produces_lmps():
    '''
    Control: the guard must not block the supported configuration.
    '''
    solver = _available_solver()
    if solver is None:
        pytest.skip('no LP/MIP solver available')

    md = ModelData(UC_CASE)
    md = solve_unit_commitment(md, solver, relaxed=True,
                               slack_type=_slack_type(),
                               network_constraints='ptdf_power_flow')

    buses = md.data['elements']['bus']
    assert len(buses) > 0
    for bus_name, bus in buses.items():
        assert bus.get('lmp') is not None, \
            'bus {} has no lmp on the supported lazy path'.format(bus_name)


# ---------------------------------------------------------------------------
# Solver capability handling on the lazy path
# ---------------------------------------------------------------------------

class _RecordingSolver(object):
    '''Minimal stand-in that records the kwargs it was called with.'''

    def __init__(self, warm_start_capable, with_reset=False):
        self._warm_start_capable = warm_start_capable
        self.seen_kwargs = None
        self.reset_calls = 0
        if with_reset:
            ## Only attach the attribute when asked, so getattr() differs
            self.reset = self._reset

    def _reset(self):
        self.reset_calls += 1

    def warm_start_capable(self):
        return self._warm_start_capable

    def solve(self, model, **kwargs):
        self.seen_kwargs = kwargs
        return 'results'


class _FakeModel(object):
    '''Just enough of a model for _lazy_ptdf_solve.'''

    class solutions(object):
        @staticmethod
        def load_from(results):
            pass

    def model(self):
        return self


def test_warm_start_capable_detection():
    assert uc._solver_warm_start_capable(_RecordingSolver(True)) is True
    assert uc._solver_warm_start_capable(_RecordingSolver(False)) is False

    class NoCapabilityMethod(object):
        pass

    ## A solver object without the method must be treated as incapable rather
    ## than raising.
    assert uc._solver_warm_start_capable(NoCapabilityMethod()) is False

    class ExplodingCapability(object):
        def warm_start_capable(self):
            raise RuntimeError('boom')

    assert uc._solver_warm_start_capable(ExplodingCapability()) is False


@pytest.mark.parametrize('persistent_solver', [True, False])
def test_warmstart_keyword_is_omitted_for_incapable_solver(persistent_solver):
    '''
    Regression: the 'warmstart' keyword must not be sent at all when the solver
    does not support it -- not even as False, because Pyomo forwards unknown
    keywords to the file writer, which rejects the key itself.
    '''
    solver = _RecordingSolver(warm_start_capable=False)
    if persistent_solver:
        solver.load_vars = lambda vars_to_load=None: None

    uc._lazy_ptdf_solve(_FakeModel(), solver, persistent_solver,
                        symbolic_solver_labels=False, solver_tee=False,
                        vars_to_load=None)

    assert 'warmstart' not in solver.seen_kwargs


@pytest.mark.parametrize('persistent_solver', [True, False])
def test_warmstart_keyword_is_passed_for_capable_solver(persistent_solver):
    '''CPLEX/Gurobi keep their warm start -- this is not a behaviour removal.'''
    solver = _RecordingSolver(warm_start_capable=True)
    if persistent_solver:
        solver.load_vars = lambda vars_to_load=None: None

    uc._lazy_ptdf_solve(_FakeModel(), solver, persistent_solver,
                        symbolic_solver_labels=False, solver_tee=False,
                        vars_to_load=None)

    assert solver.seen_kwargs.get('warmstart') is True


def test_reset_is_optional_on_the_solver(monkeypatch):
    '''
    Regression: solver.reset() used to be called unconditionally, even though
    the comment next to it says "If the solver implements it".  Solvers built
    on the newer pyomo.contrib.solver API (e.g. HiGHS) do not have reset().
    '''
    monkeypatch.setattr(uc, '_lazy_ptdf_warmstart_copy_violations',
                        lambda *args, **kwargs: None)
    monkeypatch.setattr(uc, '_lazy_ptdf_solve',
                        lambda *args, **kwargs: 'results')

    def cleanup(solver):
        return uc._lazy_ptdf_termination_cleanup(
            m=None, md=None, time_periods=[], solver=solver, ptdf_options={},
            prepend_str='', t_subset=[], persistent_solver=False,
            symbolic_solver_labels=False, solver_tee=False, vars_to_load=None,
            solve_method_options=None, duals=False, warmstart_loop=True,
            results=None)

    ## A solver without reset() must be tolerated ...
    assert cleanup(_RecordingSolver(True, with_reset=False)) == 'results'

    ## ... while a solver that has it must still be reset, exactly as before.
    with_reset = _RecordingSolver(True, with_reset=True)
    assert cleanup(with_reset) == 'results'
    assert with_reset.reset_calls == 1
