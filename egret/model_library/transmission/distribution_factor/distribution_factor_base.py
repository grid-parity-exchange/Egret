from pyomo.core.base.block import _BlockData, Block
import pyomo.environ as pe
from collections import Iterable
from pyomo.core.kernel.component_map import ComponentMap
from pyomo.core.kernel.component_set import ComponentSet
from egret.model_library.defn import DistributionFactorType, BasePointType
import sys
"""
Base class for distribution factor approximations
"""

@declare_custom_block(name='BaseDF')
class BaseDistributionFactorData(_BlockData):
    def __init__(self, component):
        _BlockData.__init__(self, component)
        self._persistent_solvers = ComponentSet()
        self._allow_changes = False

    def add_component(self, name, val):
        if self._allow_changes:
            _BlockData.add_component(self, name, val)
        else:
            raise RuntimeError('Pyomo components cannot be added to objects of type {0}.'.format(type(self)))

    def _set_input(self, distribution_factor_type=DistributionFactorType.PTDF, base_point_type=BasePointType.FLATSTART, persistent_solvers=None):
        self._persistent_solvers = persistent_solvers
        if self._persistent_solvers is None:
            self._persistent_solvers = ComponentSet()
        if not isinstance(self._persistent_solvers, Iterable):
            self._persistent_solvers = ComponentSet([self._persistent_solvers])
        else:
            self._persistent_solvers = ComponentSet(self._persistent_solvers)
        self._distribution_factor_type = distribution_factor_type
        assert self._distribution_factor_type in DistributionFactorType
        self._base_point_type = base_point_type
        assert self._base_point_type in BasePointType

    def remove(self):
        """
        Remove any auto-created vars/constraints from the distribution factor block
        """
        self._remove_from_persistent_solvers()
        comps = [pe.Block, pe.Constraint, pe.Var, pe.Set, pe.Param]
        for comp in comps:
            comps_to_del = list(self.component_objects([comp], descend_into=False))
            for _comp in comps_to_del:
                self.del_component(_comp)
        for comp in comps:
            comps_to_del = list(self.component_data_objects([comp], descend_into=False))
            for _comp in comps_to_del:
                self.del_component(_comp)

    def rebuild(self):
        """
        Remove any auto-created vars/constraints from the distribution factor block and recreate it
        """
        self._allow_changes = True
        self.remove()
        self._build()
        self._add_to_persistent_solvers()
        self._allow_changes = False

    def _build(self):
        """
        Build the auto-created vars/constraints that form the distribution factor constraint set
        """
        raise NotImplementedError('This should be implemented in the derived class.')

    def _remove_from_persistent_solvers(self):
        for i in self._persistent_solvers:
            i.remove_block(block=self)

    def _add_to_persistent_solvers(self):
        for i in self._persistent_solvers:
            i.add_block(block=self)

    def add_persistent_solver(self, persistent_solver):
        self._persistent_solvers.add(persistent_solver)

    def remove_persistent_solver(self, persistent_solver):
        self._persistent_solvers.remove(persistent_solver)

    def clear_persistent_solvers(self):
        self._persistent_solvers = ComponentSet()

    def get_abs_violation(self):
        """
        Compute the absolute value of the constraint violation given the current values of the corresponding vars.

        Returns
        -------
        float
        """
        raise NotImplementedError('This method should be implemented in the derived class.')


    def _get_pprint_string(self, relational_operator_string):
        raise NotImplementedError('This method should be implemented by subclasses.')

    def pprint(self, filename=None, ostream=None, verbose=False, prefix=""):
        if filename is not None:
            output = open(filename, 'w')
            self.pprint(ostream=output, verbose=verbose, prefix=prefix)
            output.close()
            return

        if ostream is None:
            ostream = sys.stdout

        ostream.write('{0}{1}\n'.format(prefix, self.name))

