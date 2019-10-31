#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module contains several helper functions and classes that are useful when
modifying the data dictionary
"""
import abc
import pickle
import numpy as np
import egret.model_library.transmission.tx_calc as tx_calc

from egret.model_library.defn import BasePointType, ApproximationType
from egret.common.log import logger
from math import radians

def get_ptdf_potentially_from_file(ptdf_options, branches_keys,
                                   buses_keys, interfaces=None):
    '''
    small loop to get a PTDF matrix previously pickled, 
    returns None if not found
    '''

    PTDF = None
    PTDF_pickle = None
    if ptdf_options['load_from'] is not None:
        try:
            PTDF_pickle = pickle.load(open(ptdf_options['load_from'], 'rb'))
        except:
            print("Error loading PTDF matrix from pickle file, calculating from start")

    if PTDF_pickle is not None:
        ## This may be a dict of data_utils.PTDFMatrix objects or just an object
        if isinstance(PTDF_pickle, dict):
            for key, PTDFo in PTDF_pickle.items():
                if _is_consistent_ptdfm(PTDFo, branches_keys, buses_keys):
                    PTDF = PTDFo
                    break
        ## could be a single ptdf dict
        else:
            if _is_consistent_ptdfm(PTDF_pickle, branches_keys, buses_keys):
                PTDF = PTDF_pickle

    if PTDF is not None:
        PTDF._set_lazy_limits(ptdf_options)

        ## since it's simple, recalculate the
        ## interfaces from scratch each time
        if interfaces is None:
            interfaces = dict()
        PTDF._calculate_ptdf_interface(interfaces)

    return PTDF

def write_ptdf_potentially_to_file(ptdf_options, PTDF):
    if ptdf_options['save_to'] is not None:
        pickle.dump(PTDF, open(ptdf_options['save_to'], 'wb'))

def _is_consistent_ptdfm(ptdf_mat, branches_keys, buses_keys):
    '''
    Checks the branches and buses keys for agreement when loading
    PTDF matrix from disk.

    Parameters
    ----------
    ptdf_mat : PTDFMatrix or PTDFLossesMatrix
    branches_keys : iterable of branches
    buses_keys : iterable of buses


    Returns
    ------
    bool : True if the branches_keys and buses_keys are consistent
           with those in the object ptdf_mat

    '''
    return ( set(branches_keys) == set(ptdf_mat.branches_keys) and \
             set(buses_keys) == set(ptdf_mat.buses_keys) )

class PTDFMatrix(object):
    '''
    This is a helper 
    '''
    def __init__(self, branches, buses, reference_bus, base_point,
                        ptdf_options, branches_keys = None, buses_keys = None,
                        interfaces = None):
        '''
        Creates a new _PTDFMaxtrixManager object to provide
        some useful methods for interfacing with Egret pyomo models

        Parameters
        ----------
        '''
        self._branches = branches
        self._buses = buses
        self._reference_bus = reference_bus
        self._base_point = base_point
        if branches_keys is None:
            self.branches_keys = tuple(branches.keys())
        else:
            self.branches_keys = tuple(branches_keys)
        if buses_keys is None:
            self.buses_keys = tuple(buses.keys())
        else:
            self.buses_keys = tuple(buses_keys)
        self._branchname_to_index_map = {branch_n : i for i, branch_n in enumerate(self.branches_keys)}
        self._busname_to_index_map = {bus_n : j for j, bus_n in enumerate(self.buses_keys)}

        self.branch_limits_array = np.fromiter((branches[branch]['rating_long_term'] for branch in self.branches_keys), float, count=len(self.branches_keys))
        self.branch_limits_array.flags.writeable = False

        self._base_point = base_point
        self._calculate()

        self._set_lazy_limits(ptdf_options)

        if interfaces is None:
            interfaces = dict()
        self._calculate_ptdf_interface(interfaces)

    def _calculate(self):
        self._calculate_ptdf()
        self._calculate_phi_adjust()
        self._calculate_phase_shift()

    def _calculate_ptdf(self):
        '''
        do the PTDF calculation
        '''
        ## calculate and store the PTDF matrix
        PTDFM = tx_calc.calculate_ptdf(self._branches,self._buses,self.branches_keys,self.buses_keys,self._reference_bus,self._base_point,
                                        mapping_bus_to_idx=self._busname_to_index_map)

        self.PTDFM = PTDFM

        ## protect the array using numpy
        self.PTDFM.flags.writeable = False

    def _calculate_ptdf_interface(self, interfaces):
        self.interface_keys = tuple(interfaces.keys())

        self.PTDFM_I, self.PTDFM_I_const \
                = tx_calc.calculate_interface_sensitivities(interfaces,
                                            self.interface_keys,
                                            self.PTDFM,
                                            self.phase_shift_array,
                                            self.phi_adjust_array,
                                            self._branchname_to_index_map)

        ## protect the array using numpy
        self.PTDFM_I.flags.writeable = False
        self.PTDFM_I_const.flags.writeable = False

        self._interfacename_to_index_map = \
                { i_n: idx for idx, i_n in enumerate(self.interface_keys) }

    def _calculate_phi_from_phi_to(self):
        return tx_calc.calculate_phi_constant(self._branches,self.branches_keys,self.buses_keys,ApproximationType.PTDF, mapping_bus_to_idx=self._busname_to_index_map)

    def _calculate_phi_adjust(self):
        phi_from, phi_to = self._calculate_phi_from_phi_to()
        
        ## hold onto these for line outages
        self._phi_from = phi_from
        self._phi_to = phi_to

        phi_adjust_array = phi_from-phi_to

        ## sum across the rows to get the total impact, and convert
        ## to dense for fast operations later
        self.phi_adjust_array = phi_adjust_array.sum(axis=1).T.A[0]

        ## protect the array using numpy
        self.phi_adjust_array.flags.writeable = False

    def _calculate_phase_shift(self):
        
        phase_shift_array = np.fromiter(( -(1/branch['reactance']) * (radians(branch['transformer_phase_shift'])/branch['transformer_tap_ratio']) if (branch['branch_type'] == 'transformer') else 0. for branch in (self._branches[bn] for bn in self.branches_keys)), float, count=len(self.branches_keys))

        self.phase_shift_array = phase_shift_array
        ## protect the array using numpy
        self.phase_shift_array.flags.writeable = False


    def _get_filtered_lines(self, ptdf_options):
        if ptdf_options['branch_kv_threshold'] is None:
            ## Nothing to do
            self.branch_mask = np.arange(len(self.branch_limits_array))
            self.branches_keys_masked = self.branches_keys
            self.branchname_to_index_masked_map = self._branchname_to_index_map
            self.PTDFM_masked = self.PTDFM
            self.phase_shift_array_masked = self.phase_shift_array
            self.branch_limits_array_masked = self.branch_limits_array
            return

        branches = self._branches
        buses = self._buses
        branch_mask = list()
        one = (ptdf_options['kv_threshold_type'] == 'one')
        kv_limit = ptdf_options['branch_kv_threshold']

        for i, bn in enumerate(self.branches_keys):
            branch = branches[bn]
            fb = buses[branch['from_bus']]

            fbt = True
            ## NOTE: this warning will be printed only once if we just check the from_bus
            if 'base_kv' not in fb:
                logger.warning("WARNING: did not find 'base_kv' for bus {}, considering it large for the purposes of filtering".format(branch['from_bus']))
            elif fb['base_kv'] < kv_limit:
                fbt = False

            if fbt and one:
                branch_mask.append(i)
                continue

            tb = buses[branch['to_bus']]
            tbt = False
            if ('base_kv' not in tb) or tb['base_kv'] >= kv_limit:
                tbt = True

            if fbt and tbt:
                branch_mask.append(i)
            elif one and tbt:
                branch_mask.append(i)

        self.branch_mask = np.array(branch_mask)
        self.branches_keys_masked = tuple(self.branches_keys[i] for i in self.branch_mask)
        self.branchname_to_index_masked_map = { bn : i for i,bn in enumerate(self.branches_keys_masked) }
        self.PTDFM_masked = self.PTDFM[branch_mask]
        self.phase_shift_array_masked = self.phase_shift_array[branch_mask]
        self.branch_limits_array_masked = self.branch_limits_array[branch_mask]

    def _set_lazy_limits(self, ptdf_options):
        if ptdf_options['lazy']:
            self._get_filtered_lines(ptdf_options)
            ## add / reset the relative limits based on the current options
            branch_limits = self.branch_limits_array_masked
            rel_flow_tol = ptdf_options['rel_flow_tol']
            abs_flow_tol = ptdf_options['abs_flow_tol']
            lazy_flow_tol = ptdf_options['lazy_rel_flow_tol']

            ## only enforce the relative and absolute, within tollerance
            self.enforced_branch_limits = np.maximum(branch_limits*(1+rel_flow_tol), branch_limits+abs_flow_tol)
            ## make sure the lazy limits are a superset of the enforce limits
            self.lazy_branch_limits = np.minimum(branch_limits*(1+lazy_flow_tol), self.enforced_branch_limits)


    def get_branch_ptdf_iterator(self, branch_name):
        row_idx = self._branchname_to_index_map[branch_name]
        ## get the row slice
        PTDF_row = self.PTDFM[row_idx]
        yield from zip(self.buses_keys, PTDF_row)

    def get_branch_ptdf_abs_max(self, branch_name):
        row_idx = self._branchname_to_index_map[branch_name]
        ## get the row slice
        PTDF_row = self.PTDFM[row_idx]
        return np.abs(PTDF_row).max()

    def get_branch_phase_shift(self, branch_name):
        return self.phase_shift_array[self._branchname_to_index_map[branch_name]]

    def get_bus_phi_adj(self, bus_name):
        return self.phi_adjust_array[self._busname_to_index_map[bus_name]]

    def get_branch_phi_adj(self, branch_name):
        row_idx = self._branchname_to_index_map[branch_name]
        ## get the row slice
        PTDF_row = self.PTDFM[row_idx]
        return PTDF_row.dot(self.phi_adjust_array)

    def bus_iterator(self):
        yield from self.buses_keys

    def get_interface_const(self, interface_name):
        return self.PTDFM_I_const[self._interfacename_to_index_map[interface_name]]

    def get_interface_ptdf_abs_max(self, interface_name):
        row_idx = self._interfacename_to_index_map[interface_name]
        ## get the row slice
        PTDF_I_row = self.PTDFM_I[row_idx]
        return np.abs(PTDF_I_row).max()

    def get_interface_ptdf_iterator(self, interface_name):
        row_idx = self._interfacename_to_index_map[interface_name]
        ## get the row slice
        PTDF_I_row = self.PTDFM_I[row_idx]
        yield from zip(self.buses_keys, PTDF_I_row)


class PTDFLossesMatrix(PTDFMatrix):

    def _calculate(self):
        self._calculate_ptdf()
        self._calculate_phi_adjust()
        self._calculate_phi_loss_constant()
        self._calculate_phase_shift()
        self._calculate_losses_phase_shift()

    def _calculate_ptdf(self):
        ptdf_r, ldf, ldf_c = tx_calc.calculate_ptdf_ldf(self._branches,self._buses,self.branches_keys,self.buses_keys,self._reference_bus,self._base_point,\
                                                        mapping_bus_to_idx=self._busname_to_index_map)

        self.PTDFM = ptdf_r
        self.LDF = ldf
        self.LDF_C = ldf_c

        ## protect the arrays using numpy
        self.PTDFM.flags.writeable = False
        self.LDF.flags.writeable = False
        self.LDF_C.flags.writeable = False

    def _calculate_phi_from_phi_to(self):
        return tx_calc.calculate_phi_constant(self._branches,self.branches_keys,self.buses_keys,ApproximationType.PTDF_LOSSES, mapping_bus_to_idx=self._busname_to_index_map)

    def _calculate_phi_loss_constant(self):
        phi_loss_from, phi_loss_to = tx_calc.calculate_phi_loss_constant(self._branches,self.branches_keys,self.buses_keys,ApproximationType.PTDF_LOSSES, mapping_bus_to_idx=self._busname_to_index_map)

        ## hold onto these for line outages
        self._phi_loss_from = phi_loss_from
        self._phi_loss_to = phi_loss_to

        ## sum the across the columns, which are indexed by branch
        phi_losses_adjust_array = phi_loss_from-phi_loss_to

        ## sum across the rows to get the total impact, and convert
        ## to dense for fast operations later
        self.phi_losses_adjust_array = phi_losses_adjust_array.sum(axis=1).T.A[0]

        ## protect the array using numpy
        self.phi_losses_adjust_array.flags.writeable = False

    def _calculate_phase_shift(self):
        
        phase_shift_array = np.fromiter(( tx_calc.calculate_susceptance(branch) * (radians(branch['transformer_phase_shift'])/branch['transformer_tap_ratio']) 
            if (branch['branch_type'] == 'transformer') 
            else 0. 
            for branch in (self._branches[bn] for bn in self.branches_keys)), float, count=len(self.branches_keys))

        self.phase_shift_array = phase_shift_array

        ## protect the array using numpy
        self.phase_shift_array.flags.writeable = False

    def _calculate_losses_phase_shift(self):

        losses_phase_shift_array = np.fromiter(( (tx_calc.calculate_conductance(branch)/branch['transformer_tap_ratio']) * radians(branch['transformer_phase_shift'])**2 
            if branch['branch_type'] == 'transformer' 
            else 0.
            for branch in (self._branches[bn] for bn in self.branches_keys)), float, count=len(self.branches_keys))

        self.losses_phase_shift_array = losses_phase_shift_array

        ## protect the array using numpy
        self.losses_phase_shift_array.flags.writeable = False

    def get_branch_ldf_iterator(self, branch_name):
        row_idx = self._branchname_to_index_map[branch_name]
        ## get the row slice
        losses_row = self.LDF[row_idx]
        yield from zip(self.buses_keys, losses_row)

    def get_branch_ldf_abs_max(self, branch_name):
        row_idx = self._branchname_to_index_map[branch_name]
        ## get the row slice
        losses_row = self.LDF[row_idx]
        return np.abs(losses_row).max()

    def get_branch_ldf_c(self, branch_name):
        return self.LDF_C[self._branchname_to_index_map[branch_name]]

    def get_branch_losses_phase_shift(self, branch_name):
        return self.losses_phase_shift_array[self._branchname_to_index_map[branch_name]]

    def get_bus_phi_losses_adj(self, bus_name):
        return self.phi_losses_adjust_array[self._busname_to_index_map[bus_name]]

    def get_branch_phi_losses_adj(self, branch_name):
        row_idx = self._branchname_to_index_map[branch_name]
        ## get the row slice
        losses_row = self.LDF[row_idx]
        return losses_row.dot(self.phi_losses_adjust_array)
