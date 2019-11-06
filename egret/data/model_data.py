#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
The data used for building models in EGRET is stored in a nested Python dictionary.
The format for this dictionary is as follows:

.. code-block:: python
   :linenos:

    {
    'elements':
        {
        <element-type>:
            {
            <element-name>:
                {
                <attribute-1>: <value-1>,
                <attribute-2>: <value-2>,
                ...
                }
            }
        }
    'system':
        {
        <attribute-1>, <value-1>
        }
    }

* The main level must contain the following keys: "elements" and "system".

* data["elements"] is a dictionary where the keys are the element types (e.g. "generator",
  "bus", "load"). Each of these is a dictionary where the keys are the model element names.
  Each of these is a dictionary containing the attributes and values.

* Attributes on modeling elements are typically floats, but can also be any desired
  dictionary-based data structure. There are some specialized attribute values that
  are recognized if the attribute is itself a dictionary with the key "data_type".
  Some supported data types are "time_series" and "cost_curve".

* There are a number of data types in addition to native python data.
  if the "data_type" is equal to "time_series", then the "values" key is interpreted
  as a dictionary that contains key, value pairs where the keys are the timestamps.
  These timestamps could be used to create a slice of the data object at a particular
  time. (See :py:meth:GridData.`slice_at_index`)

* data["system"] is intended to store additional system-wide requirements
  which are not part of the physical infrastructure or modeling elements.

For an example of the dictionary data structure, see the following:

.. code-block:: python
   :linenos:

    {
    'elements':
        {
        'generator':
            {
            'G1': {
                  'generator_type': 'thermal',
                  'connected_bus': 'B1',
                  'pg': 100.0,
                  'qg': 11.0,
                  'in_service': True
                  },
            'G2': {
                  'generator_type': 'solar',
                  'connected_bus': 'B1',
                  'pg': 200.0,
                  'qg': 22.0,
                  'in_service': True
                  }
             },
        'bus':
            {
            'B1': {'bus_type': 'PV'}
            'B2': {'bus_type': 'PQ'}
            'B3': {'bus_type': 'PQ'}
            }
        'branch':
            {
            'TL1': {
                   'from': 'B1',
                   'to': 'B2'
                   }
            'TL1': {
                   'from': 'B2',
                   'to': 'B3'
                   }
            }
        'load':
            {
            'L1': {
                  'connected_bus': 'B2',
                  'Pl': {
                      'data_type':'time_series',
                      'values': {0.0: 11.0, 1.0:111.0, 2.0:111.1}
                      },
                  'Ql': 11.0
                }
            }
        },
    'system': { 
              'reference_bus': 'B1',
              'reference_bus_angle': 0.0,
              },
    }

This module provides the helper class :py:class:`ModelData` that
provides utilities for interacting and manipulating this
nested dictionary structure.

.. todo::

    * Document the data-types for attributes

"""
import copy as cp
import logging
logger = logging.getLogger('egret.model_data')

class ModelData(object):
    @staticmethod
    def empty_model_data_dict():
        """
        Create an empty model_data dictionary with the high-level "elements"
        and "system" keys.

        Returns
        -------
            dict : a new model_data dictionary
        """
        return {"elements": dict(), "system": dict()}
    
    def __init__(self, data=None):
        """
        Create a new ModelData object to wrap a model_data dictionary with some helper methods.

        Parameters
        ----------
        data : dict or None
           An initial model_data dictionary if it is available, otherwise, a new model_data
           dictionary is created.
        """
        if data:
            self.data = data
        else:
            self.data = ModelData.empty_model_data_dict()

    @classmethod
    def read(cls, filename, file_type=None):
        """
        Reads data from a file into a new ModelData object

        Parameters
        ----------
        filename : str
            The path to the file
        file_type : None,str (optional)
            The specification of the file_type. Valid values are 'json', 'json.gz' for json-ed
            EGRET ModelData objects, 'm' for MATPOWER files, 'dat' for Prescient data files, and
            'pglib-uc' for json files from pglib-uc. If None, the file type is inferred from the
            extension.
        """
        valid_file_types = ['json', 'json.gz', 'm', 'dat', 'pglib-uc']
        if file_type is not None and file_type not in valid_file_types:
            raise Exception("Unrecognized file_type {}. Valid file types are {}".format(file_type, valid_file_types))
        elif file_type is None:
            ## identify the file type
            if filename[-5:] == '.json':
                file_type = 'json'
            elif filename[-8:] == '.json.gz':
                file_type = 'json.gz'
            elif filename[-2:] == '.m':
                file_type = 'm'
            elif filename[-4:] == '.dat':
                file_type = 'dat'
            else:
                raise Exception("Could not infer type of file {} from its extension!".format(filename))

        if file_type == 'json':
            import json
            with open(filename) as f:
                data = json.load(f)
        elif file_type == 'json.gz':
            import json
            import gzip
            with gzip.open(filename, 'rt') as f:
                data = json.load(f)
        elif file_type == 'm':
            from egret.parsers.matpower_parser import create_model_data_dict
            data = create_model_data_dict(filename)
        elif file_type == 'dat':
            from egret.parsers.prescient_dat_parser import create_model_data_dict
            data = create_model_data_dict(filename)
        elif file_type == 'pglib-uc':
            from egret.parsers.pglib_uc_parser import create_model_data_dict
            data = create_model_data_dict(filename)

        logger.debug("ModelData read from {}".format(filename))

        return cls(data=data)

    def elements(self, element_type, **kwargs):
        """
        A python generator that loops over modeling elements of a particular element type
        and, if requested, other sub-attributes

        This method provides a generator that loops over the model elements (dictionaries)
        defined in model_data["elements"][element_type]. If additional arguments are provided,
        then those subattributes are also tested.

        Parameters
        ----------
        element_type : str
           Desired element type.
        **kwargs : key=str named arguments
           Additional arguments provides key=value pairs to test on each element.

        Returns
        -------
            generator

        .. todo::
           * need a better error message when element_type is not found

        """
        if element_type not in self.data['elements'].keys():
            return

        if len(kwargs) == 0:
            for name, elem in self.data['elements'][element_type].items():
                yield name, elem
            return

        # additional attributes have been specified
        for name, elem in self.data['elements'][element_type].items():
            include = True
            for k, v in kwargs.items():
                if k not in elem or elem[k] != v:
                    include = False
                    break
            if include:
                yield name, elem

    def attributes(self, element_type, **kwargs):
        """
        Returns a dictionary arranged by attribute -- element-name -- value instead of
        element-name -- attribute -- value.

        This method loops over the modeling elements of a particular element type
        and, if requested, filtes using other sub-attributes. While the original 
        structure follows:

        >>> data['elements'][<element-name>][<attribute>] = <value>  # doctest: +SKIP

        the returned dictionary follows:

        >>> d[<attribute>][<element-name>] = <value>  # doctest: +SKIP

        Note that this call is actually building a dictionary (not a generator).

        Parameters
        ----------
        element_type : str
           Desired element type.
        **kwargs : key=str named arguments
           Additional arguments provides key=value pairs to test on each element.

        Returns
        -------
            dict

        .. todo::
           * need a better error message when element_type is not found

        """
        if element_type not in self.data['elements']:
            return None

        retdict = dict()
        retdict['names'] = list()

        for name, elem in self.elements(element_type=element_type, **kwargs):
            retdict['names'].append(name)
            for attrib, value in elem.items():
                if attrib not in retdict:
                    retdict[attrib] = dict()
                retdict[attrib][name] = value

        return retdict

    def clone(self):
        """
        Create a copy of this ModelData object using a deep copy on the underlying dictionary

        Returns
        -------
            ModelData
        """
        return ModelData(cp.deepcopy(self.data))

    def clone_in_service(self):
        """
        Create a copy of this ModelData object using a deep copy on the underlying dictionary,
        only returning the elements for which the in_service flag is not set to False

        Returns
        -------
            ModelData
        """
        return ModelData(_copy_only_in_service(self.data))

    def clone_at_timestamp(self, timestamp):
        """
        Creae a copy of the ModelData object using values from a single timestamp.

        Create a copy of this ModelData object, but with the following change. Whenever
        a time_series is encountered (recognized by "data_type"="time_series"), the
        attribute containing the time series is replaced with a float value corresponding
        to the specified timestamp.

        .. todo::
           This can actually be made general based on "data_type" and idx instead of specific to time_series

        Parameters
        ----------
        timestamp : str, int, or TimeStamp
           The desired timestamp to use when collapsing the representation.


        Returns
        -------
            ModelData
        """
        time_indices = self.data['system']['time_indices']

        time_index = time_indices.index(timestamp)

        gd = self.clone_at_timeindex(time_index)

        return gd

    def clone_at_timeindex(self, time_index):
        """
        Creae a copy of the ModelData object using values from a single time index value.

        Create a copy of this ModelData object, but with the following change. Whenever
        a time_series is encountered (recognized by "data_type"="time_series"), the
        attribute containing the time series is replaced with a float value corresponding
        to the specified time index.

        .. todo::
           This can actually be made general based on "data_type" and idx instead of specific to time_series

        Parameters
        ----------
        time_index : int
            The index into time series data at which to copy


        Returns
        -------
            ModelData
        """
        mdclone = ModelData(self._recurse_into_timestamp(self.data, time_index))

        ## the new model data has no time
        del mdclone.data['system']['time_indices']

        return mdclone

    def write(self, filename, file_type=None):
        """
        Dumps the ModelData object dict to the specified file.
        Optionally, the file_type can be specified if not inferred.

        Parameters
        ----------
        filename : str
            The full filename including extension and path.
        file_type : None
            If specified, the encoding used when writing the file.
            If None, it will be inferred from the file extension.
        """
        valid_file_types = ['json', 'json.gz']
        if file_type is not None and file_type not in valid_file_types:
            raise Exception("Unrecognized file_type {}. Valid file types are {}".format(file_type, valid_file_types))
        elif file_type is None:
            if filename[-5:] == '.json':
                file_type = 'json'
            elif filename[-8:] == '.json.gz':
                file_type = 'json.gz'
            else:
                logger.warning("Unrecognized file_type for file {} in ModelData.write, using 'json'".format(filename))
                file_type = 'json'

        if file_type == 'json':
            import json
            with open(filename,'w') as f:
                json.dump(self.data, f)
        elif file_type == 'json.gz':
            import json
            import gzip
            with gzip.open(filename, 'wt') as f:
                json.dump(self.data, f)
        logger.debug("ModelData written to {}".format(filename))

    def _recurse_into_timestamp(self, old_node, time_index):
        # create a new node for the new dict
        new_node = dict()
        # loop of the exisiting attributes
        for key, att in old_node.items():
            # ignore if not at dict
            if isinstance(att, dict):
                if 'data_type' in att and att['data_type'] == 'time_series':
                    vals = att['values']
                    new_node[key] = att['values'][time_index]
                else:
                    new_node[key] = self._recurse_into_timestamp(att,time_index)
            else:
                # be paranoid about other attributes (could be list, or other mutable type)
                new_node[key] = cp.deepcopy(att)
        return new_node


# TODO: These should be moved to a more general "model utilities" module
def map_items(func, d):
    return {k: func(v) for k, v in d.items()}


def zip_items(dict_lb, dict_ub):
    return {k: (dict_lb[k], dict_ub[k]) for k in dict_lb.keys()}

def _copy_only_in_service(data_dict):
    new_dd = dict()
    for key, value in data_dict.items():
        if key == 'elements':
            ## value is the elements dictionary
            new_dd[key] = dict()
            new_elements = new_dd[key]
            for elements_name, elements in value.items():
                new_elements[elements_name] = dict()
                new_element_dict = new_elements[elements_name]
                for element_name, element in elements.items():
                    if 'in_service' in element and (not element['in_service']):
                        continue
                    else:
                        new_element_dict[element_name] = cp.deepcopy(element)
        else:
            new_dd[key] = cp.deepcopy(value)
    return new_dd

