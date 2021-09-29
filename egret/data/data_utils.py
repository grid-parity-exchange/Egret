#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module defines some utilities for handling ModelData dictionaries
"""
import copy as cp

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

def _recurse_into_time_index(old_node, time_index):
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
                new_node[key] = _recurse_into_time_index(att,time_index)
        else:
            # be paranoid about other attributes (could be list, or other mutable type)
            new_node[key] = cp.deepcopy(att)
    return new_node

def _recurse_into_time_indices(old_node, time_indices):
    new_node = dict()
    for key, att in old_node.items():
        if isinstance(att, dict):
            if 'data_type' in att and att['data_type'] == 'time_series':
                vals = att['values']
                new_node[key] = { 'data_type': 'time_series',
                                  'values' : [vals[i] for i in time_indices] }
            else:
                new_node[key] = _recurse_into_time_indices(att,time_indices)
        else:
            new_node[key] = cp.deepcopy(att)
    return new_node

def _get_sub_list_indicies(master_list, sub_list):
    '''
    Finds the indices of the elements in sub_list in master_list and
    returns them as a list. Optimized for when the elements in sub_list
    are in the order they appear in master_list.
    '''
    sub_list_pos = 0
    master_list_pos = 0
    len_sub_list = len(sub_list)
    sub_index_list = list()

    while sub_list_pos < len_sub_list:
        begin_sub_list_pos = sub_list_pos
        for idx, val in enumerate(master_list):
            if val == sub_list[sub_list_pos]:
                sub_index_list.append(idx)
                sub_list_pos += 1
            if sub_list_pos >= len_sub_list:
                break
        if begin_sub_list_pos == sub_list_pos:
            raise Exception("Could not find element {} in the list {}".format(sub_list[sub_list_pos], master_list))
    return sub_index_list

def _read_from_file(filename, file_type):
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
    return data
