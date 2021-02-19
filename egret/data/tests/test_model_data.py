#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

import pytest
from egret.data.model_data import ModelData

testdata = {
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
                    'B1': {'bus_type': 'PV'},
                    'B2': {'bus_type': 'PQ'},
                    'B3': {'bus_type': 'PQ'}
                },
            'branch':
                {
                    'TL1': {
                        'from': 'B1',
                        'to': 'B2'
                    },
                    'TL1': {
                        'from': 'B2',
                        'to': 'B3'
                    }
                },
            'load':
                {
                    'L1': {
                        'connected_bus': 'B2',
                        'Pl': {
                            'data_type': 'time_series',
                            'values': [11.0, 111.0, 111.1]
                        },
                        'Ql': 11.0
                    }
                }
        },
    'system': { 
              'reference_bus': 'B1',
              'reference_bus_angle': 0.0,
              'time_keys' : [0.0, 1.0, 2.0],
              },
}

"""
testdata = {
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
                    'B1': {'bus_type': 'PV'},
                    'B2': {'bus_type': 'PQ'},
                    'B3': {'bus_type': 'PQ'}
                },
            'branch':
                {
                    'TL1': {
                        'from': 'B1',
                        'to': 'B2'
                    },
                    'TL1': {
                        'from': 'B2',
                        'to': 'B3'
                    }
                },
            'load':
                {
                    'L1': {
                        'connected_bus': 'B2',
                        'Pl': {
                            'data_type': 'time_series',
                            'values': [11.0, 111.0, 111.1]
                        },
                        'Ql': 11.0
                    }
                }
        },
    'system': { 
              'reference_bus': 'B1',
              'reference_bus_angle': 0.0,
              'time_keys' : [0.0, 1.0, 2.0],
              },
}
"""

def test_elements():
    md = ModelData(dict(testdata))

    for n,e in md.elements(element_type="generator"):
        assert n == 'G1' or n =='G2'
        if n == 'G1':
            assert e["generator_type"] == 'thermal'
        if n == 'G2':
            assert e["generator_type"] == 'solar'

    buses = dict(md.elements(element_type='bus'))
    assert len(buses) == 3
    assert 'B1' in buses
    assert 'B2' in buses
    assert 'B3' in buses
    assert buses['B1']['bus_type'] == 'PV'

    buses = dict(md.elements(element_type='bus', bus_type='PQ'))
    assert len(buses) == 2
    assert 'B2' in buses
    assert 'B3' in buses

def test_attributes():
    md = ModelData(dict(testdata))

    attr = md.attributes(element_type='generator')

    attr_cmp = {
        'names' : ['G1', 'G2'],
        'generator_type': {'G1': 'thermal', 'G2': 'solar'},
        'connected_bus': {'G1': 'B1', 'G2': 'B1'},
        'pg': {'G1': 100.0, 'G2': 200.0},
        'qg': {'G1': 11.0, 'G2': 22.0},
        'in_service': {'G1': True, 'G2': True}
    }
    assert attr == attr_cmp

def test_clone():
    md = ModelData(testdata)
    cmd = md.clone()

    assert md.data == cmd.data

def test_clone_at_time():
    md = ModelData(testdata)
    cloned_md = md.clone_at_time(2.0)

    comparison_md = md.clone()
    comparison_md.data['elements']['load']['L1']['Pl'] = 111.1
    del comparison_md.data['system']['time_keys']

    assert cloned_md.data == comparison_md.data

def test_clone_at_time_index():
    md = ModelData(testdata)
    cloned_md = md.clone_at_time_index(2)

    comparison_md = md.clone()
    comparison_md.data['elements']['load']['L1']['Pl'] = 111.1
    del comparison_md.data['system']['time_keys']

    assert cloned_md.data == comparison_md.data

def test_clone_at_time_keys():
    md = ModelData(testdata)
    cloned_md = md.clone_at_time_keys([1.0,2.0])

    comparison_md = md.clone()
    comparison_md.data['elements']['load']['L1']['Pl'] = \
            {'data_type':'time_series', 'values': [111.0, 111.1]}

    comparison_md.data['system']['time_keys'] = [1.0, 2.0]

    assert cloned_md.data == comparison_md.data

def test_clone_at_time_indices():
    md = ModelData(testdata)
    cloned_md = md.clone_at_time_keys([1,2])

    comparison_md = md.clone()
    comparison_md.data['elements']['load']['L1']['Pl'] = \
            {'data_type':'time_series', 'values': [111.0, 111.1]}

    comparison_md.data['system']['time_keys'] = [1.0, 2.0]

    assert cloned_md.data == comparison_md.data

def test_json_read_write():
    md = ModelData(testdata)
    md.write('testdata.json')

    md_read = ModelData.read('testdata.json')

    assert md.data == md_read.data

def test_json_gz_read_write():
    md = ModelData(testdata)
    md.write('testdata.json.gz')

    md_read = ModelData.read('testdata.json.gz')

    assert md.data == md_read.data

def test_init_read():
    md = ModelData(testdata)
    md.write('testdata.json')

    md_read = ModelData('testdata.json')

    assert md.data == md_read.data

def test_init_clone():
    md = ModelData(testdata)
    md_clone = ModelData(md)

    assert md.data == md_clone.data
    assert id(md.data) != id(md_clone.data)
