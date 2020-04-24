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
from egret.data.tests.test_model_data import testdata
from egret.model_library.unit_commitment.uc_utils import uc_time_helper

## these should be arbitary to mimic 
## a pyomo RangeSet
TimePeriods = [3,4,5]

time_mapper = uc_time_helper(TimePeriods)

md_testdata = ModelData(testdata)
load_attrs = md_testdata.attributes(element_type='load')

def test_None():
    assert time_mapper(None) == dict()

def test_empty_dict():
    assert time_mapper(dict()) == dict()

def test_single_item_mapping():
    expected_result = { 3:11.0, 4:111.0, 5:111.1 }
    assert time_mapper(load_attrs['Pl']['L1']) == expected_result

def test_single_item_expansion():
    expected_result = { 3:11.0, 4:11.0, 5:11.0 }
    assert time_mapper(load_attrs['Ql']['L1']) == expected_result

def test_multi_item_mapping():
    expected_result = { ('L1', 3):11.0, ('L1', 4):111.0, ('L1',5):111.1 }
    assert time_mapper(load_attrs['Pl']) == expected_result

def test_multi_item_expansion():
    expected_result = { ('L1',3):11.0, ('L1',4):11.0, ('L1',5):11.0 }
    assert time_mapper(load_attrs['Ql']) == expected_result
