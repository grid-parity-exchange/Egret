#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

from enum import Enum

class FlowType(Enum):
    CURRENT = 1
    POWER = 2

class CoordinateType(Enum):
    POLAR = 1
    RECTANGULAR = 2

class ApproximationType(Enum):
    BTHETA = 1
    BTHETA_LOSSES = 2
    PTDF = 3
    PTDF_LOSSES = 4

class BasePointType(Enum):
    FLATSTART = 1
    SOLUTION = 2

class RelaxationType(Enum):
    NONE = 1
    SOC = 2
