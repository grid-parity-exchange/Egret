#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

import sys
import logging
import egret.models.acopf

logger = logging.getLogger('egret.models.acopf')

if sys.argv[1] == 'debug':
    logger.debug(sys.argv[2])
elif sys.argv[1] == 'info':
    logger.info(sys.argv[2])
elif sys.argv[1] == 'warning':
    logger.warning(sys.argv[2])
elif sys.argv[1] == 'error':
    logger.error(sys.argv[2])
elif sys.argv[1] == 'critical':
    logger.critical(sys.argv[2])

