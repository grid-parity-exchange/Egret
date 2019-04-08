#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This is the logging configuration for EGRET.

The documentation below is primarily for EGRET developers.

Examples
========
To use the logger in your code, add the following 
after your import
.. code-block:: python
   
   import logging
   logger = logging.getLogger('egret.path.to.module')

Then, you can use the standard logging functions
.. code-block:: python
   
   logger.debug('message')
   logger.info('message')
   logger.warning('message')
   logger.error('message')
   logger.critical('message')
   
Note that by default, any message that has a logging level
of warning or higher (warning, error, critical) will be
logged.

To log an exception and capture the stack trace
.. code-block:: python

   try:
      c = a / b
   except Exception as e:
      logging.error("Exception occurred", exc_info=True)

"""
import logging
log_format = '%(message)s'
logging.basicConfig(level=logging.DEBUG,
                    format=log_format)

