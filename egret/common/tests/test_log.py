#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module tests the egret logger. Note that the pytest infrastructure captures log output,
therefore, we make an external call to bypass the pytest capturing
"""
import subprocess
import os

msg = 'test message'
msgcheck = b'test message\n'
blankmsgcheck = b''

test_dir = os.path.dirname(os.path.realpath(__file__))

def run_logging_ext(level):
    proc = subprocess.Popen(['python', 
                             os.path.join(test_dir, 'logging_ext.py'),
                             level, msg], 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.STDOUT)
    return proc.communicate()

def _test_logging():
    # these should not show output by default
    for level in ['debug']:
        out, err = run_logging_ext(level)
        assert out == blankmsgcheck

    # these should show output by defult
    for level in ['info', 'warning', 'error', 'critical']:
        out, err = run_logging_ext(level)
        assert out == msgcheck



