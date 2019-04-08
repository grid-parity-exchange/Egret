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

msg = 'test message'
msgcheck = b'test message\n'
blankmsgcheck = b''

def test_logging():
    # these should not show output by default
    for level in ['debug']:
        proc = subprocess.Popen(['python', 'logging_ext.py', level, msg], 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT)
        out, err = proc.communicate()
        assert out == blankmsgcheck

    # these should show output by defult
    for level in ['info', 'warning', 'error', 'critical']:
        proc = subprocess.Popen(['python', 'logging_ext.py', level, msg], 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT)
        out, err = proc.communicate()
        assert out == msgcheck



