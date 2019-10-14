#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

import pytest
import os
from pyomo.common.fileutils import this_file
from egret.thirdparty.get_pglib_opf import get_pglib_opf
from egret.thirdparty.get_pglib_uc import get_pglib_uc

_this_file = this_file()

# not sure we want to download files during tests
# also, this will be tested by travis since it is needed for other tests
@pytest.mark.skip('Skip testing the downloads since they will already be tested in travis')
def test_get_pglib_opf():
    # go one level before the repository folder
    # Egret-1/egret/thirdparty/tests
    repos_parent_path = \
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(
                    os.path.dirname(
                        os.path.dirname(_this_file)
                        )
                    )
                )
            )
    
    downloads_dir = os.path.join(repos_parent_path, 'downloads_test')
    pglib_dir = os.path.join(downloads_dir, 'pglib-opf-master')
    get_pglib_opf(downloads_dir)
    assert os.path.isfile(os.path.join(pglib_dir, 'LICENSE'))

    
@pytest.mark.skip('Skip testing the downloads since they will already be tested in travis')
def test_get_pglib_uc():
    # go one level before the repository folder
    # Egret-1/egret/thirdparty/tests
    repos_parent_path = \
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(
                    os.path.dirname(
                        os.path.dirname(_this_file)
                        )
                    )
                )
            )
    
    downloads_dir = os.path.join(repos_parent_path, 'downloads_test')
    pglib_dir = os.path.join(downloads_dir, 'pglib-uc-master')
    get_pglib_uc(downloads_dir)
    assert os.path.isfile(os.path.join(pglib_dir, 'LICENSE'))

