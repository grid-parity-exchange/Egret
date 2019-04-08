import pytest
import os
from pyomo.common.fileutils import this_file
from egret.thirdparty.get_pglib import get_pglib

_this_file = this_file()

# not sure we want to download files during tests
# also, this will be tested by travis since it is needed for other tests
@pytest.mark.skip('Skip testing the downloads since they will already be tested in travis')
def test_get_pglib():
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
    get_pglib(downloads_dir)
    assert os.path.isfile(os.path.join(pglib_dir, 'LICENSE'))

    
