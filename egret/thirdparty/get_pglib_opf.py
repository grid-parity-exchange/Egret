#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This script downloads the Power Grid Lib - Optimal Power Flow benchmark library
as a ZIP archive from GitHub at the following url:
https://github.com/power-grid-lib/pglib-opf/archive/master.zip

These benchmarks provide commonly used test cases for OPF formulations.
For more information see:
https://power-grid-lib.github.io

The ZIP archive is extracted to pglib-opf-master in the folder specified by the 
keyword argument download_dir (or the current working directory if left None)

This module can also be executed at the command line, and the benchmarks will
be extracted to pglib-opf-master in the current working directory.

To run the script, execute the following at a terminal prompt:
   
   > python get_pglib_opf.py

"""

import os
import logging
from zipfile import ZipFile
import pyomo.common.fileutils as futil
import pyomo.common.download as dload
# we specifically list the files to prevent extracting unexpected data
from egret.thirdparty.pglib_opf_files import pglib_files_to_extract as files

logger = logging.getLogger('egret.thirdparty.get_pglib_opf')

_pglib_zip_url = 'https://github.com/power-grid-lib/pglib-opf/archive/master.zip'

def get_pglib_opf(download_dir=None):
    print("\n\n##################################################################################\n"
                    "# This script downloads the Power Grid Lib - Optimal Power Flow benchmark library\n"
                    "# as a ZIP archive from GitHub at the following url:\n"
                    "# https://github.com/power-grid-lib/pglib-opf/archive/master.zip\n"
                    "##################################################################################\n"
)

    if download_dir is None:
        download_dir = os.path.join(os.getcwd())

    zipfile_dest = os.path.join(download_dir, 'pglib-opf-master.zip')
    pglib_dir = os.path.join(download_dir, 'pglib-opf-master')
    license_file = os.path.join(pglib_dir, 'LICENSE')

    if not os.path.isdir(download_dir):
        try:
            os.mkdir(download_dir)
        except OSError:
            logger.error("***\nFailed to create directory: {}\n"
                         "when trying to download pglib opf data\n***".format(download_dir)
                         )
            raise

    try:
        downloader = dload.FileDownloader()
        downloader.set_destination_filename(zipfile_dest)
        logger.info('... downloading from: {}'.format(_pglib_zip_url))
        downloader.get_binary_file(_pglib_zip_url)
    except:
        logger.error("***\nFailed to download: {}\n***".format(_pglib_zip_url))
        raise

    print('complete')
    print('... extracting files from zip archive')
    try:
        zf = ZipFile(zipfile_dest, 'r')
        zf.extractall(download_dir, files)
    except:
        logger.error("***\nFailed to extract files from {}\n***".format(zipfile_dest))
        raise

    print('complete')
    print('\n\n###################################################################################\n'
                    '# Successfully downloaded the Power Grid Lib - Optimal Power Flow benchmark library\n'
                    '# from: {}\n'
                    '# to the following location:\n'
                    '# {}\n'
                    '# THE LICENSE FILE FOR THESE BENCHMARKS IS SHOWN BELOW\n'
                    "# \n"
                    "# For more information see:\n"
                    "# https://power-grid-lib.github.io\n"
                    '###################################################################################\n\n'
                    .format(_pglib_zip_url, pglib_dir))


    # show the license file and the download message
    with open(license_file, 'r') as fd:
        print(fd.read())
    print('\n\n####################################################################################\n'
                    '# Successfully downloaded the Power Grid Lib - Optimal Power Flow benchmark library\n'
                    '# from: {}\n'
                    '# to the following location:\n'
                    '# {}\n'
                    '# THE LICENSE FILE FOR THESE BENCHMARKS IS SHOWN ABOVE\n'
                    "# \n"
                    "# For more information see:\n"
                    "# https://power-grid-lib.github.io\n"
                    '####################################################################################\n\n'
                    .format(_pglib_zip_url, pglib_dir))

if __name__ == '__main__':
    get_pglib_opf()

