"""
This script downloads the Power Grid Lib - Optimal Power Flow benchmark library
as a ZIP archive from GitHub at the following url:
https://github.com/power-grid-lib/pglib-opf/archive/master.zip

These benchmarks provide commonly used test cases for OPF formulations.
For more information see:
https://power-grid-lib.github.io

The ZIP archive is extracted and placed in the following folder relative to the egret
repository location:

./egret/thirdparty/downloads/pglib-opf-master

To run the script, execute the following at a terminal prompt:
   
   > python get_pglib.py

"""

import os
import logging
from zipfile import ZipFile
import pyomo.common.fileutils as futil
import pyomo.common.download as dload
from egret.thirdparty.pglib_files import pglib_files_to_extract as files

logger = logging.getLogger('egret.thirdparty.get_pglib')

_pglib_zip_url = 'https://github.com/power-grid-lib/pglib-opf/archive/master.zip'

_this_file = futil.this_file()
_thirdparty_dir = os.path.dirname(_this_file)
_download_dir = os.path.join(_thirdparty_dir, 'downloads')
_zipfile_dest = os.path.join(_download_dir, 'pglib-master.zip')
_pglib_dir = os.path.join(_download_dir, 'pglib-opf-master')
_license_file = os.path.join(_pglib_dir, 'LICENSE')

def get_pglib():
    logger.critical("\n\n##################################################################################\n"
                    "# This script downloads the Power Grid Lib - Optimal Power Flow benchmark library\n"
                    "# as a ZIP archive from GitHub at the following url:\n"
                    "# https://github.com/power-grid-lib/pglib-opf/archive/master.zip\n"
                    "##################################################################################\n"
)

    if not os.path.isdir(_download_dir):
        try:
            os.mkdir(_download_dir)
        except OSError:
            logger.error("***\nFailed to create directory: {}\n"
                         "when trying to download pglib data\n***".format(_download_dir)
                         )
            raise

    try:
        downloader = dload.FileDownloader()
        downloader.set_destination_filename(_zipfile_dest)
        logger.info('... downloading from: {}'.format(_pglib_zip_url))
        downloader.get_binary_file(_pglib_zip_url)
    except:
        logger.error("***\nFailed to download: {}\n***".format(_pglib_zip_url))
        raise

    logger.critical('complete')
    logger.critical('... extracting files from zip archive')
    try:
        zf = ZipFile(_zipfile_dest, 'r')
        zf.extractall(_download_dir, files)
    except:
        logger.error("***\nFailed to extract files from {}\n***".format(_zipfile_dest))
        raise

    logger.critical('complete')
    logger.critical('\n\n###################################################################################\n'
                    '# Successfully downloaded the Power Grid Lib - Optimal Power Flow benchmark library\n'
                    '# from: {}\n'
                    '# to the following location:\n'
                    '# {}\n'
                    '# THE LICENSE FILE FOR THESE BENCHMARKS IS SHOWN BELOW\n'
                    "# \n"
                    "# For more information see:\n"
                    "# https://power-grid-lib.github.io\n"
                    '###################################################################################\n\n'
                    .format(_pglib_zip_url, _pglib_dir))


    # show the license file and the download message
    with open(_license_file, 'r') as fd:
        logger.critical(fd.read())
    logger.critical('\n\n####################################################################################\n'
                    '# Successfully downloaded the Power Grid Lib - Optimal Power Flow benchmark library\n'
                    '# from: {}\n'
                    '# to the following location:\n'
                    '# {}\n'
                    '# THE LICENSE FILE FOR THESE BENCHMARKS IS SHOWN ABOVE\n'
                    "# \n"
                    "# For more information see:\n"
                    "# https://power-grid-lib.github.io\n"
                    '####################################################################################\n\n'
                    .format(_pglib_zip_url, _pglib_dir))

if __name__ == '__main__':
    get_pglib()

