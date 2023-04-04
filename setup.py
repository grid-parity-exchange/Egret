#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

from setuptools import setup, find_packages
from pathlib import Path

DISTNAME = 'gridx-egret'
VERSION = '0.5.6.dev0'
PACKAGES = find_packages()
EXTENSIONS = []
DESCRIPTION = 'EGRET: Electrical Grid Research and Engineering Tools.'
AUTHOR = 'Michael Bynum, Anya Castillo, Carl Laird, Bernard Knueven and Jean-Paul Watson'
MAINTAINER_EMAIL = 'carldlaird@users.noreply.github.com'
LICENSE = 'Revised BSD'
URL = 'https://github.com/grid-parity-exchange/Egret'

setuptools_kwargs = {
    'zip_safe': False,
    'scripts': [],
    'include_package_data': True,
    'install_requires': ['pyomo>=6.4', 'numpy', 'pytest', 'pandas',
                         'matplotlib', 'seaborn', 'scipy', 'networkx',
                         'coramin==0.1.1'],
    'python_requires' : '>=3.7, <4',
}

this_directory = Path(__file__).parent
long_description = (this_directory / 'README.md').read_text()

setup(name=DISTNAME,
      version=VERSION,
      packages=PACKAGES,
      ext_modules=EXTENSIONS,
      description=DESCRIPTION,
      author=AUTHOR,
      maintainer_email=MAINTAINER_EMAIL,
      license=LICENSE,
      long_description=long_description,
      long_description_content_type='text/markdown',
      url=URL,
      **setuptools_kwargs)
