#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

from setuptools import setup, find_packages
from distutils.core import Extension

DISTNAME = 'gridx-egret'
VERSION = '0.1.0'
PACKAGES = find_packages()
EXTENSIONS = []
DESCRIPTION = 'EGRET: Electrical Grid Research and Engineering Tools.'
LONG_DESCRIPTION = open('README.md').read()
AUTHOR = 'Michael Bynum, Anya Castillo, Carl Laird, Bernard Knueven and Jean-Paul Watson'
MAINTAINER_EMAIL = 'carldlaird@users.noreply.github.com'
LICENSE = 'Revised BSD'
URL = 'https://github.com/grid-parity-exchange/Egret'

setuptools_kwargs = {
    'zip_safe': False,
    'install_requires': [],
    'scripts': [],
    'include_package_data': True,
    'install_requires' : ['pyomo>=5.7.1', 'numpy', 'pytest', 'pandas', \
                            'matplotlib', 'seaborn', 'scipy', 'networkx'],
    'python_requires' : '>=3.7, <4',
}

setup(name=DISTNAME,
      version=VERSION,
      packages=PACKAGES,
      ext_modules=EXTENSIONS,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      maintainer_email=MAINTAINER_EMAIL,
      license=LICENSE,
      url=URL,
      **setuptools_kwargs)
