![EGRET GitHub CI](https://github.com/grid-parity-exchange/Egret/workflows/EGRET%20GitHub%20CI/badge.svg)

## EGRET Overview

EGRET is a Python-based package for electrical grid optimization based on the Pyomo optimization modeling language. EGRET is designed to be friendly for performing high-level analysis (e.g., as an engine for solving different optimization formulations), while also providing flexibility for researchers to rapidly explore new optimization formulations.

Major features:
* Solution of Unit-Commitment problems
* Solution of Economic Dispatch (optimal power flow) problems (e.g., DCOPF, ACOPF)
* Library of different problem formulations and approximations
* Generic handling of data across model formulations
* Declarative model representation to support formulation development

EGRET is available under the BSD License (see LICENSE.txt)

### Installation

* EGRET is a Python package and therefore requires a Python installation. We recommend using Anaconda with the latest Python (https://www.anaconda.com/distribution/).
* These installation instructions assume that you have a recent version of Pyomo installed, in addition to a suite of relevant solvers (see www.pyomo.org for additional details).
* Download (or clone) EGRET from this GitHub site.
* From the main EGRET folder (i.e., the folder containing setup.py), use a terminal (or the Anaconda prompt for Windows users) to run setup.py to install EGRET into your Python installation - as follows:

   pip install -e .

### Requirements

* Python 3.7 or later
* Pyomo version 5.7.1 or later
* pytest
* Optimization solvers for Pyomo - specific requirements depends on the models being solved. EGRET is tested with Gurobi or CPLEX for MIP-based problems (e.g., unit commitment) and Ipopt (with HSL linear solvers) for NLP problems.

We additionally recommend that EGRET users install the open source CBC MIP solver. The specific mechanics of installing CBC are platform-specific. When using Anaconda on Linux and Mac platforms, this can be accomplished simply by:

   conda install -c conda-forge coincbc

The COIN-OR organization - who developers CBC - also provides pre-built binaries for a full range of platforms on https://bintray.com/coin-or/download.

### Testing the Installation

To test the functionality of the unit commitment aspects of EGRET, execute the following command from the EGRET models/tests sub-directory:

   pytest test_unit_commitment.py

If EGRET can find a commerical MIP solver on your system via Pyomo, EGRET will execute a large test suite including solving several MIPs to optimality. If EGRET can only find an open-source solver, it will execute a more limited test suite which mostly relies on solving LP relaxations. Example output is below.

=================================== test session starts ==================================<br/>
platform darwin -- Python 3.7.7, pytest-5.4.2, py-1.8.1, pluggy-0.13.0<br/>
rootdir: /home/some-user/egret<br/>
collected 21 items<br/>
<br/>
test_unit_commitment.py s....................                                       [100%]<br/>
<br/>
========================= 20 passed, 1 skipped in 641.80 seconds =========================<br/>


### How to Cite EGRET in Your Research

If you are using the unit commitment functionality of EGRET, please cite the following paper: 

On Mixed-Integer Programming Formulations for the Unit Commitment Problem
Bernard Knueven, James Ostrowski, and Jean-Paul Watson.
INFORMS Journal on Computing (Ahead of Print)
https://pubsonline.informs.org/doi/10.1287/ijoc.2019.0944
