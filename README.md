## EGRET Overview

EGRET is a Python-based package for electrical grid optimization based on the Pyomo optimization modeling language. EGRET is designed to be friendly for performing high-level analysis (e.g., as an engine for solving different optimization formulations), while also providing flexibility for researchers to rapidly explore new optimization formulations.

Major features:
* Solution of Unit-Commitment problems
* Solution of Economic Dispatch (optimal power flow) problems (e.g., DCOPF, ACOPF)
* Library of different problem formulations and approximations
* Generic handling of data across model formulations
* Declarative model representation to support formulation development

EGRET is available under the BDS License (see LICENSE.txt)

### Installation

* EGRET is a Pythonb package and requires a python installation. We recommend using Anaconda with the latest Python (https://www.anaconda.com/distribution/).
* These installation instructions assume that you have a recent version of Pyomo installed with the required solvers (see www.pyomo.org).
* Download (or clone) EGRET from this GitHub site
* From the main egret folder (folder containing setup.py), use a terminal (or the Anaconda prompt for Windows users), and run setup.py to install EGRET into your Python installation.

   python setup.py install

### Requirements
* Pyomo version 5.6 or later
* Optimization solvers for Pyomo - specific requirements depends on the models being solved, however, EGRET is tested with GUROBI or CPLEX for MIP-based problems (e.g., unit commitment) and IPOPT (with HSL linear solvers) for NLP problems.
