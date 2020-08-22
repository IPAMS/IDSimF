# Ion Dynamics Simulation Framework - IDSimF

IDSimF, the *Ion Dynamics Simulation Framework*, is an open source framework for the simulation of non relativistic dynamics of molecular ions. The primary application of IDSimF is in the domain of mass spectrometry and ion mobility spectrometry.

## Overview

IDSimF is a collection of modules and programs which allows to simulate the trajectories of molecular ions at non relativistic conditions. It provides individual solver programs for individual simulation tasks and a modular framework for the quick development of new solvers. 

IDSimF was inspired by other simulation frameworks of a similar structure, particularly the open source fluid dynamics framework [OpenFOAM](https://openfoam.org/). IDSimF is developed primarily in C++. 

The primary application of IDSimF is the simulation of the dynamics of molecular ions in the domain of Mass Spectrometry (MS), Ion Mobility Spectrometry (IMS) and similar methods.  

### Features

IDSimF currently provides as main features: 

* Charged particle trajectory integration in electric fields
* Import of electrostatics solutions and electric fields from other codes (e.g. [Comsol Multiphysics](https://www.comsol.com/) and [SIMION](https://simion.com/>))
* Electric particle - particle interaction (space charge)
* Interaction of ions with neutral background gas
* Particle based chemical kinetics / simulation of reactive ions

### What IDSimF is not

IDSimF is developed according to the requirements in molecular mass spectrometry. It is therefore *not at all* a full featured plasma code: 

* Currently only electrostatics are considered, no magnetic fields or magnetic couplings between particles are calculated. However, simple magnetic interactions (e.g. for ion cyclotron resonance simulations) can be added by explicit force terms in the trajectory integration.
* Oscillating / RF-Fields are simulated with an electrostatic approximation. Wave phenomena / RF radiation effects are not simulated. 

### Contributers

Currently, the main contributor to IDSimF is the [Institute for Pure and Applied Mass Spectrometry (IPAMS)](https://www.ipams.uni-wuppertal.de/en/>) and the [Physical and Theroetical Chemistry](https://www.ptc.uni-wuppertal.de/)  of the [University of Wuppertal](https://www.uni-wuppertal.de/). 

## Requirements

IDSimF requires a C++17 compatible C++ compiler and [HDF5](https://www.hdfgroup.org/solutions/hdf5) C++ libraries from the HDF group. It uses [CMake](https://cmake.org) as build tool.

## Installation 

Refer to the "Installation" section of the documentation for details how to compile and install IDSimF. 

## Documentation

The documentation and user guide is available online at https://idsimf.readthedocs.io/ 

The main documentation is written with [Sphinx](https://www.sphinx-doc.org/). It is located in `/documentation` and can be compiled to html locally with `make html`. See the documentation for details of the simulation applications and modules. 

## IDSimPy

The analysis of IDSimF results can be performed with [IDSimPy](https://github.com/IPAMS/IDSimPy). IDSimPy is a separated Python package for the preprocessing of IDSimF input and postprocessing of IDSimF result data.

## License
License: IDSimF is licensed under GNU General Public License v3.0, see LICENSE


