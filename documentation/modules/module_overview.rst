.. _modules-overview:

================
Modules Overview
================

The modules are the building blocks to create an IDSimF simulation application. Every module opens its own namespace. There are the following modules in the `modules` directory: 


:doc:`Core <core>`: Core data structures and base facilities 
    Bundles basic data structures like vectors and basic facilities required in most other modules like random generators.

:doc:`SpaceCharge <space_charge>`: Space Charge simulation
    Provides algorithms / methods to simulate electric particle-particle interaction (space charge). 

:doc:`Integration <integration>`: Numerical integration of particle motion
    Provides algorithms / methods to simulate electric particle-particle interaction (space charge). 

:doc:`ParticleSimulation <particlesimulation>`: Facilities for particle trajectory simulation
    Provides algorithms, data structures and methods for charged particle trajectory simulations. This module contains for example trajectory integrators, file readers and writers, potential arrays and interpolated scalar and vector fields. 

:doc:`CollisionModel <collisionmodels>`: Background gas collision modeling
    Provides models to describe the collisional interaction between charged particles and neutral background gas. 

:doc:`RS <rs>`: Reaction simulation and particle chemical kinetics modeling
    Provides particle based modeling of chemical kinetics and reaction dynamics of charged particles. 

There is an additional application utilities module :cpp:any:`AppUtils` in `applications/util` which bundles supporting building blocks for actual simulation applications: 

:doc:`AppUtils <apputils>`: Application utils for building of actual simulation applications

Parts of the in source docstrings are not finished yet. In the meantime, refer to the actual source code for details in such cases. 