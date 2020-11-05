.. _modules-overview:

================
Modules Overview
================

The modules are the building blocks to create an IDSimF simulation application. Every module opens its own namespace. There are the following modules: 

``Core``: Core data structures and base facilities 
    Bundles basic data structures like vectors and basic facilities required in most other modules like random generators.

``BTree``: Barnes-Hut tree
    Provides Barnes-Hut tree implementations which are used to simulate electric particle-particle interaction (space charge). 

``ParticleSimulation``: Facilities for particle trajectory simulation
    Provides algorithms, data structures and methods for charged particle trajectory simulations. This module contains for example trajectory integrators, file readers and writers, potential arrays and interpolated scalar and vector fields. 

``CollisionModel``: Background gas collision modeling
    Provides models to describe the collisional interaction between charged particles and neutral background gas. 

``RS``: Reaction simulation and particle chemical kinetics modeling
    Provides particle based modeling of chemical kinetics and reaction dynamics of charged particles. 


Parts of the in source docstrings are not finished yet. In the meantime, refer to the actual source code for details. 