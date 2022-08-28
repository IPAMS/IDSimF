.. _modules-particlesimulation:

==========================
Module: particleSimulation
==========================

The `ParticleSimulation` module provides algorithms, data structures, input/output and data transformation methods required for actual charged particle trajectory simulations. 


Data Structures / Simulation Objects
====================================

.. doxygenclass:: ParticleSimulation::InterpolatedField
    :members:
    :undoc-members:

.. doxygenclass:: ParticleSimulation::SampledWaveform
    :members:
    :undoc-members:

.. doxygenclass:: ParticleSimulation::SimionPotentialArray
    :members:
    :undoc-members:




Particle Start Zones
====================

Particles / Ions are started from particle start zones. :cpp:class:`ParticleSimulation::ParticleStartZone` is the abstract super class of all particle start zones. All particle start zones can generate random start positions in the start zone with the :cpp:any:'getRandomParticlePosition' method. A random set of particles in the particle start zone can be generated with the :cpp:any:`getRandomParticlesInStartZone` method.

Currently, there are two instantiable particle start zones, a box with faces parallel to the main axes :cpp:class:`ParticleSimulation::BoxStartZone`, and a cylindrical start zone which can be rotated and shifted :cpp:class:`ParticleSimulation::CylinderStartZone`. 

.. doxygenclass:: ParticleSimulation::ParticleStartZone
    :members:
    :undoc-members:

.. doxygenclass:: ParticleSimulation::BoxStartZone
    :members:
    :undoc-members:

.. doxygenclass:: ParticleSimulation::CylinderStartZone
    :members:
    :undoc-members:    

Utilities
=========

`PSim_util.hpp / .cpp` bundles a set of utility functions in the :cpp:any:`ParticleSimulation::util` namespace: 

.. doxygennamespace:: ParticleSimulation::util
   :undoc-members:

`PSim_math.hpp / .cpp` bundles some math functions: 

.. doxygenfile:: PSim_math.hpp