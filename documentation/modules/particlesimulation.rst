.. _modules-particlesimulation:

==========================
Module: particleSimulation
==========================

The `ParticleSimulation` module provides algorithms, data structures, input/output and data transformation methods required for actual charged particle trajectory simulations. 


Time Integrators
================

A central aspect of IDSimF is time integration of the particle motion, which is done by time integrators. 


:cpp:class:`ParticleSimulation::AbstractTimeIntegrator` is the abstract super class of all time integrators, which provides a basic common interface for time integration: 

.. doxygenclass:: ParticleSimulation::AbstractTimeIntegrator
    :members:
    :undoc-members:


-------------------
Velocity Integrator
-------------------

:cpp:class:`ParticleSimulation::VelocityIntegrator` is a simple integrator which uses the current velocity in a time step to move the particle. Thus, it ignores particle inertia. Furthermore, space charge and background gas interactions are ignored by this integrator. 

.. doxygenclass:: ParticleSimulation::VelocityIntegrator
    :members:
    :undoc-members:


------------------
Verlet Integrators
------------------

:cpp:class:`ParticleSimulation::VerletIntegrator` and :cpp:class:`ParticleSimulation::ParallelVerletIntegrator` are a serial (non parallelized) and parallelized version of a time integrator implementing the Verlet integration scheme. Space charge and background gas interactions are considered with verlet integrators. 

.. doxygenclass:: ParticleSimulation::VerletIntegrator
    :members:
    :undoc-members:


.. doxygenclass:: ParticleSimulation::ParallelVerletIntegrator
    :members:
    :undoc-members:


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


File Writers
============

File writer classes are used to export data to persistent files. 


------------------------------------------
Primary Simulation Result Data File Writer
------------------------------------------

:cpp:class:`ParticleSimulation::TrajectoryHDF5Writer` writes trajectory HDF5 files, which is the current primary trajectory data export format of IDSimF.

.. doxygenclass:: ParticleSimulation::TrajectoryHDF5Writer
    :members:
    :undoc-members:


:cpp:class:`ParticleSimulation::Scalar_writer` writes tables of scalar values from simulations: 

.. doxygenclass:: ParticleSimulation::Scalar_writer
    :members:
    :undoc-members:


-----------------------------
Additional Result File Writer
-----------------------------

Additional file writer provide additional export file formats. 

.. note:: 

    The additional file writer are currently not well maintained. 

.. doxygenclass:: ParticleSimulation::TrajectoryExplorerJSONwriter
    :members:
    :undoc-members:

.. doxygenclass:: ParticleSimulation::SimpleVTKwriter
    :members:
    :undoc-members:


------------------------------
Special Simulation File Writer
------------------------------

There are some file writers for special simulation requirements: 

.. doxygenclass:: ParticleSimulation::InductionCurrentWriter
    :members:
    :undoc-members:

.. doxygenclass:: ParticleSimulation::IdealizedQitFFTWriter
    :members:
    :undoc-members:

.. doxygenclass:: ParticleSimulation::AverageChargePositionWriter
    :members:
    :undoc-members:        


File Readers
============

File readers import data from persistent files 

:cpp:class:`ParticleSimulation::HDF5Reader` is a general reader for HDF5 files. 

.. doxygenclass:: ParticleSimulation::HDF5Reader
    :members:
    :undoc-members:

.. doxygenclass:: ParticleSimulation::IonCloudReader
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