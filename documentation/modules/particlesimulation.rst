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

