.. _modules-integration:

===================
Module: Integration
===================

The Integration module bundles numerical integrators for equations of motion / particle motion.

Abstract Interface classes
==========================

.. doxygenclass:: Integration::AbstractTimeIntegrator
    :members:
    :undoc-members:


Simple Integrators
==================

.. doxygenclass:: Integration::VelocityIntegrator
    :members:
    :undoc-members:

Velocity Verlet (Leapfrog type) Integrators
===========================================

Velocity verlet integration is a leapfrog type integration algorithm. IDSimF has different velocity verlet integrators with different space charge / particle-particle interaction calculation schemes 
(Full sum, serial / parallel BTree, FMM with different FMM libraries)


.. doxygenclass:: Integration::VerletIntegrator
    :members:
    :undoc-members:

.. doxygenclass:: Integration::ParallelVerletIntegrator
    :members:
    :undoc-members:

.. doxygenclass:: Integration::FMMVerletIntegrator
    :members:
    :undoc-members:

.. doxygenclass:: Integration::FullSumVerletIntegrator
    :members:
    :undoc-members:


Runge-Kutta 4 Integrators
=========================

Integrators with the classical Runge-Kutta 4 algorithm with different space charge / particle-particle interaction calculation schemes (Full sum, parallel BTree)


.. doxygenclass:: Integration::ParallelRK4Integrator
    :members:
    :undoc-members:

.. doxygenclass:: Integration::FullSumRK4Integrator
    :members:
    :undoc-members:
