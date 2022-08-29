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

Verlet Integrators with BTree
=============================

.. doxygenclass:: Integration::VerletIntegrator
    :members:
    :undoc-members:

.. doxygenclass:: Integration::ParallelVerletIntegrator
    :members:
    :undoc-members:


Verlet Integrators with FMM
===========================

.. doxygenclass:: Integration::FMMVerletIntegrator
    :members:
    :undoc-members: