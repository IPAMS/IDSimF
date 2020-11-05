.. _modules-rs:

================================
Module: RS - Reaction Simulation
================================

Particle based modeling of chemical kinetics and reaction dynamics of charged particles. 

Base Classes / Simulation System
================================

.. doxygenclass:: RS::Simulation
    :members:
    :undoc-members:

.. doxygenclass:: RS::SimulationConfiguration
    :members:
    :undoc-members:

.. doxygenclass:: RS::ReactiveParticle
    :members:
    :undoc-members:


Reactions
=========

RS provides different reaction types. Each reaction type is a subclass of :cpp:class:`RS::AbstractReaction`

.. doxygenclass:: RS::AbstractReaction
    :members:
    :undoc-members:

.. doxygenclass:: RS::StaticReaction
    :members:
    :undoc-members:

.. doxygenclass:: RS::StaticThermalizingReaction
    :members:
    :undoc-members:

.. doxygenclass:: RS::SimpleCollisionStepReaction
    :members:
    :undoc-members:

.. doxygenclass:: RS::VantHoffReaction
    :members:
    :undoc-members:

.. doxygenclass:: RS::FieldDependentVantHoffReaction
    :members:
    :undoc-members:


Data Import / Export
====================

.. doxygenclass:: RS::ConfigFileParser
    :members:
    :undoc-members:

.. doxygenclass:: RS::ConcentrationFileWriter
    :members:
    :undoc-members:


Utilities
=========

`RS_util.hpp / .cpp` bundles a set of utility functions in the :cpp:any:`RS::util` namespace: 

.. doxygennamespace:: RS::util
   :undoc-members:

`RS_constants.hpp` bundles a set of constants: 

.. doxygenfile:: RS_constants.hpp