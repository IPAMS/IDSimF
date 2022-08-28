.. _modules-core:

============
Module: Core
============

The Core module bundles basic data structures, particularly the basic spatial vector class, and basic facilities required in most other modules. 

Basic 3D Vector: :cpp:class:`Core::Vector`
==========================================

:cpp:class:`Core::Vector` is the fundamental spatial vector class used in IDSimF: 

.. doxygenclass:: Core::Vector
    :members:
    :undoc-members:


Particle Class
==============

IDSimF implements particle based simulation methods. :cpp:class:`Core::Particle` is the current base class for all simulated particles in IDSimF. It bundles the basic characteristics of simulated charged particles in IDSimF.

.. doxygenclass:: Core::Particle
    :members:
    :undoc-members:



Random Generators / Random Distributions
========================================

--------------
Random sources
--------------

The Core module provides basic randomness / random bit sources: 

.. doxygenclass:: Core::RandomSource
    :members:
    :undoc-members:

.. doxygenclass:: Core::MersenneBitSource
    :members:
    :undoc-members:

.. doxygenclass:: Core::TestBitSource
    :members:
    :undoc-members:


------------
Random pools
------------

Randomness generators are organized in pools for multithreaded applications: 

.. doxygenclass:: Core::AbstractRandomGeneratorPool
    :members:
    :undoc-members:

.. doxygenclass:: Core::RandomGeneratorPool
    :members:
    :undoc-members:

.. doxygenclass:: Core::TestRandomGeneratorPool
    :members:
    :undoc-members:

--------------------
Random distributions
--------------------

The random generators can return random distributions, which generates random variables with a defined statistical distribution: 

.. doxygenclass:: Core::RandomDistribution
    :members:
    :undoc-members:

.. doxygenclass:: Core::UniformRandomDistribution
    :members:
    :undoc-members:

--------------------
Test distributions
--------------------

Similarly to the test random generator, there are test random distributions which generate a small, predictable sequence of pre calculated values for testing purposes: 

.. doxygenclass:: Core::UniformTestDistribution
    :members:
    :undoc-members:

.. doxygenclass:: Core::NormalTestDistribution
    :members:
    :undoc-members:



Physical Constants
==================

`Core_constants.hpp` defines a set of physical constants widely used across IDSimF:  

.. doxygenfile:: Core_constants.hpp


Math Functions
==============

`Core_math.hpp` defines a set of general mathematical functions used elsewhere in the framework

.. doxygenfunction:: Core::degToRad
.. doxygenfunction:: Core::radToDeg
.. doxygenfunction:: Core::cartesianToPolar
.. doxygenfunction:: elevationRotate
.. doxygenfunction:: azimuthRotate
