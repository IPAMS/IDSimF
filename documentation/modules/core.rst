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

Random Generators / Random Distributions
========================================

The Core module provides production and test random generators: 

.. doxygenclass:: Core::AbstractRandomGenerator
    :members:
    :undoc-members:

.. doxygenclass:: Core::RandomGenerator
    :members:
    :undoc-members:

.. doxygenclass:: Core::TestRandomGenerator
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

.. doxygenclass:: Core::NormalRandomDistribution
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