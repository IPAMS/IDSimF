.. _modules-collisionmodels:

=======================
Module: collisionModels
=======================

The collisionModels module provides models to describe the collisional interaction between charged particles and neutral background gas. 

Collision models
================

------------
Base Classes
------------

:cpp:class:`CollisionModel::AbstractCollisionModel` is the abstract base class for all collision models. 

.. doxygenclass:: CollisionModel::AbstractCollisionModel
    :members:
    :undoc-members:


:cpp:class:`CollisionModel::MultiCollisionModel` combines multiple collision models, primarily :cpp:class:`CollisionModel::HardSphereModel`, to model background gas mixtures in the hard sphere model. 

.. doxygenclass:: CollisionModel::MultiCollisionModel
    :members:
    :undoc-members:


---------------------------
Hard Sphere Collision Model
---------------------------

:cpp:class:`CollisionModel::HardSphereModel` implements a low pressure, collision model based on hard sphere collisions: 

.. doxygenclass:: CollisionModel::HardSphereModel
    :members:
    :undoc-members:


-------------------------------------
Statistical Diffusion Collision Model
-------------------------------------

:cpp:class:`CollisionModel::StatisticalDiffusionModel` implements a high pressure statistical collision model: 

.. doxygenclass:: CollisionModel::StatisticalDiffusionModel
    :members:
    :undoc-members:

:cpp:class:`CollisionModel::CollisionStatistics` provides normalized collision statistics / density functions for statistical diffusion collision models:

.. doxygenclass:: CollisionModel::CollisionStatistics
    :members:
    :undoc-members:

Utilities
=========

`CollisionModel_util.hpp / .cpp` bundles a set of utility functions in the :cpp:any:`CollisionModel::util` namespace: 

.. doxygennamespace:: CollisionModel::util
   :undoc-members:


`CollisionModel_MathFunctions.hpp / .cpp` bundles some math functions: 

.. doxygenfile:: CollisionModel_MathFunctions.hpp