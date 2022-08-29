.. _modules-space_charge:

===================
Module: SpaceCharge
===================

The SpaceCharge module bundles different approaches to efficiently simulate the electric interaction between particles (space charge). 

The module contains submodules implementing individual calculation approaches, generalized interface classes are located in the module itself: 

.. doxygenclass:: SpaceCharge::FieldCalculator
    :members:
    :undoc-members:

.. doxygenclass:: SpaceCharge::GenericSpaceChargeSolver
    :members:
    :undoc-members:



Sub Module: BTree
=================

The BTree sub module provides an Barnes-Hut tree implementation for space charge calculation.

----------------------
Tree Node Base Classes 
----------------------

Barnes-Hut trees are in fact tree structures which are formed by a set of nodes and connections (edges) between them. :cpp:class:`BTree::AbstractNode` is an abstract base class for all tree nodes, it bundles all characteristics and methods all nodes have. :cpp:class:`BTree::GenericBaseNode` bundles generic basic methods which are applicable for all concrete node types, e.g. particle insertion into the node or particle removal from the node. 

.. doxygenclass:: BTree::AbstractNode
    :members:
    :undoc-members:

.. doxygenclass:: BTree::GenericBaseNode
    :members:
    :undoc-members:

----------------------------
Parallelized Barnes-Hut Tree
----------------------------

A partly parallelized version of a Barnes-Hut tree is realized in :cpp:class:`BTree::ParallelTree`, which uses :cpp:class:`BTree::ParallelNode` as node class. 

.. doxygenclass:: BTree::ParallelTree
    :members:
    :undoc-members:


.. doxygenclass:: BTree::ParallelNode
    :members:
    :undoc-members:

----------------------
Serial Barnes-Hut Tree
----------------------

A legacy, non parallelized, simple version of a Barnes-Hut tree is :cpp:class:`BTree::Tree`, which uses :cpp:class:`BTree::Node` as node class. This tree version is mostly used for benchmarks or comparison purposes. 

.. doxygenclass:: BTree::Tree
    :members:
    :undoc-members:


.. doxygenclass:: BTree::Node
    :members:
    :undoc-members:


Sub Modules: Fast Multipole Methods (FMM)
=========================================

Two submodules provide interfaces to external fast multipole method (FMM) libraries to allow fast space charge calculations with FMM. 

------------------
exafmm-t Interface 
------------------

Interface to `exafmm-t <https://exafmm.github.io/exafmm-t>`_ for coulombic interaction (space charge) calculations. This modules a build with the path to exafmm-t passed via ``EXAFMMT_PATH`` variable in the ``cmake`` configuration, as described in the :doc:`installation guide </installation/prerequisites_compilation>`.

.. doxygenclass:: ExaFMMt::FMMSolver
    :members:
    :undoc-members:


---------------
FMM3D Interface 
---------------

Interface to `FMM3D <https://fmm3d.readthedocs.io/en/latest/index.html>`_ for coulombic interaction (space charge) calculations. This modules a build with the path to exafmm-t passed via ``FMM_3D_PATH`` variable in the ``cmake`` configuration, as described in the :doc:`installation guide </installation/prerequisites_compilation>`

.. doxygenclass:: FMM3D::FMMSolver
    :members:
    :undoc-members:
