.. _modules-btree:

=============
Module: BTree
=============

The BTree module provides Barnes-Hut tree implementations which are used to efficiently simulate the electric interaction between particles (space charge). 

Particle Class
==============

IDSimF implements particle based simulation methods. :cpp:class:`BTree::Particle` is the current base class for all simulated particles in IDSimF. It bundles the basic characteristics of simulated charged particles in IDSimF. 

.. doxygenclass:: BTree::Particle
    :members:
    :undoc-members:

Tree Node Base Classes 
======================

Barnes-Hut trees are in fact tree structures which are formed by a set of nodes and connections (edges) between them. :cpp:class:`BTree::AbstractNode` is an abstract base class for all tree nodes, it bundles all characteristics and methods all nodes have. :cpp:class:`BTree::GenericBaseNode` bundles generic basic methods which are applicable for all concrete node types, e.g. particle insertion into the node or particle removal from the node. 

.. doxygenclass:: BTree::AbstractNode
    :members:
    :undoc-members:

.. doxygenclass:: BTree::GenericBaseNode
    :members:
    :undoc-members:


Parallelized Barnes-Hut Tree
============================

A partly parallelized version of a Barnes-Hut tree is realized in :cpp:class:`BTree::ParallelTree`, which uses :cpp:class:`BTree::ParallelNode` as node class. 

.. doxygenclass:: BTree::ParallelTree
    :members:
    :undoc-members:


.. doxygenclass:: BTree::ParallelNode
    :members:
    :undoc-members:

Serial Barnes-Hut Tree
======================

A legacy, non parallelized, simple version of a Barnes-Hut tree is :cpp:class:`BTree::Tree`, which uses :cpp:class:`BTree::Node` as node class. This tree version is mostly used for benchmarks or comparison purposes. 

.. doxygenclass:: BTree::Tree
    :members:
    :undoc-members:


.. doxygenclass:: BTree::Node
    :members:
    :undoc-members: