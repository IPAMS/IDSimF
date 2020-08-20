.. _usersguide-general_structure:

===========================
General structure of IDSimF
===========================


IDSimF is a framework implemented in c++ with a structure similar to `OpenFOAM <https://openfoam.org/>`_. The core of IDSimF is a set of `modules` which provide the functionality to simulate ion trajectories, e.g. integration algorithms, data structures which represent particles or collision models which describe the interaction between simulated particles and a background gas. 

Programs which solve a concrete simulation problem is an `application` in IDSimF terms. Every application is itself a stand alone c++ program which uses the IDSimF modules for the simulation task. The applications bundled with IDSimF are described in the :doc:`applications user guide <applications>`.

 
File and directory structure of the IDSimF repository
=====================================================

In addition to `modules` and `applications`, the IDSimF repository contains further directories, which contain for example integrated libraries, tests or helping scripts written in different programming languages. 

The structure of the IDSimF repository is:


* ``applications``: Simulation applications, solvers implementing concrete simulations. The applications are described in the :doc:`applications users guide <applications>`.
* ``modules``: Framework modules of IDSimF. They are described in the :doc:`module documentation <../modules/module_overview>`.
* ``libs``: Libraries from other open source projects which are integrated in IDSimF. 
* ``tests``: Unit Tests and Benchmarks.
* ``julia``: Modules for the `Julia Programming Language <https://julialang.org/>`_, primarily a monte carlo simulator to generate statistical input for :doc:`collision models <../modules/collisionmodels>`
* ``python``: Some helping / support scripts written in `Python <https://www.python.org/>`_, primarily generator scripts to generate some of the test input used in the unit tests.
* ``documentation``: The source files for this documentation
* ``notebooks``: `Jupyter <https://jupyter.org/>`_ notebooks with additional tests and test result analysis
