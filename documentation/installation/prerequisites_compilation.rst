.. _installation-compilation:

=============================
Prerequisites and Compilation
=============================

IDSimF is a pure command line tool written in C++ which requires a POSIX compatible environment. It is primarily developed and used on Linux and MacOS systems. However, compilation is possible on recent Windows systems with relatively moderate effort within the Linux Subsystem for Windows.

Prerequisites
=============

---
Git
---

`Git <https://git-scm.com/>`_ is used as version control system for IDSimF. Thus it is also used to clone the `IDsimF source code repository <https://github.com/IPAMS/IDSimF>`_ from `Github <https://github.com/>`_ to a local system where IDSimF should be installed. 

-----
CMake
-----

IDSimF uses `CMake <https://cmake.org/>`_ as build tool. IDSimF requires at least CMake version 3.9. 

--------
Compiler
--------

IDSimF is developed in C++17 and uses `OpenMP <https://www.openmp.org/>`_, thus a recent compiler with OpenMP support is required. 

Compilation of the current IDSimF release was tested with 

* clang version 14 (installed with Macports) on MacOS 12 (Monterey)
* gcc version 12 (installed with Macports) on  MacOS 12 (Monterey)
* gcc version 9.4 on Ubuntu 20 LTS
* gcc version 8.2.1 on openSUSE Leap 42.3 (with ``USE_CPP_FSLIB`` flag on and activated ``Exafmm-t`` and ``FMM3D`` external FMM libraries)
* gcc version 9.3.0. on Windows 10 in an Ubuntu 20.04 LTS within the Windows Subsystem for Linux

IDSimF was tested on systems with x86-64 and ARM (Apple M1) CPUs, currently there are no obvious restrictions in IDSimF regarding the CPU architecture.

Libraries and Dependencies
==========================

------------------
Required Libraries
------------------


IDSimF currently depends on one external libraries which have to be available on the target system: 

* `HDF5 library <https://www.hdfgroup.org/solutions/hdf5>`_ which provides fast access to HDF5 files

When the HDF5 library is installed on the target system, usually through a package management tool of the target system, CMake should find it without further configuration. 

There are currently two required libraries which are packaged into the source distribution of IDSimF: 

* `Catch2 <https://github.com/catchorg/Catch2>`_ unit test framework
* `jsoncpp <https://github.com/open-source-parsers/jsoncpp>`_ JSON file parsing library


----------------------------------------------
Optional Fast Multipole Method (FMM) Libraries
----------------------------------------------

The `Fast Multipole Method (FMM) <https://en.wikipedia.org/wiki/Fast_multipole_method>`_ is an approach to calculate the forces and therefore interactions between large numbers of particles in an an efficient way. It is used for many applications, e.g. gravitational interactions or coulombic interactions, therefore there are multiple open source libraries which provide FMM capabilities. 

Beside the internal fast coulomb solver, based on a `Barnes Hut Tree algorithm <https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation>`_, IDSimF can currently use two FMM libraries for the calculation of coulombic inter-particle interactions:

* `exafmm-t <https://exafmm.github.io/exafmm-t>`_ 
* `FMM3D <https://fmm3d.readthedocs.io/en/latest/index.html>`_ 

To prepare a build with activated FMM libraries, the locations of the libraries have to be provided to CMake in the configuration step, via the ``EXAFMMT_PATH`` and ``FMM_3D_PATH`` variables. Details are described in the compilation / installation guides below.


Guides: Cloning the IDSimF repository and compiling IDSimF 
==========================================================


The basic compilation / installation process has the same basic steps on all platforms: Checkout of the source files, preparing the system / installation of dependencies and compilation. While the compilation process is basically the same on all platforms, the preparation and setup process for compilation differs on the operating systems. 

.. note::
    The compilation / installation guides assume that you have familiarity with usage of the shell / command line on your target system and you have administrator rights on the target system. 

IDSimF can be built on all major operating systems. The following guides show in detail how to clone the IDSimF repository, configure the build process with cmake and build IDSimF. 


.. rubric:: Linux

.. toctree::
    :maxdepth: 1

    compiling_ubuntu18lts
    compiling_ubuntu20lts


.. rubric:: Windows

.. toctree::
    :maxdepth: 1

    compiling_windows


.. rubric:: MacOS

.. toctree::
    :maxdepth: 1

    compiling_macos
