.. _installation-compilation-macos:

---------------------
Building guide: MacOS
---------------------

IDSimF is developed and tested on MacOS extensively. This guide is describing compilation on MacOS Monterey with the open source package management system `Macports <https://www.macports.org/>`_. However, compilation of IDSimF on other MacOS versions and with different package management systems, particularly `Homebrew <https://brew.sh>`_, should be possible in a similar way to this guide. 

This installation guide assumes that you have `Macports <https://www.macports.org/>`_ installed. Xcode and the Xcode commandline tools are required for Macports and are therefore assumed to be installed too. This guide also assumes that you have basic familiarity with the command line / terminal and the ``port`` command line tool of Macports. 

Update package sources, install ``git`` and ``cmake``:
------------------------------------------------------

First the macports installation and the package sources should be updated: 

.. code-block:: console

    sudo port selfupdate

Install ``git`` with macports:

.. code-block:: console

    sudo port install git

Install ``cmake`` with macports:

.. code-block:: console

    sudo port install cmake

This should make ``git`` and ``cmake`` available. This can be verified by checking their versions: 

.. code-block:: console

    $ git --version
    git version 2.28.0

    $ cmake --version
    cmake version 3.23.3


Install a C++17 compatible C++ compiler :
-----------------------------------------

IDSimF is written in C++17, therefore a recent compiler fully supporting the C++17 language standard has to be used to compile IDSimF. The c++ compiler of `gnu compiler collection (gcc) <https://gcc.gnu.org/>`_  gnu compiler collection (``g++``
) in major version 12 is readily available with macports and is compatible with IDSimF. Install it with 

.. code-block:: console

    sudo port install gcc12

Alternatively, the `llvm <https://www.llvm.org/>`_ c++ frontend `clang <https://clang.llvm.org/>`_  in major version 14 is also available and can be used to compile IDSimF. Install it with 

.. code-block:: console

    sudo port install clang-14


Clone the IDSimF repository 
---------------------------

.. include:: default_git_clone.rst


Configuration and building with ``cmake``
---------------------------------------------


Preparing the build
...................

IDSimF uses ``cmake`` as helping tool for configuration and compilation. ``cmake`` allows a so called "out of source build". This creates a separated "binary tree" in a separated build folder, in which the compilation of executable binaries takes place without interfering with the cloned sources. 

To do an out of source build, change into the cloned IDSimF folder and create a build folder, for example ``build`` in it and change into it: 

.. code-block:: console
    
    cd IDSimF
    mkdir build
    cd build

Basically ``cmake`` prepares a build tree in the current folder if it is called with an source folder as argument. However, since we do not want to build with the default compiler of the system, we have to set some options for ``cmake``. This is done with optional arguments of the form ``-D<OPTION NAME>=<VALUE>``.

In the build folder, prepare a build / binary tree for the compilation with ``g++``
 with 

.. code-block:: console

    cmake .. -DCMAKE_CXX_COMPILER=g++-mp-12 -DCMAKE_BUILD_TYPE=Release

The build options are: 

    + ``-DCMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-12`` sets the used c++ compiler to the ``g++-12`` compiler with `OpenMP <https://www.openmp.org>`_ parallelization support. Macports installs the compiler at ``/opt/local/bin/`` but also sets the ``PATH`` variable, so that the path to the compiler has not to be set explicitly. 
    + ``-DCMAKE_BUILD_TYPE=Release`` sets the build type to `Release` which means, that optimizations are switched on and debugging information is removed from the compiled binary. This results in a significantly faster binary than when building with debugging information switched on (``-DCMAKE_BUILD_TYPE=Debug``) which is required to analyze the compiled binary with debugging tools. 

Configuring a build with ``clang`` is very similar:

.. code-block:: console

    cmake .. -DCMAKE_CXX_COMPILER=clang++-mp-14 -DCMAKE_BUILD_TYPE=Release

Here the c++ compiler with `OpenMP <https://www.openmp.org>`_ parallelization support, installed with the `clang-14` package, ``clang++-mp-14`` is used as compiler. 

Building
........

.. include:: default_build.rst

Test the build
..............

After compilation has finished without problems, the IDSimF build can be tested by :doc:`running tests or benchmarks <testing_installation>`.

