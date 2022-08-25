.. _installation-compilation-windows:

-----------------------
Building guide: Windows
-----------------------

Building with Windows Subsystem for Linux and Ubuntu 20 LTS
-----------------------------------------------------------

Alternatively to Cygwin, Windows 10 and Windows Server 2019 are able to run Linux executables natively with an compatibility layer developed by Microsoft, the `Windows Subsystem for Linux (WSL) <https://docs.microsoft.com/en-us/windows/wsl/>`_. Different Linux distributions are available from the Microsoft store for installation into the WSL, including `Ubuntu 20.04 LTS <https://ubuntu.com/download/desktop>`_. After setup of WSL and installation of Ubuntu 20.04 LTS (Focal Fossa) within WSL, the build process of IDSimF is the same as in Ubuntu 20.04 LTS.

.. note::
    This guide assumes that you have basic familiarity with Linux / Unix shells and you can invoke commands, change directories, list directory contents etc. within a shell. 

Install Windows Subsystem for Linux (WSL) and Ubuntu 20.04 LTS
..............................................................

Install Windows Subsystem for Linux (WSL) according to the `installation guide <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_. Install Ubuntu 20.04 LTS from the `Microsoft Store <https://www.microsoft.com/en-us/p/ubuntu-2004-lts/9n6svws3rx71>`_. 

After the installation of Ubuntu 20.04 (Focal Fossa) in WSL, you should have an "Ubuntu" tile in your Windows menu, which starts a new shell within the Ubuntu installation. 


Build and test IDSimF within Ubuntu 20.04 LTS in WSL
....................................................

The build process of IDSimF in the Ubuntu system within WSL is exactly the same as in a natively installed Ubuntu system. Refer to the :doc:`Ubuntu 20.04 LTS installation guide <compiling_ubuntu20lts>` for the remaining build process.



Building with Cygwin 
--------------------

`Cygwin <https://cygwin.com/>`_  Cygwin is an UNIX (Posix compatible) programming and runtime environment. It provides basically a Linux distribution which runs natively under Windows. IDSimF can be built and run in the Cygwin environment. 

.. note::
    This guide assumes that you have basic familiarity with Linux / Unix shells and you can invoke commands, change directories, list directory contents etc. within a shell. 

Install Cygwin and required Cygwin packages
...........................................

Use the Cygwin setup program to install Cygwin and the following packages: 

+ ``git``
+ ``cmake``
+ ``g++```

+ ``gnu make``
+ ``hdf5``

This makes all required tools available in the Cygwin environment. You can check if ``git`` ``cmake`` and a c++ compiler is available in the Cygwin environment by checking the versions of the individual tools in a cygwin shell: 

.. code-block:: console

    $ git --version
    git version 2.28.0

    $ cmake --version
    cmake version 3.14.5

    $ g++ --version
    g++ (GCC) 9.3.0
    Copyright (C) 2019 Free Software Foundation, Inc.


Clone the IDSimF repository 
...........................

Within a Cygwin shell, clone the IDSimF repository to your local machine with git from GitHub: 

.. code-block:: console
    
    git clone https://github.com/IPAMS/IDSimF.git

This clones the IDSimF repository to a local folder "IDSimF". 


Configuration and building with cmake
.....................................

Preparing the build
+++++++++++++++++++

.. include:: default_cmake_configuration.rst

Optional FMM Libraries
++++++++++++++++++++++

The optional FMM libraries, `exafmm-t <https://exafmm.github.io/exafmm-t>`_ and `FMM3D <https://fmm3d.readthedocs.io/en/latest/index.html>`_ were not yet tested with cygwin. It should be possible to compile them within the cygwin environment if all dependencies of those libraries are available within cygwin.
    
Building
++++++++

.. include:: default_build.rst


Test the build
++++++++++++++

After compilation has finished without problems, the IDSimF build can be tested by :doc:`running tests or benchmarks <testing_installation>`.
