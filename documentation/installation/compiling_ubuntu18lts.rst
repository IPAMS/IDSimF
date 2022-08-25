.. _installation-compilation-ubuntu18lts:

-----------------------------------------------------
Building guide: Linux - Ubuntu 18 LTS (Bionic Beaver)
-----------------------------------------------------

IDSimF was tested on Ubuntu 18 LTS extensively. Therefore this guide uses Ubuntu 18 LTS as example how to compile IDSimF on a Linux system. The configuration and compilation process should be very similar on other Debian derived Linux distributions. Configuration and compilation on other distribution families should also be similar, if the required dependencies are available on the target platform. 


Update package sources, Check / install Git and cmake: 
------------------------------------------------------

First the package sources should be updated: 

.. code-block:: console

    sudo apt update

Usually ``git`` should be installed already, however check if git is really installed: 

.. code-block:: console

    git --version 

prints the installed ``git`` version. If ``git`` is not installed, install it with 

.. code-block:: console

    sudo apt install git 

``cmake`` is usually not installed. Install it with 

.. code-block:: console

    sudo apt install cmake


Install a modern C++ compiler:
------------------------------

The default C++ compiler in Ubuntu is too old to support all used C++17 features. Therefore, a modern / recent compiler has to be installed. The c++ compiler of `gnu compiler collection (gcc) <https://gcc.gnu.org/>`_  gnu compiler collection (``g++``
) in major version 8 is readily available on Ubuntu 18 and is compatible with IDSimF. Install it with 

.. code-block:: console

    sudo apt install g++-8

Alternatively the `llvm <https://www.llvm.org/>`_ c++ frontend `clang 10 / llvm <https://clang.llvm.org/>`_ is also available and can be used to compile IDSimF. Install it with 

.. code-block:: console

    sudo apt install clang-10

Note that you have to install *both* packages (g++-8 and clang-10) to use ``clang``, since the clang-10 package is missing some header files. 

Clone the IDSimF repository 
---------------------------

Clone the IDSimF repository to your local machine with ``git`` from GitHub: 

.. code-block:: console
    
    git clone https://github.com/IPAMS/IDSimF.git

This clones the IDSimF repository to a local folder ``IDSimF``. 


Configuration and building with cmake
--------------------------------------


Preparing the build
...................

IDSimF uses ``cmake`` as helping tool for configuration and compilation. ``cmake`` allows a so called "out of source build". This creates a separated "binary tree" in a separated build folder, in which the compilation of executable binaries takes place without interfering with the cloned sources. 

To do an out of source build, change into the cloned IDSimF folder and create a build folder, for example ``build`` in it and change into it: 

.. code-block:: console
    
    cd IDSimF
    mkdir build
    cd build

Basically ``cmake`` prepares a build tree in the current folder if it is called with an source folder as argument. However, since we do not want to build with the default compiler of the system, we have to set some options for ``cmake``. This is done with optional arguments of the form ``-D<OPTION NAME>=<VALUE>``.

In the build folder, prepare build / binary tree with 

.. code-block:: console

    cmake .. -DCMAKE_CXX_COMPILER=/usr/bin/g++-8 -DCMAKE_BUILD_TYPE=Release -DUSE_CPP_FSLIB=on

The build options mean the following: 

    + ``-DCMAKE_CXX_COMPILER=/usr/bin/g++-8`` sets the used c++ compiler to the ``g++-8`` which is installed at ``/usr/bin/``. 
    + ``-DCMAKE_BUILD_TYPE=Release`` sets the build type to `Release` which means, that optimizations are switched on and debugging information is removed from the compiled binary. This results in a significantly faster binary than when building with debugging information switched on (``-DCMAKE_BUILD_TYPE=Debug``) which is required to analyze the compiled binary with debugging tools. 
    + ``-DUSE_CPP_FSLIB=on`` switches on explicitly linking against the filesystem library. This is explicitly required by `gcc` until gcc 9. With gcc 8, the build will fail without this linking. 

Configuring the build with ``clang`` is very similar: : 

.. code-block:: console

    cmake .. -DCMAKE_CXX_COMPILER=/usr/bin/clang++-10 -DCMAKE_BUILD_TYPE=Release

Here the c++ compiler installed with the `clang-10` package, ``clang++-10``, is used as compiler. Explicitly linking the filesystem library is not necessary with ``clang``. 

Optional FMM Libraries
......................

.. include:: default_exafmm_t.rst

.. include:: default_fmm3d.rst

Building
........

After configuration, the individual built targets (IDSimF modules, simulation applications and tests) can be built with ``cmake --build <path to a target>``. Since the root of the build tree is also a target for the whole project, all build targets in IDSimF are built serially (with no parallelization in the build process) with

.. code-block:: console
    
    cmake --build .

``cmake`` supports parallelized builds since version 3.12, but Ubuntu 18 LTS installs an older ``cmake`` version. However, the native build tool used by ``cmake`` to actually build is `make <https://en.wikipedia.org/wiki/Make_(software)>`_  which itself supports parallelized builds. ``cmake`` is able to pass options to the native build tool. A parallelized build with ``cmake`` and ``make`` is done with 

.. code-block:: console

    cmake --build . -- -j <number of parallel jobs>

Alternatively, ``make`` can also be used directly: 

.. code-block:: console

    make -j <number of parallel jobs>


Test the build
..............

After compilation has finished without problems, the IDSimF build can be tested by :doc:`running tests or benchmarks <testing_installation>`.
