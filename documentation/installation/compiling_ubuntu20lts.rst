.. _installation-compilation-ubuntu20lts:

---------------------------------------------------
Building guide: Linux - Ubuntu 20 LTS (Focal Fossa)
---------------------------------------------------

IDSimF can be compiled on Ubuntu 20 LTS (Focal Fossa) with comparably small effort. The configuration and compilation process should be very similar on other Debian derived Linux distributions. Configuration and compilation on other distribution families should also be similar, if the required dependencies are available on the target platform. 


Update package sources, check / install ``git`` and ``cmake``: 
--------------------------------------------------------------


First the package sources should be updated: 

.. code-block::

    sudo apt update

Usually git should be installed already, however check if ``git`` is really installed: 

.. code-block::

    git --version 

prints the installed git version. If ``git`` is not installed, install it with 

.. code-block::

    sudo apt install git 

``cmake`` is usually not installed. Install it with 

.. code-block::

    sudo apt install cmake


Install GCC C++ compiler:
-------------------------

The c++ compiler of `gnu compiler collection (gcc) <https://gcc.gnu.org/>`_  gnu compiler collection (``g++``
) in major version 9 is readily available on Ubuntu 20 LTS and is fully compatible with IDSimF. Install it with 

.. code-block::

    sudo apt install g++

Clone the IDSimF repository 
---------------------------

.. include:: default_git_clone.rst


Configuration and building with cmake
--------------------------------------


Preparing the build
...................

.. include:: default_cmake_configuration.rst

    
Building
........

.. include:: default_build.rst


Test the build
..............

After compilation has finished without problems, the IDSimF build can be tested by :doc:`running tests or benchmarks <testing_installation>`.
