After configuration, the individual built targets (IDSimF modules, simulation applications and tests) can be built with ``cmake --build <path to a target>``. Since the root of the build tree is also a target for the whole project, all build targets in IDSimF are built serially (with no parallelization in the build process) with

.. code-block::
    
    cmake --build .
    
``cmake`` supports parallelized builds which use multiple cpu cores. To perform a parallelized build, pass the ``-j`` option with the number of parallel jobs to the cmake build command: 

.. code-block:: shell

    cmake --build . -j <number of parallel jobs>

For example, to build all IDSimF targets with 8 cores use

.. code-block:: shell

    cmake --build . -j 8