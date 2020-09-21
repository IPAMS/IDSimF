.. _testing-installation:

=======================
Testing an IDSimF build
=======================

IDSimF can be tested after an out of source compilation by running unit tests, running benchmarks or running simple test simulations in the build tree. 


Running unit tests
==================

IDSimF uses the `Catch2 <https://github.com/catchorg/Catch2>`_ test framework to verify the correct function of its modules with unit tests. The unit tests for the individual modules of IDSimF are located in ``tests/modules``. After compilation, every subfolder in ``tests/modules`` contains an test runner executable binary named ``test_<modulename>``. To run the tests of a module, execute the test runner binary in its folder. For example, run the tests for the `core` module with 

.. code-block:: console
    
    cd tests/modules/core
    ./test_core

this should produce an output similar to 

.. code-block:: none

    ===============================================================================
    All tests passed (4035 assertions in 4 test cases)

All module unit tests should pass. Note that some unit tests are statistical, use actual random numbers and thus fail with a slight probability. In case of tests failing, rerun the test runner a couple of times to check if they are consistently failing which would be a defect of IDSimF. 


Running benchmarks
==================

IDSimF has a couple of benchmark programs which also test aspects of IDSimF but have a longer runtime than the usual unit tests. The benchmarks are located in ``tests/benchmarks``. 
Most of the executables do not expect command line options and can be run without further arguments. 

For example the benchmark in ``tests/benchmarks/interpolatedField/runtime_benchmark`` performs a test of sampling in interpolated vector and scalar fields. It should produce an output similar to 

.. code-block:: console

    $ ./benchmark_interpolated_field_runtime
    Own implementation: sum:1.4e+08 vec:1.4e+08 2e+08 4e+07 elapsed:3.0186 for 40000000 samples

Detailed information about the available benchmarks can be found in the :doc:`benchmarks section of the users guide <../usersguide/benchmarks>`


Running simulations with example input 
======================================

The actual simulation applications / solvers can also be used with example input to test IDSimF. The simulation applications are located in ``applications``. They are grouped by topic into subdirectories. 

Simulation applications expect a simulation run configuration file as first argument and the simulation run name which should be used for writing the output files as second argument. Every application directory has an ``example`` directory with at least one example simulation run configuration for the respective simulation application. 

For example, ``applications/basic/BT-spaceChargeMinimalParallelSim`` implements a simple parallelized simulation of the electric particle-particle interaction (space charge) in a sphere of ions. To run it with an example configuration run change into the simulation application directory 

.. code-block:: console

    cd applications/basic/BT-spaceChargeMinimalParallelSim

and run a test simulation with 

.. code-block:: console

    ./BT-spaceChargeMinimalParallelSim example/benchmarkRun_short.json testRun

A test simulation should calculate 1500 time steps and a result file ``testRun_trajectories.hd5`` should be written. The result file contains the simulated ion trajectories. It can be analyzed with `IDSimPy <https://github.com/IPAMS/IDSimPy>`_, a Python based data analysis and preprocessing package. 
