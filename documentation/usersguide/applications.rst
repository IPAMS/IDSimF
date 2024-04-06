.. _usersguide-applications:

====================
IDSimF: Applications
====================

Applications in IDSimF are programs which use the IDSimF modules to solve an actual simulation problem. They are located in the ``applications`` directory and are grouped into subfolders by topic. The compiled applications are stand alone programs, since the IDSimF modules are typically linked statically into the application executables. 

Running applications 
====================


Applications are typically called with an configuration file defining the parameters of a simulation run as first argument and a simulation name which is used as base name for the written result files as second argument. 

For example, running the simulation application ``BT-spaceChargeMinimalParallelSim``, which implements a simple ion cloud with space charge interaction, using a simulation configuration ``simConfiguration.json`` with 

.. code-block:: console

    $ ./BT-spaceChargeMinimalParallelSim simConfiguration.json simulationRun001

creates a result file ``simulationRun001_trajectories.hd5`` and a simulation log file ``simulationRun001.log``. The trajectory file contains the simulated ion trajectory data in a structured `HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ file which can be read and analyzed with `IDSimPy <https://github.com/IPAMS/IDSimPy>`_. The log file contains the status information presented during the simulation to the user, which is also printed to the terminal. 

Which types of result files are produced depends on the simulation type, but all simulation applications generate a simulation log file. 

-------------------------------
Optional command line arguments
-------------------------------

* The number of threads in multithreaded applications can be controlled with the ``--n_threads <number of threads>`` (alias to ``-n <number of threads>``) option.  
* A help message with the command line arguments for a simulation application is printed with the ``--help`` switch. 


Simulation run configurations
=============================

The simulation configurations are provided as `JSON <https://www.json.org/>`_ [#fn_json_additional]_ formatted text files. They contain a set of key-value pairs in a dictionary structure, which set the values of the parameters in the respective solver applications. The values of the simulation parameters can be numeric, character strings like file paths or keywords, or arrays of values. Since simulation configuration files are plain text files and human readable, simulation parameters can be changed by editing the simulation configuration file with a text editor, e.g. `Visual Studio Code <https://code.visualstudio.com/>`_. 

Strictly, the JSON standard allows no comments. However, IDSimF simulation configuration files can contain comments analog to the comments in the c-family of programming languages (c, c++, javascript etc.). Therefore, 

.. code-block:: javascript

    {
        /* Block comment
        Multiline comments also can be set as a block comment */
        "sim_time_steps":1500,   //inline comment
        "n_ions":[300,300],             
        "ion_masses":[35,37]            
    }

is a valid IDSimF simulation configuration input. 

If IDSimF configuration files should be parsed with an external program, the comments has to stripped from the file to get standard compliant JSON. One option for this is to preprocess the configuration files content with a minify tool for JSON, e.g. the Python package `JSON  minify <https://github.com/getify/JSON.minify>`_. 


Simulation Applications Documentation
=====================================

This section gives an overview of the simulation application categories and the individual simulation applications in IDSimF. Every application has a subsection where the simulation run configuration parameters are described in detail.

--------------------------------------------
``basic``: Basic / Idealized Ion interaction 
--------------------------------------------

:doc:`spaceChargeMinimalSim <applications/spaceChargeMinimalSim>`: Defined particle cloud interacting by space charge
    Minimal simulation of an ensemble of charged particles interacting by space charge. Particle ensemble is predefined by an input file. Space charge is calculated with a Barnes-Hut tree. 
:doc:`spaceChargeSimpleSim <applications/spaceChargeSimpleSim>`: Random particle particles interacting by space charge
    Minimal simulation of a random ensemble of charged particles interacting by space charge. Particle ensemble initialized as a random box. Space charge is calculated with a serial version of a Barnes-Hut tree. 
:doc:`staticSimionPASim <applications/staticSimionPASim>`: Particle cloud interacting by space charge in an arbitrary electrode geometry
    A random ensemble of charged particles in an arbitrary electric potential distribution / electrode geometry defined by SIMION potential arrays, interacting by electric particle-particle force (space charge). 

.. toctree::
    :maxdepth: 1
    :hidden:

    applications/spaceChargeMinimalSim
    applications/spaceChargeSimpleSim
    applications/staticSimionPASim

-------------------------------------
``ionMobility``: Ion mobility devices
-------------------------------------

:doc:`IMSSim <applications/IMSSim>`: Simulation of drift tube Ion Mobility Spectrometry
    Simulation of a drift tube Ion Mobility Spectrometry Device, including background gas interaction, ion chemistry and space charge. 

:doc:`DMSSim <applications/DMSSim>` Simulation of planar Differential Ion Mobility 
    Simulation of a Differential Ion Mobility (DMS) separation device with idealized planar electrodes, including background gas interaction, ion chemistry and space charge. 

:doc:`DMSSimplifiedSim <applications/DMSSimplifiedSim>`: Simplified / idealized simulation of planar Differential Ion Mobility
    Simplified simulation of a planar electrode Differential Ion Mobility (DMS) separation device: Search for simulated compensation voltage (CV) of a chemically active ion ensemble in DMS.
    
:doc:`TWIMSSim <applications/TWIMSSim>`: Simulation of Traveling Wave Ion Mobility Spectrometry
    Simulation of a Traveling Wave Ion Mobility Spectrometry Device, including background gas interaction, ion chemistry and space charge.  

.. toctree::
    :maxdepth: 1
    :hidden:

    applications/IMSSim
    applications/DMSSim
    applications/DMSSimplifiedSim
    applications/TWIMSSim


------------------------------------------------------
``ionTransfer``: Ion transfer devices / ion guides 
------------------------------------------------------

:doc:`generalQuadSim <applications/generalQuadSim>`: Transfer (RF-only) quadrupole with arbitrary electrode geometry
    Ion trajectories in an RF only quadrupole device with arbitrary electrode geometry considering space charge and background gas collisions. The electrode geometry is defined by SIMION potential arrays.

.. toctree::
    :maxdepth: 1
    :hidden:

    applications/generalQuadSim

--------------------------------------------------------
``ionCollision``: Ion activation / ion collision devices
--------------------------------------------------------

:doc:`quadrupoleCollisionCellSim <applications/quadrupoleCollisionCellSim>`: Quadrupolar collsion cell 
    A quadrupolar collision cell with hard sphere collisions between ions and background gas, space charge and variable electrode geometry given by SIMION potential arrays. 

.. toctree::
    :maxdepth: 1
    :hidden:

    applications/quadrupoleCollisionCellSim

------------------------------
``ionTraps``: Ion trap devices
------------------------------

:doc:`QITSim <applications/QITSim>`: Idealized Quadrupole Ion Trap (QIT) with space charge and background gas collisions
    Ion trajectory simulation of an idealized quadrupole ion trap (QIT) considering space charge effects and hard sphere collisions with neutral background gas particles. The electric field of the ion trap is defined by idealized, analytical equations and FFT detection by determining the mirror charge on the cap electrodes is supported.

:doc:`reactiveQITSim <applications/reactiveQITSim>`: Idealized Quadrupole Ion Trap (QIT) with space charge, background gas collisions and reactive ions
    Ion trajectory simulation of an idealized quadrupole ion trap (QIT) considering space charge effects, neutral background particle collisions and chemical reactions of the trapped ions. The electric field of the ion trap is defined by idealized, analytical equations. 

:doc:`generalTrapSim <applications/generalTrapSim>`: RF ion trap with arbitrary ion geometry, space charge and background gas collisions
    Simulation of ion trajectories in an RF ion trap device with arbitrary geometry considering space charge and hard sphere collisions with neutral background gas particles. The electrode geometry is defined by SIMION potential arrays.

.. toctree::
    :maxdepth: 1
    :hidden:

    applications/QITSim
    applications/reactiveQITSim
    applications/generalTrapSim

--------------------------------
``chemistry``: Chemical Kinetics
--------------------------------

:doc:`idealIsothermReactorSim <applications/idealIsothermReactorSim>`: Chemical kinetics in ideally mixed isotherm chemical reactor 
    Particle based kinetics simulation of the reaction of an ensemble of reactive particles with reaction partners in a background gas in an ideally mixed isotherm reactor. 


.. toctree::
    :maxdepth: 1
    :hidden:

    applications/idealIsothermReactorSim


Guide: How to build an IDSimF Application?
==========================================

.. toctree::
    :maxdepth: 1
    
    application_development

.. rubric:: Footnotes

.. [#fn_json_additional] See also `the JSON wikipedia page <https://en.wikipedia.org/wiki/JSON>`_ for further information. 
