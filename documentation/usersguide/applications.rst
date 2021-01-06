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

creates a result file ``simulationRun001_trajectories.hd5``. This result file contains the simulated ion trajectory data in a structured `HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ file which can be read and analyzed with `IDSimPy <https://github.com/IPAMS/IDSimPy>`_.

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

:doc:`BT-spaceChargeMinimalSim <applications/BT_spaceChargeMinimalSim>`: Particle cloud interacting by space charge (serial version)
    Minimal serial (non parallelized) simulation of an ensemble of charged particles interacting by space charge. Particle ensemble is predefined by an input file. Space charge is calculated with a Barnes-Hut tree. 
:doc:`BT-spaceChargeMinimalParallelSim <applications/BT_spaceChargeMinimalParallelSim>`: Particles interacting by space charge (parallel version)
    Minimal parallelized simulation of an ensemble of charged particles interacting by space charge. Particle ensemble is predefined by an input file or a random box of particles is used as start condition. Space charge is calculated with a parallelized Barnes-Hut tree. 
:doc:`BT-spaceChargeSimpleSim <applications/BT_spaceChargeSimpleSim>`: Random particle particles interacting by space charge (serial version)
    Minimal serial (non parallelized) simulation of a random ensemble of charged particles interacting by space charge. Particle ensemble initialized as a random box. Space charge is calculated with a serial version of a Barnes-Hut tree. 
:doc:`BT-staticSimionPASim <applications/BT_staticSimionPASim>`: Particle cloud interacting by space charge in an arbitrary electrode geometry
    A random ensemble of charged particles in an arbitrary electric potential distribution / electrode geometry defined by SIMION potential arrays, interacting by electric particle-particle force (space charge). 

.. toctree::
    :maxdepth: 1
    :hidden:

    applications/BT_spaceChargeMinimalSim
    applications/BT_spaceChargeMinimalParallelSim
    applications/BT_spaceChargeSimpleSim
    applications/BT_staticSimionPASim

-------------------------------------
``ionMobility``: Ion mobility devices
-------------------------------------

:doc:`BT-RS-IMSSim <applications/BT_RS_IMSSim>`: Simulation of drift tube Ion Mobility Spectrometry
    Simulation of a drift tube Ion Mobility Spectrometry Device, including background gas interaction, ion chemistry and space charge. 

:doc:`BT-RS-DMSSim <applications/BT_RS_DMSSim>` Simulation of planar Differential Ion Mobility 
    Simulation of a Differential Ion Mobility (DMS) separation device with idealized planar electrodes, including background gas interaction, ion chemistry and space charge. 

:doc:`RS-DMSSimplifiedSim <applications/BT_DMSSimplifiedSim>`: Simplified / idealized simulation of planar Differential Ion Mobility
    Simplified simulation of a planar electrode Differential Ion Mobility (DMS) separation device: Search for simulated compensation voltage (CV) of a chemically active ion ensemble in DMS. 

.. toctree::
    :maxdepth: 1
    :hidden:

    applications/BT_RS_IMSSim
    applications/BT_RS_DMSSim
    applications/RS_DMSSimplifiedSim


------------------------------------------------------
``ionTransfer``: Ion transfer devices / ion guides 
------------------------------------------------------

:doc:`BT-generalQuadSim <applications/BT_generalQuadSim>`: Transfer (RF-only) quadrupole with arbitrary electrode geometry
    Ion trajectories in an RF only quadrupole device with arbitrary electrode geometry considering space charge and background gas collisions. The electrode geometry is defined by SIMION potential arrays.

.. toctree::
    :maxdepth: 1
    :hidden:

    applications/BT_generalQuadSim

--------------------------------------------------------
``ionCollision``: Ion activation / ion collision devices
--------------------------------------------------------

:doc:`BT-quadrupoleCollisionCellSim <applications/BT_quadrupoleCollisionCellSim>`: Quadrupolar collsion cell 
    A quadrupolar collision cell with hard sphere collisions between ions and background gas, space charge and variable electrode geometry given by SIMION potential arrays. 

.. toctree::
    :maxdepth: 1
    :hidden:

    applications/BT_quadrupoleCollisionCellSim

------------------------------
``ionTraps``: Ion trap devices
------------------------------

:doc:`BT-QITSim <applications/BT_QITSim>`: Idealized Quadrupole Ion Trap (QIT) with space charge and background gas collisions
    Ion trajectory simulation of an idealized quadrupole ion trap (QIT) considering space charge effects and hard sphere collisions with neutral background gas particles. The electric field of the ion trap is defined by idealized, analytical equations and FFT detection by determining the mirror charge on the cap electrodes is supported.

:doc:`BT-RS-reactiveQITSim <applications/BT_RS_reactiveQITSim>`: Idealized Quadrupole Ion Trap (QIT) with space charge, background gas collisions and reactive ions
    Ion trajectory simulation of an idealized quadrupole ion trap (QIT) considering space charge effects, neutral background particle collisions and chemical reactions of the trapped ions. The electric field of the ion trap is defined by idealized, analytical equations. 

:doc:`BT-generalTrapSim <applications/BT_generalTrapSim>`: RF ion trap with arbitrary ion geometry, space charge and background gas collisions
    Simulation of ion trajectories in an RF ion trap device with arbitrary geometry considering space charge and hard sphere collisions with neutral background gas particles. The electrode geometry is defined by SIMION potential arrays.

.. toctree::
    :maxdepth: 1
    :hidden:

    applications/BT_QITSim
    applications/BT_RS_reactiveQITSim
    applications/BT_generalTrapSim

--------------------------------
``chemistry``: Chemical Kinetics
--------------------------------

:doc:`RS-idealIsothermReactorSim <applications/RS_idealIsothermReactorSim>`: Chemical kinetics in ideally mixed isotherm chemical reactor 
    Particle based kinetics simulation of the reaction of an ensemble of reactive particles with reaction partners in a background gas in an ideally mixed isotherm reactor. 


.. toctree::
    :maxdepth: 1
    :hidden:

    applications/RS_idealIsothermReactorSim


Guide: How to build an IDSimF Application?
==========================================

.. toctree::
    :maxdepth: 1
    
    application_development

.. rubric:: Footnotes

.. [#fn_json_additional] See also `the JSON wikipedia page <https://en.wikipedia.org/wiki/JSON>`_ for further information. 
