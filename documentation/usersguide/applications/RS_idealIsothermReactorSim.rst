.. _application-RS-idealIsothermReactorSim:

==========================
RS-idealIsothermReactorSim
==========================

Particle based kinetics simulation of the reaction of an ensemble of reactive particles with reaction partners in a background gas in an ideally mixed isotherm reactor. 

The chemical kinetics is simulated with :doc:`RS <../rs>`. The reaction system, chemical species parameters and reactions, is defined in an RS configuration file specified in ``reaction_configuration``. 


Simulation configuration description
====================================

``sim_time_steps`` : integer
    Number of simulation time steps

``dt_s`` : float 
    Time step length in seconds

``n_particles`` : Vector of integers
    Number of particles of the ``discrete`` chemical substances defined in the reaction configuration. The order in this vector is the same as the order of ``discrete`` substances defined in the reaction configuration. 

    Example: 
    If the ``[SUBSTANCES]`` block in the reaction configuration is 

    .. code-block:: 

        [SUBSTANCES]
        Cl_1 discrete 19 1 3.57e-4  4.00000000e-10
        Cl_2 discrete 37 1 2.76e-4  5.17391304e-10
        Cl_3 discrete 55 1 2.35e-4  6.07659574e-10

    the ``n_ions`` vector ``[100, 50, 10]`` will initalize the simulation with 100 particles of ``Cl_1``, 50 of ``Cl_2`` and 10 of ``Cl_3``. 


``background_temperature_K`` : float
    Isotropic background temperature in K. 

``reaction_configuration`` : file path 
    Path to a RS configuration file, defining the chemical reaction system for the simulation. The file path is relative to the simulation run config file. 

``concentrations_write_interval`` : integer
    Interval, in time steps, between the writes of the species concentration to the concentrations result file.
