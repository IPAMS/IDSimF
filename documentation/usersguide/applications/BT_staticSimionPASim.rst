.. _application-BT-staticSimionPASim:

====================
BT-staticSimionPASim
====================

Simulates the particle trajectories of an ensemble of charged particles in an arbitrary electric potential distribution and electrode geometry, defined by a SIMION potential array. The particles interact by electric particle-particle force (space charge) which is calculated by an serial version of a Barnes-Hut tree. 


Simulation configuration description
====================================

``sim_time_steps`` : integer
    Number of simulation time steps

``trajectory_write_interval`` : integer
    Interval, in time steps, between writes to the trajectory result file

``dt`` : float
    Time step length in seconds

``space_charge_factor`` : float
    Multiplication factor for particle-particle interaction (space charge)

``ion_cloud_init_file`` : file path
    Path to an ion cloud initialization / definition file

``potential_array_file``: file path
    Path to an potential array file defining the electrode geometry / electric potentials in the simulation domain