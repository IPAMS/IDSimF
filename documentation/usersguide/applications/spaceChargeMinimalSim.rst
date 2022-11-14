.. _application-spaceChargeMinimalSim:

=====================
spaceChargeMinimalSim
=====================

Simulates the electric particle-particle interactions (space-charge) in an defined ensemble of charged particles. The particle ensemble is predefined by an input file which contains the initial positions, initial velocities, charge and mass of the particles in the simulated particle ensemble. Space charge is calculated with a serial version of a Barnes-Hut tree. 


Simulation configuration parameters 
===================================

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

