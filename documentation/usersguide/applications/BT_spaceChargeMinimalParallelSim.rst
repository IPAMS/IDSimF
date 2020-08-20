.. _application-BT-spaceChargeMinimalParallelSim:

================================
BT-spaceChargeMinimalParallelSim
================================

Minimal parallelized simulation of an ensemble of charged particles interacting by space charge. The particle ensemble can be defined by two ways: 
* The particle ensemble is predefined by an input file which contains the initial positions, initial velocities, charge and mass of the particles in the simulated particle ensemble or
* A random cube with edge length of 3 mm around the origin of the world coordinate system is used for particle initialization. The cube consists of a set of ions with configurable masses. 

Space charge is calculated with a parallelized Barnes-Hut tree. 


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


The ion ensemble is either defined by an input file: 

``ion_cloud_init_file`` : file path
    Path to an ion cloud initialization / definition file 

or a random cube of ions is defined by 

``ion_masses`` : vector of float 
    a vector of ion masses in u

``n_ions`` : vector of integers
    Number of ions with the masses defined in ``ion_masses``
