.. _application-BT-spaceChargeSimpleSim:

=======================
BT-spaceChargeSimpleSim
=======================

Simple simulation of an random ensemble of charged particles interacting by space charge. The particle ensemble is defined by a random cube with edge length of 3 mm around the origin of the world coordinate system. The cube consists of a set of ions with configurable masses. 

Space charge is calculated with a serial version of a Barnes-Hut tree. 


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

``ion_masses`` : vector of float 
    a vector of ion masses in u

``n_ions`` : vector of integers
    Number of ions with the masses defined in ``ion_masses``
