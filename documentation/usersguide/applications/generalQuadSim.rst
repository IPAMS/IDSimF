.. _application-generalQuadSim:

==============
generalQuadSim
==============

Simulates ion trajectories in an RF only quadrupole device with arbitrary electrode geometry considering space charge and interaction with a flowing background gas. The electrode geometry and resulting electric field is defined by a set of SIMION potential arrays. The pressure / density and flow velocity of the background gas can be defined by spatially resolved fields, e.g. from an CFD simulation of the simulated device. 

.. warning::
    This application has to be revised. It can fail with an exception under certain conditions. 

.. note:: 
    Currently, no collisions of simulated particles with electrode surfaces are considered in this application and the RF frequency of the device is currently fixed at 1 MHz. 


Simulation configuration description
====================================

``sim_time_steps`` : integer
    Number of simulation time steps

``trajectory_write_interval`` : integer
    Interval, in time steps, between writes to the trajectory result file.

``dt`` : float
    Time step length in seconds.

``n_ions`` : vector of integers
    Number of ions with the masses defined in ``ion_masses``

``ion_masses`` : vector of float 
    Vector of ion masses in amu.

``V_rf`` : float
    RF amplitude in V.

``V_entrance`` : float
    DC amplitude for the entrance electrode into the quadrupole in V.

``P_factor`` : float
    Scaling factor for the pressure field. 

``space_charge_factor`` : float
    Multiplication factor for particle-particle interaction (space charge) 

``collision_gas_mass_amu`` : float
    Molecular mass of the particles of the background gas in amu.

``collision_gas_diameter_angstrom`` : float
    Effective collision diameter of the particles of the background gas in angström. 

``background_temperature`` : 
    Isotropic temperature of the background gas in K. 

``start_q_length_mm`` : float
    Length of the ion start zone in axial direction in mm.

``entrance_aperture_mm`` : float
    With of the quadratic cross section of the ion start zone in mm. 

``q_start_mm`` : float
    Begin of the ion start zone in axial (``x``) direction in mm.

``restart_q_length_mm`` : float
    Length of the ion start zone for restarted ions in axial (``x``) direction in mm. 

``max_q_length_mm`` : float
    End position in axial (``x``) direction, where ions are terminated  in mm.

``max_r_mm`` : float
    Maximal distance in radial direction from the center line of the quad where ions are terminated. 

--------------------------------------------------
Paths to gas parameter fields and potential arrays
--------------------------------------------------

.. note::
    All file paths are relative to the simulation configuration file. 

``rho_field_file`` : File path 
    Path to background gas density (rho) field.

``flow_field_file`` : File path
    Path to background gas flow velocity field.

The potential arrays have to have the same geometric extend and are assumed to be normalized. The total potential at a location is calculated by a linear combination of the potentials: 

.. math::

    U = U_{\text{RF norm.}} \cdot V_{\text{RF}} + U_{\text{entrance norm.}} \cdot V_{\text{entrance}}

``electric_field_rf_file`` : File path
    Path to potential array for RF electrodes with normalized bipolar field for RF electrodes. 

``electric_field_entrance_file`` : File path
    Path to the normalized potential array for the entrance electrode. 
