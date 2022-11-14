.. _application-quadrupoleCollisionCellSim:

==========================
quadrupoleCollisionCellSim
==========================

Simulates a quadrupolar collision cell with hard sphere collisions between ions and background gas, space charge and variable electrode geometry given by SIMION potential arrays. The ions are accelerated into the collision cell which is filled with a collision gas. The ions are confined radially in a quadrupolar RF field. 

.. note::
    
    * The time step length is *not* adapted to the gas dynamic parameters of the hard sphere model model. For a valid modeling, the time step length should be significantly shorter than the mean time between ion-neutral collisions.


Simulation configuration description
====================================

``sim_time_steps`` : ion_start_geometry
    Number of simulation time steps

``record_mode`` : string
    Particle data record mode: Defines if a simple or an extended set of particle attributes is recorded in the trajectory files.

    Valid modes are: 

    * ``full`` mode records 
        * velocity
        * electric field acting on the particles and 
        * space charge force acting on the particles
    * ``simple`` records only the velocity

    The recorded particle attributes are recorded as separated (:math:`x,y,z`) components. 

``dt`` : float
    Time step length in seconds 

``trajectory_write_interval`` : integer
    Interval, in time steps, between writes to the trajectory result file.

``space_charge_factor`` : float
    Multiplication factor for particle-particle interaction (space charge).

``frequency_rf`` : float
    Frequency of the RF trapping field in Hz. 

``V_rf`` : float
    RF amplitude in V. 

``collision_gas_mass_amu`` : float
    Molecular mass of the particles of the background gas in amu.

``collision_gas_diameter_angstrom`` : float
    Effective collision diameter of the particles of the background gas in angström.

``background_gas_temperature_K`` : float
    Isotropic temperature of the background gas in K. 

``background_gas_pressure_Pa`` : float
    Pressure of the background / collision gas in Pa. 

``termination_mode`` : string
    Particle termination mode: Defines if simulated particles should be terminated or restarted when leaving the simulation domain or hitting an electrode.

    Valid modes are: 

    * ``terminate`` Particles are terminated and not restarted 
    * ``restart`` Particles are restarted in their start zone

    .. note::

        Currently ``restart`` can not be used with an particle cloud start file, only with an particle start zone definition. 
    
``simulation_domain_boundaries`` : vector of vector of floats
    Defines the outer boundaries of the simulation domain around the coordinate system origin, where ions are terminated. Is defined as vector of three two component vectors, defining the minimum and maximum in the spatial dimensions: 
    
    .. code::
        
        [[x low, x high], [y low, y high], [z low, z high]]

----------------------------------------------
Potential array / electric field configuration 
----------------------------------------------

The electrode geometry and the resulting electric potential is defined by SIMION potential arrays. Every potential array contains the potential distribution resulting from an individual electrode on a normalized potential and every other electrode on ground. The total potential is calculated by a linear combination of the individual normalized potentials. 

.. Warning::
    The potential arrays are assumed to be generated with the *fast adjustment* feature of SIMION with the same physical extent, resolution and coordinate system origin. 

``potential_arrays`` : Vector of files paths 
    Paths to a set of potential potential arrays defining the geometry and the electric fields. 

``potential_array_scale`` : float
    Geometric scaling factor for the potential arrays.

``dc_potentials`` : Vector of float
    Invariant (DC) potentials on the electrodes defined by the potential arrays in ``potential_arrays`` in V. 

``rf_potential_factors`` : Vector of float
    Factors defining the applied RF amplitude on the electrodes defined by the potential arrays in ``potential_arrays``. 

    The applied voltage on the individual electrode is 

    .. math::

        U = U_{\text{DC}} + \cos(\omega \cdot t) \cdot  V_{\text{RF}} \cdot F_{\text{RF}}

    with

    * :math:`t` the current time in the simulation
    * :math:`V_{\text{RF}}` given by ``V_rf``
    * the angular frequency :math:`\omega`, given by :math:`2\pi\cdot` ``frequency_rf``
    * :math:`F_{\text{RF}}` given by ``rf_potential_factors``.


.. include:: includes/ion_definition_reading_options.rst

