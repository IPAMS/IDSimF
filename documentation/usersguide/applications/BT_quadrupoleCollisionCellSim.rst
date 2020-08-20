.. _application-BT-quadrupoleCollisionCellSim:

=============================
BT-quadrupoleCollisionCellSim
=============================

Simulates a quadrupolar collision cell with hard sphere collisions between ions and background gas, space charge and variable electrode geometry given by SIMION potential arrays. The ions are accelerated into the collision cell which is filled with a collision gas. The ions are confined radially in a quadrupolar RF field. 

.. note::
    * Currently, the collision of ions with electrode surfaces is not determined and ions are not terminated on electrode surfaces. 
    
    * The time step length is *not* adapted to the gas dynamic parameters of the hard sphere model model. For a valid modeling, the time step length should be significantly shorter than the mean time between ion-neutral collisions.


Simulation configuration description
====================================

``sim_time_steps`` : ion_start_geometry
    Number of simulation time steps

``dt`` : float
    Time step length in seconds 

``trajectory_write_interval`` : integer
    Interval, in time steps, between writes to the trajectory result file.

``space_charge_factor`` : float
    Multiplication factor for particle-particle interaction (space charge).

``frequency_rf`` : float
    Frequency of the RF trapping field in Hz. 

``n_ions`` : Vector of integers
    Number of ions of the ionic species defined by the the masses defined in ``ion_masses``, the collision diameters in ``ion_collision_gas_diameters_angstrom`` and the charges in ``ion_charges``. 

``ion_masses`` : Vector of float 
    Ion masses in amu. 

``ion_charges`` : Vector of float
    Ion charges in elementary charges. 

``ion_collision_gas_diameters_angstrom`` : Vector of float
    Effective hard sphere collision diameters of the ionic species in angström. 

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


-----------------------
Ion start configuration
-----------------------

The simulated particle ensemble can be initialized in different start geometries, which have different configuration parameters: 

``ion_start_geometry`` : Keyword:[``box``, ``cylinder``, ``vector``]
    Sets the ion start geometry mode.

    ``box`` : Ion start zone is a box 
        The ion start zone is a box, oriented parallel to the axes. The box is defined by: 
        
        ``ion_start_box_center_position_m`` : vector of 3 floats
            Position of the box center in m. 

        ``ion_start_box_size_m`` : vector of 3 floats
            Size of the box in ``x``, ``y``, ``z`` direction in m. 

    ``cylinder`` : Ion start zone is a cylinder in ``x`` direction
        The ion start zone is a cylinder parallel to the ``x`` axis, with its center at the origin of the coordinate system. 

        ``ion_start_cylinder_radius_m``: float
            Radius of the cylinder around the ``x`` axis. 

        ``ion_start_cylinder_length_m`` : float
            Length from the origin of the cylinder to the cylinder end in ``x`` direction. The cylinder is in total 2 * ``ion_start_cylinder_length_m`` long. 

    ``vector`` : single vector
        All particles start at a single position given by

        ``ion_start_vector_position_m``: vector of 3 floats
            Position to start the particles at. 

        Note that starting two particles at the exact same position will lead to an exception. Therefore, this mode can only used with single particles. 