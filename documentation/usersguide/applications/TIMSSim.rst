.. _application-TIMSSim:

========
TIMSSim
========

Simulates the trajectories of charged particles in a Trapped Ion Mobility Spectrometry (TIMS) Device, including background gas interaction, ion chemistry and space charge.
A TIMS device consists of a set of electrodes forming three regions: the entrance funnel, the TIMS tunnel and the exit funnel. Ions are trapped by an electric field which counteract the drag force of a flow of gas. The magnitude of the axial electric field is progressively decreased, allowing for ions to be diluted from the device.

The electrode geometry and potentials are defined through SIMION potential arrays. An RF voltage is applied to the TIMS analyzer electrodes to provide radial ion confinement. In addition to this an axial electric field gradient is created by applying DC potentials to the electrodes.
Ions start uniformly distributed in a configurable start zone. The interaction between background gas and the ions can be described with different collision models. Electric ion-ion interaction (space charge) can be modeled with a parallelized version of a Barnes-Hut tree. 

The modeled ions are chemically reactive. The temperature dependent reactions of the individual chemical species are simulated with RS. The reactions are defined in an :ref:`RS configuration file <usersguide-rs-configurations>`. To simulate chemically inert ions, an RS configuration without reactions between the species has to be defined.


Simulation configuration description
====================================


``sim_time_steps`` : integer
    Number of simulation time steps

``dt_s`` : float
    Time step length in seconds 

``concentrations_write_interval`` : integer
    Interval, in time steps, between the writes of the species concentration .

``trajectory_write_interval`` : integer
    Interval, in time steps, between writes to the trajectory result file.

``trajectory_write_velocities`` : boolean
    if ``true``: Particle velocities are written to the auxiliary data in the trajectory result file. 

``space_charge_factor`` : float
    Multiplication factor for particle-particle interaction (space charge).

``reaction_configuration`` : file path 
    Path to a RS configuration file, defining the chemical reaction system for the simulation. This file path is interpreted relatively to the simulation run configuration file.

``n_particles`` : vector of integers
    Number of particles of the ``discrete`` chemical substances defined in the reaction configuration. The order in this vector is the same as the order of ``discrete`` substances defined in the reaction configuration. 

    Example: 
    If the ``[SUBSTANCES]`` block in the reaction configuration is 

    .. code-block:: none

        [SUBSTANCES]
        Cl_1 discrete 19 1 3.57e-4  4.00000000e-10
        Cl_2 discrete 37 1 2.76e-4  5.17391304e-10
        Cl_3 discrete 55 1 2.35e-4  6.07659574e-10

    the ``n_ions`` vector ``[100, 50, 10]`` will initalize the simulation with 100 particles of ``Cl_1``, 50 of ``Cl_2`` and 10 of ``Cl_3``. 

``start_box_dimensions_mm`` : Vector of floats
    x-, y-, and z-dimensions (in mm) of the ion start zone.
    
``start_box_position_mm`` : Vector of floats
    x-, y-,  and z-coordinates (in mm) of the corner of the ion start zone.
    
``simulation_domain_boundaries`` : Vector of vector of floats
    Defines the outer boundaries of the simulation domain around the coordinate system origin, where ions are terminated. Is defined as vector of three two component vectors, defining the minimum and maximum in the spatial dimensions: 
    
        .. code::
        
            [[x low, x high], [y low, y high], [z low, z high]] 
    
``confining_RF_amplitude_V`` : float
    Peak-to-peak amplitude (in V) of the confining voltage meant to reduce radial ion drift.

``confining_RF_frequency_Hz`` : float
    Frequency (in Hz) of the confining voltage.
    
``potential_gradient_filename`` : file path
    Path to a CSV file defining the potential gradient applied to the electrodes. The file is expected to be single column with one value each row.
    
``potential_gradient``: Vector of floats
    Potential gradient that is the be applied to the electrodes. This gradient can be defined either through a CSV file as described above or directly through a vector in the input JSON file. The application will always check for a file path to a CSV first. If none is given, it will then check for a vector containing the gradient. If both a file path and a vector are provided, the file path will take precedence.
    
-----------------------------
Potential Array Configuration
-----------------------------

The electrode geometry and potentials of the TWIMS electrode stack are defined through SIMION potential arrays. As TWIMS utilizes two different voltage patterns applied to the ring electrodes, the transient wave pattern and a confining RF-voltage, two sets of potential arrays are required.

The Traveling Wave pattern can be created across a variable number of electrodes that is defined by the number of potential arrays given, i.e. if 8 potential array files are given it is assumed that a single wave pattern encompasses 8 ring electrodes. 
The confining field is meant to prevent radial ion drift and is created by applying opposing phases of a RF-voltage to two adjacent electrodes, thus two potential arrays are required for the RF-field.

``potential_array_scale`` : float
    Geometric scaling factor for the potential arrays specified in ``potential_arrays``.
    
``potential_arrays`` : Vector of file paths
    Paths to the SIMION potential array files defining the electrode geometry of the stack and the potentials of the transient wave. Typically, SIMION potential arrays generated with the *fast adjust* option are used for potential definition. 

    The potential arrays have to have the same geometric extend and are assumed to be normalized. The total potential at a location is calculated by a linear combination of the individual potentials. 

    The file paths are relative to the simulation run configuration file. 

-----------------------------------------------
Collision models and background gas interaction 
-----------------------------------------------

The simulation has different modes to model the interactions between ions and the background gas which are suitable for different background gas pressure ranges. 

The collision model mode is controlled by the ``collision_model`` parameter: 

``collision_model`` : keyword [``SDS``, ``HS``, ``MD``, ``none``]
    Sets the used collision / background gas interaction model: 

    * ``SDS``: Statistical Diffusion Simulation model
    * ``HS``: Hard Sphere model
    * ``MD``: Molecular Dynamics model
    * ``none``: No background gas interaction (mostly for testing purposes)

``background_pressure_Pa`` : float 
    Isotropic pressure of the neutral background gas in Pascal.
    
``background_temperature_K``: float
    Background gas temperature in Kelvin.
    
``collision_gas_velocity_x_ms-1`` : float
    Uniform background gas velocity in ``x`` direction in m per second.

``collision_gas_mass_amu`` : float
    Molecular mass of the particles of the background gas in amu.

``collision_gas_diameter_nm`` : float 
    Effective collision diameter of the particles of the background gas in nm. 
