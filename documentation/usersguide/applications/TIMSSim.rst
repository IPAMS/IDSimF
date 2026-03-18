.. _application-TIMSSim:

========
TIMSSim
========

Simulates the trajectories of charged particles in a Trapped Ion Mobility Spectrometry (TIMS) Device, including background gas interaction, ion chemistry and space charge.
A TIMS device consists of a set of electrodes forming three regions: the entrance funnel, the TIMS tunnel and the exit funnel. Ions are trapped by an electric field which counteract the drag force of a flow of gas. The magnitude of the axial electric field is progressively decreased, allowing for ions to be diluted from the device.

The electrode geometry and potentials are defined through SIMION potential arrays. 

There are four different electric potentials which can be applied to the electrodes: 

    * ``DC``: Static voltages which are not changing throughout the whole simulation
    * ``RF``: High frequency RF field (sine) for radial ion confinement
    * ``Scan ramp``: Potential which is linearly ramped to a final value, for mobility scan 
    * ``Gate``: Potential to gate ions into TIMS region, gate voltage is switched on after a delay time and is held for a gate open time
  
Those electric potentials are configured in a separate CSV file for every individual potential array. See :ref:`tims-voltage-config` for details.

Ions start uniformly distributed in a configurable start zone. The interaction between background gas and the ions can be described with different collision models. Electric ion-ion interaction (space charge) can be modeled with a parallelized version of a Barnes-Hut tree. 

The modeled ions are chemically reactive. The temperature dependent reactions of the individual chemical species are simulated with RS. The reactions are defined in an :ref:`RS configuration file <usersguide-rs-configurations>`. To simulate chemically inert ions, an RS configuration without reactions between the species has to be defined.


Simulation configuration description
====================================

------------------
General Parameters
------------------

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
    

.. _tims-voltage-config:

------------------------------------------
Potential Array and Voltages Configuration
------------------------------------------

The electrode geometry and electric field shape of the TIMS analyzer are defined through SIMION potential arrays.  The potential arrays and actual voltages applied to the individual electrodes are defined in a separate CSV file (given in  ``potential_configuration``). 

This file has the following form: 

.. code-block:: shell

    # potential name; static DC potential; RF factor; ramp gradient factor; gate factor
    tims_geometry.pa1;  0;    1; 0.1; 0
    tims_geometry.pa2; -2.5; -1; 0.2; 0
    tims_geometry.pa3; -5.0;  1; 0.1; 0
    tims_geometry.pa4; -7.5; -1; 0.2; 0

The individual columns are: 

    * ``potential name``: Name / Paths of the potential array (PA) file. The file paths are relative to the simulation run configuration file. 
    * ``static DC potential`` :math:`U_{\text{dc}}`: Static (DC), absolute potential for this PA 
    * ``RF factor`` :math:`f_{\text{RF}}`: RF amplitude for this PA, given as factor relative to the absolute RF amplitude defined in ``confining_RF_amplitude_V``
    * ``ramp gradient factor`` :math:`f_{\text{Ramp}}`: Amplitude of the linear scan ramp for this PA, given as factor relative to the absolute scan ramp amplitude defined by ``gradient_voltage_V`` and ``gradient_ramp_velocity_V/ms``
    * ``gate factor`` :math:`f_{\text{Gate}}`: Amplitude of the gate voltage for this PA, given as factor relative to the absolute gate voltage defined by ``gate_voltage_V``

Typically, SIMION potential arrays generated with the *fast adjust* option are used for the electrode definitions. The potential arrays have to have the same geometric extend and are assumed to be normalized. The total potential at a location is calculated by a linear combination of the individual potentials. 

The total voltage for a PA, :math:`U_{\text{PA}}`, is calculated by: 

.. math::

    U_{\text{PA}} = U_{\text{dc}} + f_{\text{RF}} \, U_{RF}(t) + f_{\text{Ramp}} \, U_{Ramp}(t)+ f_{\text{Gate}} \, U_{Gate}(t)

with the time dependent RF, ramp and gate voltages :math:`U_{RF}(t)`, :math:`U_{Ramp}(t)`, :math:`U_{Gate}(t)`. 

``potential_array_scale`` : float
    Geometric scaling factor for the potential arrays specified in ``potential_arrays``.

``potential_configuration`` : file path 
    Path to the potential definition CSV file

RF Field Configuration
----------------------

``confining_RF_amplitude_V`` : float
    Peak-to-peak amplitude (in V) of the confining voltage meant to reduce radial ion drift.

``confining_RF_frequency_Hz`` : float
    Frequency (in Hz) of the confining voltage.


Scan Ramp  / Gradient Configuration
-----------------------------------

The scan ramp / gradient is a linear potential ramp, which starts after a delay time and increases with a fixed rate to the end value. 

``gradient_start_time_s`` : float
    Delay time until the scan gradient starts (in seconds)

``gradient_ramp_velocity_V/ms`` : float
    Slope of the scan gradient / scan ramp (in volts per milliseconds)

``gradient_voltage_V`` : float
    End voltage of the scan gradient ramp (volt)

Gate Configuration
------------------

The gate voltage is applied after a delay time and is kept for a gate open duration time 

``gate_open_time_s`` : float
    Delay time until the gate is opened (in seconds)

``gate_open_duration_s`` : float
    Opening time of the gate (in seconds)

``gate_voltage_V`` : float
    Voltage applied in the "gate open" state (volt)


-----------------------------------------------
Collision models and background gas interaction 
-----------------------------------------------

The simulation has different modes to model the interactions between ions and the background gas which are suitable for different background gas pressure ranges. The model can use external flow fields (mostly from fluid dynamic simulations) or simple uniform or parabolic flow profiles can be assumed. The external flow fields are imported from flow data in SIMION potential array (PA) files. Such files can be generated from CFD solver data with conversion scripts (e.g. `comsol_to_pa.lua`) provided by the SIMION distribution. 

.. note::

    It is planned to integrate other file formats for CFD solution import.

The collision model mode is controlled by the ``collision_model`` parameter: 

``collision_model`` : keyword [``SDS``, ``HS``, ``MD``, ``none``]
    Sets the used collision / background gas interaction model: 

    * ``SDS``: Statistical Diffusion Simulation model
    * ``HS``: Hard Sphere model
    * ``MD``: Molecular Dynamics model
    * ``none``: No background gas interaction (mostly for testing purposes)
  
``collision_gas_mass_amu`` : float
    Molecular mass of the particles of the background gas in amu.

``collision_gas_diameter_nm`` : float 
    Effective collision diameter of the particles of the background gas in nm. 
  
``flow_mode``: keyword [``uniform``, ``parabolic``]
    Sets the background gas flow mode: 

    * ``uniform``: Uniform flow velocity in ``x`` direction (default value)
    * ``parabolic``: Parabolic flow velocity profile in ``x`` direction
    * ``static_field``: External, static, flow profile (mostly from CFD solution)
    
    The parameter is *optional*, if it is omitted, ``uniform`` flow profile is assumed.

Uniform / parabolic background gas flow profiles
------------------------------------------------

Uniform or simple parabolic flow profiles are defined by the following parameters:

``background_pressure_Pa`` : float 
    Isotropic pressure of the neutral background gas in Pascal.
    
``background_temperature_K``: float
    Background gas temperature in Kelvin.

``background_velocity_x_ms-1`` : float
    Background gas velocity in ``x`` direction (meter per second). For uniform flow this is the uniform background gas velocity in ``x`` direction. For parabolic flow, this is the average flow velocity. The maximum flow velocity in the center of the parabolic flow profile is a factor of two higher (2 * background_velocity_x_ms-1).

``flow_profile_maximum_radius_m`` : float
    *For parabolic flow profile only*: With of the parabolic flow profile in ``y``-``z`` direction (in meter).

External flow profiles (CFD solutions)
--------------------------------------

Imported flow fields can be full 3d or 2d axialsymmetric. In the 2d axial symmetric case, the ``x`` axis is the symmetry axis and the flow PAs have to have the correct symmetry. 

``flow_field``: Vector of file paths
    Flow velocity component fields. For 3d case: ``x``, ``y`` and ``z`` components in three separated potential array files. For a 2d axial symmetric case ``x`` and ``r`` (the radial flow component) as two separatd potential array files. 

    The PAs are assumed to be in meter per second (m/s).

``pressure_field``: Vector of file paths
    Pressure field in Pascal (Pa).

``temperature_field``: Vector of file paths
    Temperature field in Kelvin (K).


.. note::

    The pressure and temperature fields have to be provided in a vector (within square brackets) even if they are single PA files. 
