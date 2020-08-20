.. _application-BT-RS-DMSSim:

============
BT-RS-DMSSim
============

Simulation of ion motion and chemical induced differential mobility in a Differential Ion Mobility (DMS) separation device with idealized planar electrodes, including background gas interaction, ion chemistry and space charge: Ions are drifting in a gap between two planar electrodes on which an asymmetric RF bisinusoidal high voltage waveform, the Separation Voltage (SV) is applied. Due to different ion mobilities in the high field phase and the low field phase of the SV, ions are drifting in the direction perpendicular to the elelctrodes. This drift is countered by an DC Compensation Voltage (CV), which is scanned to get an differential mobility spectrum. 

There are multiple physical and chemical mechanisms leading to differential mobility of ions. This application is currently primarily intended to simulate the chemically induced differential mobility of clustered small ions. The varying electric acceleration in the high and low field phase of the SV induce a variation of the effective ion temperature of the clustered ions. Therefore, the average cluster size is reduced in the high field phase, since high effective temperature leads to declustering.  If the clusters can reform / regrow in the low field phase, this leads to changing average sizes of the clusters between high and low field phases, which in turn leads to different mobility in the high and low field phases. 

The temperature dependent reactions of the individual species are simulated with RS. Thus, the reactions are defined in an :ref:`RS configuration file <usersguide-rs-configurations>`. 

Currently, the mobility of the individual chemical species is **not** dependent on the field or the effective ion temperature. 

The ions are transported through the gap between the electrodes of and planar DMS device by a gas flow. This gas flow can be modeled with different collision models and velocity profiles. The electrodes are parallel to the ``x`` axis, while the separation axis is the ``z`` axis. 

The ions are initialized at random positions in a variable start zone. 


Simulation configuration description
====================================

``sim_time_steps`` : integer
    Number of simulation time steps

``sim_time_steps_per_sv_oscillation`` : integer
    Number of time steps for one complete oscillation of the separation voltage (SV). The time step length in the simulation is set accordingly. 

``sv_frequency_s-1`` : float 
    Frequency of the SV in Hz.

``concentrations_write_interval`` : integer
    Interval, in time steps, between the writes of the species concentration .

``trajectory_write_interval`` : integer
    Interval, in time steps, between writes to the trajectory result file.

``space_charge_factor`` : float
    Multiplication factor for particle-particle interaction (space charge).

``reaction_configuration`` : file path 
    Path to a RS configuration file, defining the chemical reaction system for the simulation. 

``n_ions`` : vector of integers
    Number of particles of the ``discrete`` chemical substances defined in the reaction configuration. The order in this vector is the same as the order of ``discrete`` substances defined in the reaction configuration. 

    Example: 
    If the ``[SUBSTANCES]`` block in the reaction configuration is 

    .. code-block:: 

        [SUBSTANCES]
        Cl_1 discrete 19 1 3.57e-4  4.00000000e-10
        Cl_2 discrete 37 1 2.76e-4  5.17391304e-10
        Cl_3 discrete 55 1 2.35e-4  6.07659574e-10

    the ``n_ions`` vector ``[100, 50, 10]`` will initalize the simulation with 100 particles of ``Cl_1``, 50 of ``Cl_2`` and 10 of ``Cl_3``. 

``start_width_x_mm`` : float 
    Size (in mm) of the ion start zone in ``x`` direction, which is the long axis of the DMS device, parallel to the planar electrodes. 

``start_width_y_mm`` : float
    Size (in mm) of the ion start zone in ``y`` direction.

``start_width_z_mm`` : float
    Size (in mm) of the ion start zone in ``z`` direction, which is the separation axis. 

``electrode_distance_mm`` : float
    Distance (in mm) of the electrodes in the separation (``z``) axis.

``electrode_length_mm`` : float
    Length of the electrodes (in mm) in the long axis of the DMS device (``x`` axis)

``cv_mode`` : Keyword: [``static``, ``auto``]
    Compensation voltage mode: 

    * ``static``: A static CV, defined by ``cv_Vmm-1`` is applied. 
    * ``auto``: Automatic mode to find an optimized CV for maximum transmission through the DMS. The CV is adapted after every SV cycle to minimize the total drift of the simulated particle ensemble in the separation (``z``) axis. The CV is written to a csv file ``<resultbasename>_cv.csv`` after every SV cycle in this mode. ``cv_Vmm-1`` defines the initial CV value in this mode.

``cv_relaxation_parameter`` : 
    Dampening parameter for the minimization process when ``cv_mode`` is ``auto``.

``sv_Vmm-1`` : 
    Separation voltage in V per mm.

``cv_Vmm-1`` : 
    Compensation voltage in V per mm. If ``cv_mode`` is ``auto``, initial CV voltage for the minimization process. 

``collision_model`` : Keyword: [``SDS``, ``none``]
    Sets the used collision / background gas interaction model: 

    * ``SDS``: Statistical Diffusion Simulation model
    * ``none``: No background gas interaction (mostly for testing purposes)

``flow_mode`` : Keyword: [``uniform``, ``parabolic``]
    Sets the gas flow profile for the neutral background gas. 

    * ``uniform``: Uniform flow profile over the whole separation gap. The flow velocity is defined by ``collision_gas_velocity_x_ms-1``. 
    * ``parabolic``: Parabolic flow profile with an average velocity :math:`V_{\text{avg}}` defined by ``collision_gas_velocity_x_ms-1``. The flow velocity vanishes at the electrodes and becomes :math:`2 \cdot V_{\text{avg}}` in the center of the separation gap. 

``background_temperature_mode`` : Keyword:[``isotherm``, ``linear_gradient``]
    Sets the background gas temperature mode. 

    ``isotherm`` : Isotherm mode 
        The background gas temperature is independent of the position and set by

        ``background_temperature_K`` : float 
            Background gas temperature in ``isotherm`` mode
        
    ``linear_gradient`` : Linear temperature gradient
        The background gas temperature is dependent on the position, a linear gradient of the temperature along the long axis of the DMS device is assumed. The temperature gradient is defined by 

        ``background_temperature_start_K`` : float 
            Start temperature of the temperature gradient at the begin of the DMS device, in K.
            
        ``background_temperature_stop_K`` : float
            End temperature of the temperature gradient at the end of the DMS device, in K. 

``background_pressure_Pa`` : float 
    Isotropic pressure of the neutral background gas in Pascal.

``collision_gas_velocity_x_ms-1`` : float
    * Uniform background gas flow velocity in ``x`` direction if ``flow_mode`` is ``uniform``.
    * Average background gas flow velocity in ``x`` direction if ``flow_mode`` is ``parabolic``.

    (in m per second)

``collision_gas_mass_amu`` : float
    Molecular mass of the particles of the background gas in amu.

``collision_gas_diameter_nm`` : float 
    Effective collision diameter of the particles of the background gas in nm. 