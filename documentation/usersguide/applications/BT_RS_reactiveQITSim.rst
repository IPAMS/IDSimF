.. _application-BT-RS-reactiveQITSim:

====================
BT-RS-reactiveQITSim
====================

Ion trajectory simulation in an ideal Quadrupole Ion Trap (QIT) considering chemical reactions of the simulated ions. 

The field definition and basic theory is the same as in 
:doc:`BT_QITSim`. In addition to the trajectory integration considering space charge and hard sphere background gas collisions, chemical reactions are simulated with :doc:`RS <../rs>`. Since the simulation is coupled to RS, most of the particle parameters are defined in the RS reaction system file referenced in ``reaction_configuration``. 

The application allows different reaction types (static, temperature dependent, etc.) in the chemical kinetics simulation. Particularly, *collision* based reactions can be simulated: For this reaction type, possible reaction events takes place, when ions collide with neutral background gas particles in the hard sphere collision simulation. 

Simulation configuration description
====================================

``sim_time_steps`` : integer
    Number of simulation time steps

``trajectory_write_interval`` : integer
    Interval, in time steps, between writes to the trajectory result file.

``concentrations_write_interval`` : integer
    Interval, in time steps, between the writes of the species concentration.

``fft_write_interval`` : integer 
    Interval, in time steps, between samples for the FFT result file, which records a simulated transient for Fourier transformation. 

``fft_mode`` : Keyword:[``off``, ``unresolved``, ``mass_resolved``]
    Selects the mode in which the FFT result file is written. The FFT result records the average velocity of the simulated particle ensemble in the ``z`` direction as an approximate induced mirror charge on the cap electrodes. 

    * ``off``: FFT analysis is switched off, no FFT file is written
    * ``unresolved``: The FFT transient signal is calculated from the whole simulated particle ensemble. This is the signal which would be detectable in a physical experiment. 
    * ``mass_resolved``: The FFT transient signal is recorded for every ionic species, with distinct molecular mass, separately. This allows to see individual transients for the species. 

``dt`` : float 
    Time step length in seconds

``n_ions`` : Vector of integers
    Number of particles of the ``discrete`` chemical substances defined in the reaction configuration. The order in this vector is the same as the order of ``discrete`` substances defined in the reaction configuration. 

    Example: 
    If the ``[SUBSTANCES]`` block in the reaction configuration is 

    .. code-block:: none

        [SUBSTANCES]
        Cl_1 discrete 19 1 3.57e-4  4.00000000e-10
        Cl_2 discrete 37 1 2.76e-4  5.17391304e-10
        Cl_3 discrete 55 1 2.35e-4  6.07659574e-10

    the ``n_ions`` vector ``[100, 50, 10]`` will initalize the simulation with 100 particles of ``Cl_1``, 50 of ``Cl_2`` and 10 of ``Cl_3``. 

``start_width_m`` : float
    Ions are initialized randomly within a cubic box around the trap center. This is the edge length of this box in m. 

``background_temperature_K`` : float
    Temperature of the background gas in K. 

``space_charge_factor`` : float
    Multiplication factor for particle-particle interaction (space charge).

``max_ion_radius`` : float
    Maximum radius a simulated particle can have before it is terminated in the simulation, in m. 

------------------------------------------
Chemistry and background gas configuration
------------------------------------------

``reaction_configuration`` : file path 
    Path to a RS configuration file, defining the chemical reaction system for the simulation. 

The interaction between the simulated ions and the neutral background gas is simulated with hard sphere collisions. The background gas can be a mixture of individual components.

``collision_gas_names`` : Vector of strings
    Names of the background gas components. 

    .. note::
        Due to the coupling to the chemical reaction model, the background gas componentes have to exist as ``isotropic`` components in the RS configuration. 

``collision_gas_masses_amu`` : Vector of float
    Molecular masses of the particles of the background gas mixture components in amu.

``collision_gas_diameters_angstrom`` : Vector of float
    Effective collision diameters of the particles of the background gas components in Angström. 

``partial_pressures_Pa`` : vector of float 
    Partial pressures of the individual components of the background gas mixture in Pascal.

.. warning::
    There is *no* check if the reaction configuration and the parameters specified here are physically and chemically consistent. One has to check if the speciefied system is actually possible. 


---------------------------
Trap geometry configuration
---------------------------

``geometry_mode`` : Keyword:[``default``, ``scaled``, ``variable``]
    Selects trap geometry mode. The trap geometry defines the electric field in the trap and is defined by the ring electrode radius ``r_0``  and the cap distance ``z_0``. 

    ``default`` : Default trap with :math:`r_0= 10 \text{mm}`
        The default geometry is a typical small commerical QIT with :math:`r_0 = 10 \text{mm}` and :math:`z_0 = 7 \text{mm}` which is approximately (within 2%) fulfilling the ideal relationship :math:`r_0^2 = 2 z_0^2`. 

    ``scaled`` : Scaled default trap 
        The default trap geometry scaled by a factor ``geometry_scale``:

        ``geometry_scale`` : float
            Geometric scaling factor for a scaled default trap. 

    ``variable`` : Variable geometry
        Fully variable geometry, :math:`r_0` and :math:`z_0` can be configured freely: 

        ``r_0``: float 
            :math:`r_0` in meter. 
        
        ``z_0``: float 
            :math:`z_0` in meter. 


------------------------
Trap field configuration
------------------------

``f_rf`` : float
    Frequency of the RF trapping field in Hz. 

Trap field RF voltage
---------------------

The RF trap field voltage can be static or can be ramped during the simulation. 

Static field mode: 

``rf_V`` : float
    Ground to peak trap field amplitude in V. 

Ramped field mode: 

``rf_ramp_start_V`` : float 
    Ground to peak field amplitude at the start of the amplitude ramp, in V. 

``rf_ramp_stop_V`` : float
    Final ground to peak field amplitude at the end of the amplitude ramp, in V. 

``rf_ramp_waiting_timesteps`` : integer
    Number of time steps to wait on ``rf_ramp_start_V`` before starting the amplitude ramp. 

Ion excitation field
--------------------

The trapped ions can be excited by an bipolar field applied to the cap electrodes. 

``excite_mode`` : Keyword:[``off``, ``rect_pulse``, ``waveform``, ``continuous_sine``]

    Selects the ion excitation mode. 

    ``off`` : No excitation 
        No bipolar field is applied at all. 

    ``rect_pulse`` : Rectangular excitation pulse 
        Applies a rectangular excitation pulse to the cap electrodes at the beginning of the simulation. The amplitude of the pulse is defined by ``excite_potential``, the duration is defined by 

        excite_pulse_length : float 
            Length of the excitation pulse in seconds. 
    
    ``waveform`` : Excitation with sampled waveform 
        Applies an excitation with a given sampled waveform read from a waveform file.

        ``excite_waveform_csv_file`` : File path 
            File path to a file with a sampled excitation waveform. 
            
            The waveform file contains one sample per time step and is *not* looped, it is replayed only once at the begin of the simulation run. 
            The sampled waveform is assumed to be normalized, the waveform data is multiplied with "excite_pulse_potential" to calculate the applied excitation potential. 
    
            This file path is relative to the simulation run configuration file. 

    ``continuous_sine`` : Continuous sinusoidal excitation
        Applies a continuous sinusoidal excitation on the cap electrodes with a frequency which is a fraction of the main RF frequency. The ground to peak amplitude of the applied excitation field is defined by ``excite_potential``. 
        
        The frequency is defined by 

        ``excite_divisor`` : float 
            Frequency divisor. The excitation field frequency is the main trap field RF frequency devided by this devisor. 


``excite_potential`` : float 
    Excitation potential / excitation scaling multiplicator. See ``excite_mode`` for details. 






