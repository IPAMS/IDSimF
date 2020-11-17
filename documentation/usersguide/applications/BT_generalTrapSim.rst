.. _application-BT-generalTrapSim:

=================
BT-generalTrapSim
=================

Ion trajectory simulation in an RF ion trap device with arbitrary geometry considering space charge and hard sphere collisions with neutral background gas particles. The electrode geometry is defined by SIMION potential arrays. 


FT detection
============

There are different detection modes in ion traps. Commonly, ions are ejected mass selectively from the trap, e.g. by resonant excitation, and detected by a particle detector. 

Alternatively, the trapped ions can be excited with an excitation signal on the cap electrodes, which induces an oscillation of the ions in the trap. This oscillation can be detected by the induced mirror charge on detection electrodes. A mass spectrum can be calculated from the recorded mirror charge signal ("transient") by a Fourier transformation. 

The simulation app supports the simulation of this FT detection mode.


Simulation configuration description
====================================

``integrator_mode`` : Keyword:[``verlet``, ``parallel_verlet``]
    Selects the trajectory integrator

    * ``verlet``: Serial (non parallelized) Verlet integrator
    * ``parallel_verlet``: Parallelized Verlet integrator

``sim_time_steps`` : integer
    Number of simulation time steps

``trajectory_write_interval`` : integer
    Interval, in time steps, between writes to the trajectory result file.

``dt`` : float
    Time step length in seconds 

``space_charge_factor`` : float
    Multiplication factor for particle-particle interaction (space charge).

``background_pressure_Pa`` : float
    Pressure of the neutral background as in Pa. 

``background_temperature_K`` : float
    Temperature of the background gas in K. 

``collision_gas_mass_amu`` : float
    Molecular mass of the particles of the background gas in amu.

``collision_gas_diameter_angstrom`` : float
    Effective collision diameter of the particles of the background gas in angström.

``simulation_domain_boundaries`` : vector of vector of floats
    Defines the outer boundaries of the simulation domain around the coordinate system origin, where ions are terminated. Is defined as vector of three two component vectors, defining the minimum and maximum in the spatial dimensions: 
    
    .. code::
        
        [[x low, x high], [y low, y high], [z low, z high]]

--------------------------------------
Ion / simulated particle configuration
--------------------------------------

The particles to simulate can be defined in the simulation configuration file or a predefined particle ensemble can be used which is given as ion cloud file in CSV format. 

Ion Cloud File
--------------

A predefined ion configuration can be specified by 

``ion_cloud_init_file`` : file path
    Path to an ion cloud initialization / definition file 

Ion definition in simulation configuration file
-----------------------------------------------

If no ion cloud file is used, the following configuration parameters define the ion ensemble to simulate: 

``n_ions`` : vector of integers
    Number of ions of the ionic species defined by the the masses defined in ``ion_masses``. 

``ion_masses`` : vector of float 
    Ion masses in amu. 

``ion_charges`` : Vector of float
    Ion charges in elementary charges.     

``ion_collision_gas_diameters_angstrom`` : Vector of float
    Effective hard sphere collision diameters of the ionic species in angström. 

``ion_time_of_birth_range_s`` : float
    Time range in which ions are generated, in seconds. The specified number of ions are generated uniformly in this time range.

Ion start configuration
.......................

The initial positions of the simulated ions can be a cubic box or a cylinder in ``x`` direction. The center of the ion start zone is specified by ``ion_start_base_position_m``.

``ion_start_geometry`` : Keyword:[``box``, ``cylinder``]
    Sets the ion start geometry.

    ``box`` : Ion start zone is a box 
        The ion start zone is a cubic box of 3 mm edge length around ``ion_start_base_position_m``, randomly filled with particles. 

    ``cylinder`` : Ion start zone is a cylinder in ``x`` direction
        The ion start zone is a cylinder parallel to the ``x`` axis, with its center at ``ion_start_base_position_m``. The cylinder is defined by

        ``ion_start_cylinder_radius_m``: float
            Radius of the cylinder around the ``x`` axis in m. 

        ``ion_start_cylinder_length_m`` : float
            Distance from the origin of the cylinder to the cylinder ends in ``x`` direction. The cylinder is therefore in total 2 * ``ion_start_cylinder_length_m`` long. 

``ion_start_base_position_m`` : Vector of 3 floats
    Base position of the ion start zone in m. 


--------------------------------------------------------
Potential Array Configuration / Trap Field Configuration 
--------------------------------------------------------

``f_rf`` : float
    Frequency of the RF trapping field in Hz. 

``V_rf`` : float
    Ground to peak amplitude of the RF trapping field in V. 

``potential_arrays`` : Vector of file paths
    Paths to the SIMION potential array files defining the electric potentials and electrode geometry in the trap. Typically, SIMION potential arrays generated with the *fast adjust* option are used for potential definition. 

    The potential arrays have to have the same geometric extend and are assumed to be normalized. The total potential at a location is calculated by a linear combination of the individual potentials. 

    The file paths are relative to the simulation run configuration file. 

``potential_array_scaling`` : float
    Geometric scaling factor for the potential arrays specified in ``potential_arrays``. 

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


Excitation Field Configuration
------------------------------

A ion excitation field can be applied on selected electrodes. There are two modes of ion excitation: 

* Pulsed excitation with a rectangular excitation pulse of defined length and amplitude applied at the begin of the simulation run.
* Excitation with a given sampled waveform read from a waveform file. 

``excite_potential_factors`` : Vector of floats
    Excitation field factors for the individual electrodes, specified by the potential arrays in ``potential_arrays``. The excitation potential is multiplied with these factors to get the excitation field applied to the individual electrodes.  

``excite_waveform_csv_file`` : File path 
    File path to a file with a sampled excitation waveform. If this parameter is present, the excitation mode is "sampled waveform". 
    
    The waveform file contains one sample per time step and is *not* looped, it is replayed only once at the begin of the simulation run. 
    The sampled waveform is assumed to be normalized, the waveform data is multiplied with "excite_pulse_potential" to calculate the applied excitation potential. 
    
    This file path is relative to the simulation run configuration file. 

``excite_pulse_potential`` : float 
    * When in excitation pulse mode (``excite_waveform_csv_file`` not present in simulation run configuration): Amplitude of rectangular excitation pulse in volt. 
    * When in sampled waveform excitation mode: Multiplication factor for sampled waveform data specified by ``excite_waveform_csv_file`` in volt. 

 ``excite_pulse_length`` : float 
    Length of the rectangular excitation pulse in pulsed excitation mode in seconds. 

--------------------------
FT detection configuration
--------------------------

The FFT result records the total induced mirror current of the simulated ion ensemble on a set of detection electrodes, specified by ``detection_potential_factors``. 

``fft_write_interval`` : integer 
    Interval, in time steps, between samples for the FFT result file, which records a simulated transient for Fourier transformation. 

``fft_write_mode`` : Keyword:[``unresolved``, ``mass_resolved``]
    Selects the mode in which the FFT result file is written. 

    * ``unresolved``: The FFT transient signal is calculated from the whole simulated particle ensemble. This is the signal which would be detectable in a physical experiment. 
    * ``mass_resolved``: Currently not implemented 

``detection_potential_factors`` : Vector of float
    Mirror charge detection factors. The induced mirror current on the electrodes which are described by the potential arrays in ``potential_arrays``, is multiplied by this factor to get the contribution to the total induced mirror current. 
