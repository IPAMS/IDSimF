.. _application-RS-DMSSimplifiedSim:

===================
RS-DMSSimplifiedSim
===================

Simplified simulation of ion motion and chemical induced differential mobility in a Differential Ion Mobility (DMS) separation device with idealized planar electrodes: Ions are drifting in a gap between two planar electrodes on which an asymmetric RF bisinusoidal high voltage waveform, the Separation Voltage (SV), is applied. Due to different ion mobilities in the high field phase and the low field phase of the SV, ions are drifting in the direction perpendicular to the electrodes. This drift is countered by an DC Compensation Voltage (CV), which is scanned to get an differential mobility spectrum. 

There are multiple physical and chemical mechanisms leading to differential mobility of ions. The varying electric acceleration in the high and low field phase of the SV induce a variation of the effective ion temperature of the clustered ions. Therefore, the average cluster size is reduced in the high field phase, since high effective temperature leads to declustering. If the clusters can reform / regrow in the low field phase, this leads to changing average sizes of the clusters between high and low field phases, which in turn leads to different mobility in the high and low field phases. 

This application is intended to simulate compensation voltage due to chemically induced differential mobility of clustered small ions in a simplified way. This application is not a full trajectory simulation. Instead, the temperature and field dependent reactions of the individual species are simulated with RS. The drift of an ion in the separation dimension in a time step of length :math:`dt` is approximated simply by multiplying the ion mobility  :math:`K` with the time step length and the field strength :math:`E`:

    .. math:: 

        dz = K \cdot E(t) \cdot dt

The simulated compensation voltage (CV) can be found by varying the CV so that no net motion of the ions in the separation axis remains. 

The chemical reaction configuration is defined in an :ref:`RS configuration file <usersguide-rs-configurations>`. Currently, the mobility of the individual chemical species is assumed to be static and thus **not** dependent on the field or the effective ion temperature. The ions are initialized at random positions in a variable start zone. 


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


``reaction_configuration`` : file path 
    Path to a RS configuration file, defining the chemical reaction system for the simulation. 

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

``start_width_x_mm`` : float 
    Size (in mm) of the ion start zone in ``x`` direction, which is the long axis of the DMS device, parallel to the planar electrodes. 

``start_width_y_mm`` : float
    Size (in mm) of the ion start zone in ``y`` direction.

``start_width_z_mm`` : float
    Size (in mm) of the ion start zone in ``z`` direction, which is the separation axis. 

``cv_mode`` : Keyword: [``static``, ``auto``]
    Compensation voltage mode: 

    * ``static``: A static CV, defined by ``cv_Vmm-1`` is applied. 
    * ``auto``: Automatic mode to find an optimized CV for maximum transmission through the DMS. The CV is adapted after every SV cycle to minimize the total drift of the simulated particle ensemble in the separation (``z``) axis. The CV is written to a csv file ``<resultbasename>_cv.csv`` after every SV cycle in this mode. ``cv_Vmm-1`` defines the initial CV value in this mode.

``cv_relaxation_parameter`` : 
    Dampening parameter for the minimization process when ``cv_mode`` is ``auto``.

``sv_mode`` : 
    Separation waveform mode:

    * ``bi_sin``: "Classical" DMS wave form: A superposition of two sine waves, one with the base SV frequency and one with the double frequency, resulting in in a high and low field phase. 
    * ``square``: Ideal square waveform with 1/2 ratio between low and high field phase.
    * ``clipped_sin``: A clipped sine SV waveform, with half an sine wave as high field phase and a clipped phase as low field phase.


``sv_Vmm-1`` : 
    Separation voltage in V per mm.

``cv_Vmm-1`` : 
    Compensation voltage in V per mm. If ``cv_mode`` is ``auto``, initial CV voltage for the minimization process. 


``background_temperature_K`` : float 
    Background gas temperature

``background_pressure_Pa`` : float 
    Isotropic pressure of the neutral background gas in Pascal.
