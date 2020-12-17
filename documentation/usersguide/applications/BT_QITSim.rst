.. _application-BT-QITSim:

=========
BT-QITSim
=========

Model theory / description 
==========================

Ion trajectory simulation of an idealized quadrupole ion trap (QIT) considering space charge effects and hard sphere collisions with neutral background gas particles. 

.. _field-definition:

----------------
Field Definition
----------------

The electric field of the ion trap is defined by idealized, analytical equations. 
A quadrupole ion trap consists of two cap electrodes and a ring electrode, with hyperbolically inner surfaces which form an rotational symmetric quadrupolar field. The distance between the trap center and the cap electrodes is denoted as :math:`z_0`, the ring electrode radius :math:`r_0`.

Most commercial QITs have grounded end caps and apply DC and RF voltages to the ring electrode. In this case, the total electric field :math:`\Phi_{\text{total}}`, including non quadrupolar, higher field orders can be described by a series expansion [Eades1993]_:

.. math::
    \Phi_{\text{total}} = \Phi_0 + \Phi_0 \left( A_1 \Phi_1 + A_2 \Phi_2 + A_3 \Phi_3 + A_4 \Phi_4 + \dots \right)

Given the ideal geometric relationship for an ideal quadrupole ion trap :math:`r_0^2 = 2 z_0^2`, the electric potential in the center of the trap :math:`\Phi_0` is given by: 

.. math::
    \Phi_0 = \frac{U+V \cos{(\Omega t)}}{2}


with :math:`U=` DC potential on the cap electrodes, :math:`V=` RF potential between ring and caps, :math:`\Omega=` angular frequency of the RF potential, :math:`t=` time.

The linear electric potential between the caps is given by: 

.. math::
    \Phi_1 = \frac{z}{r_0}

with the Position of the ion :math:`z` in :math:`z`-Direction. 

The quadrupolar potential is given by: 

.. math:: 

    \Phi_2 = \frac{r^2 - 2 z^2}{r_0^2}

with the Position of the ion :math:`r` in :math:`r`-Direction.

Higher orders are given by: 

.. math:: 
    \begin{align}
    \Phi_3 &= \frac{2 z^3 - 3 z r^2}{r_0^3} \\
    \Phi_4 &= \frac{8 z^4 - 24 z^2 r^2 + 3r^4}{r_0^4} \\
    \Phi_5 &= \frac{8 z^5 - 40 z^3 r^2 + 15 z r^4}{r_0^5} \\
    \Phi_6 &= \frac{16 z^6 - 120 z^4 r^2 + 90 z^2 r^4 - 5r^6}{r_0^6} \\
    \end{align}

.. [Eades1993] Eades, D.M., Johnson, J. V, Yost, R.A.: Resonance in a quadrupole ion trap. J. Am. Soc. Mass Spectrom. 4, 917–929 (1993).

------------
FT detection
------------
There are different detection modes in quadrupole ion traps. Most commercial instruments eject the trapped ions by an resonant excitation with an additional bipolar RF potential on the cap electrodes and ramping the trap RF potential. 

Alternatively, the trapped ions can be excited with an excitation signal on the cap electrodes, which induces an oscillation of the ions around the trap center. This oscillation can be detected by the induced mirror charge on the cap electrodes. A mass spectrum can be calculated from the recorded mirror charge signal ("transient") by a Fourier transformation. 

The simulation app supports the simulation of this FT detection mode.

Simulation configuration description
====================================

``integrator_mode`` : Keyword:[``verlet``, ``parallel_verlet``]
    Selects the trajectory integrator

    * ``verlet``: Serial (non parallelized) verlet integrator
    * ``parallel_verlet``: Parallelized verlet integrator

``sim_time_steps`` : integer
    Number of simulation time steps

``trajectory_write_interval`` : integer
    Interval, in time steps, between writes to the trajectory result file.

``fft_write_interval`` : integer 
    Interval, in time steps, between samples for the FFT result file, which records a simulated transient for Fourier transformation. 

``fft_write_mode`` : Keyword:[``unresolved``, ``mass_resolved``]
    Selects the mode in which the FFT result file is written. The FFT result records the average velocity of the simulated particle ensemble in the ``z`` direction as an approximate induced mirror charge on the cap electrodes. 

    * ``unresolved``: The FFT transient signal is calculated from the whole simulated particle ensemble. This is the signal which would be detectable in a physical experiment. 
    * ``mass_resolved``: The FFT transient signal is recorded for every ionic species, with distinct molecular mass, separately. This allows to see individual transients for the species. 

``dt`` : float
    Time step length in seconds 

``space_charge_factor`` : float
    Multiplication factor for particle-particle interaction (space charge).

``max_ion_radius`` : float
    Maximum radius a simulated particle can have before it is terminated in the simulation, in m. 


.. include:: includes/ion_definition_reading_options.rst


----------------------------
Background gas configuration
----------------------------


``background_pressure_Pa`` : float
    Pressure of the neutral background as in Pa. 

``background_temperature_K`` : float
    Temperature of the background gas in K. 

``collision_gas_mass_amu`` : float
    Molecular mass of the particles of the background gas in amu.

``collision_gas_diameter_angstrom`` : float
    Effective collision diameter of the particles of the background gas in angström.

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

``field_mode`` : Keyword:[``basic``, ``higher_orders``]
    Selects model for the RF trap field

    ``basic`` : Pure quadrupolar trap field
        An ideal quadrupolar field is assumed as RF trap field. 

    ``higher_orders`` : Trap field with higher field orders
        Besides the quadrupolar field, hexapolar and octapolar field orders as defined in :ref:`field-definition` are considered too. The relative amount of higher field orders is defined by: 
        
        ``field_higher_orders_coeffs`` : Vector of two floats
            The relative amount of hexapolar and octapolar field order in the RF trapping field. 

``rf_waveform_csv_file`` : File path 
    File path to an file with an sampled RF waveform. If this parameter is set, the waveform of the trap RF is read from the specified file, which contains one sample per simulation timestep. The waveform is looped if the file contains fewer samples than ``sim_time_steps``. If this parameter is not set, a :math:`\cos{\left(\omega t \right)}` is used as trap field waveform. 

    Note that ``f_rf`` has no effect if a sampled RF waveform is used. This file path is relative to the simulation run configuration file. 


Trap field RF voltage
---------------------

The RF trap field voltage can be static or can be ramped during the simulation. 

A static trap field voltage is set by ``V_rf``:

``V_rf`` : float
    Static ground to peak trap field amplitude in volt. 

If ``V_rf_start`` and ``V_rf_end`` are set, the trap field voltage is ramped linearly during the simulation from ``V_rf_start`` to ``V_rf_end``: 

``V_rf_start`` : float 
    Ground to peak trap field amplitude in volt at the begin of the trap field voltage ramp.

``V_rf_end`` : float
    Ground to peak trap field amplitude in volt at the end of the trap field voltage ramp.


Ion excitation field
--------------------

The ion excitation field is applied bipolar on the cap electrodes. There are two modes of ion excitation: 

* Pulsed excitation with a rectangular excitation pulse of defined length and amplitude applied at the begin of the simulation run.
* Excitation with a given sampled waveform read from a waveform file. 

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