.. _application-BT-RS-IMSSim:

============
BT-RS-IMSSim
============

Simulates the trajectories of charged particles in a drift tube Ion Mobility Spectrometry (IMS) Device, including background gas interaction, ion chemistry and space charge. 
A drift tube IMS consists of a gas filled drift region with a uniform electric field present. Short pulses of ions are gated into this drift region and drift through the drift region due to the electric field. Since different ions have different electric mobilities, ions separate during the drift through the drift region. The drift times are measured to get a Mobility Spectrum. 

Ions start uniformly distributed in a configurable start zone. An ideally uniform electric field is assumed in the drift region. The interaction between background gas and the ions can be described with different collision models. Electric ion-ion interaction (space charge) can be modeled with a parallelized version of a Barnes-Hut tree. 

The modeled ions are chemically reactive. The temperature dependent reactions of the individual chemical species are simulated with RS. The reactions are defined in an :ref:`RS configuration file <usersguide-rs-configurations>`. 

Note that the simulated species, particularly the electric mobilities, are defined in the RS configuration. Therefore, to simulate chemically inert ions, an RS configuration without reactions between the defined species has to be defined. 


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

    .. code-block:: 

        [SUBSTANCES]
        Cl_1 discrete 19 1 3.57e-4  4.00000000e-10
        Cl_2 discrete 37 1 2.76e-4  5.17391304e-10
        Cl_3 discrete 55 1 2.35e-4  6.07659574e-10

    the ``n_ions`` vector ``[100, 50, 10]`` will initalize the simulation with 100 particles of ``Cl_1``, 50 of ``Cl_2`` and 10 of ``Cl_3``. 

``electric_field_mag_Vm-1`` : float
    Magnitude of the electric drift field in the drift region in V per m. The electric drift field is modeled as ideal uniform electric field in ``x`` direction. 

``start_width_yz_mm`` : float
    Size (in mm) of the ion start in ``y`` and ``z`` direction, perpendicular to the drift axis of the IMS. 

``start_width_x_mm`` : float
    Size (in mm) of the ion start zone in ``x`` direction, which is the drift axis of the IMS. The ion start zone begins at ``x=0``. 

``stop_position_x_mm`` : float
    Position along the drift axis of the IMS (``x`` direction) where the drifting ions terminate. 

``background_temperature_K`` : float
    Isotropic temperature of the background gas in K. 

-----------------------------------------------
Transport models and background gas interaction 
-----------------------------------------------

The simulation has different modes to model the ion migration through the drift region and different modes to model the background gas interaction which are suitable for different background gas pressure ranges. 

The transport model mode is controlled by the ``transport_model_type`` parameter: 

``transport_model_type`` : keyword [``btree_SDS``, ``btree_HS``, ``simple``, ``no_transport``]
    ``btree_SDS`` : Simulation with full trajectory integration, space charge and high pressure gas collision model (SDS)
        Simulation with full trajectory integration and activated space charge modeling with a parallelized Barnes-Hut method. The gas interaction is described by the Statistical Diffusion Simulation (SDS) gas interaction model, which is suitable for high background gas pressure. The assumption of SDS is essentially, that many ion-neutral collsions takes place between the time steps of the simulation. 

        The SDS collision model allows no background gas mixtures. Thus, the gas definition parameters described below, which are vectors of parameters, must have only one entry each with SDS.  

        .. note::
            The time step length is *not* adapted to the gas dynamic parameters of the gas interaction model and there is *no warning* if the assumptions of the SDS model are violated in a simulation. 

    ``btree_HS`` : Simulation with full trajectory integration, space charge and Hard Sphere (SDS) gas collision model
        Simulation with full trajectory integration and activated space charge modeling with a parallelized Barnes-Hut method. The gas interaction is described by the Hard Sphere (HS) gas collision model, which is suitable for low background gas pressure. 

        The HS model is able to model mixtures of gases, which are defined by the gas parameter vectors described below. 

        .. note::
            The time step length is *not* adapted to the gas dynamic parameters of the hard sphere model model. For a valid modeling, the time step length should be significantly shorter than the mean time between ion-neutral collisions.

    ``simple`` : Simple transport without gas interaction and space charge
        Simple transport mode without full trajectory integration. The ion migration distance :math:`dx` in a time step of length :math:`dt` in a field :math:`E` is calculated in this mode from the local ion mobility :math:`K_{\text{l}}` by

        .. math::

            dx = K_{\text{l}} \cdot E \cdot dt

        No diffusion or space charge effects are calculated in this mode. 

    ``no_transport`` : No transport modeling, chemical kinetics only. 
        No transport simulation takes place at all, only chemical reactions of the particle ensemble with background gas components are simulated. 

``background_partial_pressures_Pa`` : vector of float 
    Partial pressures of the individual components of the background gas mixture in Pascal. Note that with SDS background gas interaction model, only one background gas component is allowed. 

``collision_gas_masses_amu`` : vector of float
    Molecular masses of the particles of the background gas mixture components in amu. Note that with SDS background gas interaction model, only one background gas component is allowed. 

``collision_gas_diameters_angstrom`` : vector of float
    Effective collision diameters of the particles of the background gas components in Angström. Note that with SDS background gas interaction model, only one background gas component is allowed. 
