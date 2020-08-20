.. _usersguide-rs:

====================================
RS - Reaction Simulation Users Guide
====================================

The Reaction Simulation (RS) module simulates the chemical reaction kinetics of a simulated particle ensemble with a straight forward Monte Carlo approach: In every time step of the simulation, for every particle of the simulated particle ensemble, the local reaction probability of the individual particle is calculated. Then it is decided randomly, weighted by the calculated reaction probability, if a reaction of the particle takes place. Currently, this process is performed serially for the possible reactions of a particle, which induces an increasing statistical bias if the reaction probability of individual reactions becomes comparably high. 

In RS terms, the *reaction system* is the set of chemical species and reactions with their characteristics. There are different types of chemical species and reactions available to model chemical systems. They are described in detail below. 

RS models chemical kinetics with an ensemble of simulated particles. However, RS is currently not able to simulate reactions between particles or to simulate reactions which changes the number of simulated particles, since particle-particle collisions are currently not considered in RS. Therefore RS is currently suited to simulate the reactions of diluted species in a large excess of their reaction partners. This is particularly the situation of reactive ions in a neutral background gas at sufficiently high background gas pressure. 

.. _usersguide-rs-configurations:

RS Configurations
=================
The reaction system is controlled by configuration files ("RS configurations"), which are plain text files with simple structured format.

RS configurations are segmented in a ``[SUBSTANCES]`` part, which describes the chemical species in the simulation and a ``[REACTIONS]`` part, which specifies the set of simulated reactions of the species. The two segments are lists of chemical species and reactions specifying their parameters, one entry per line. 

The following example illustrates the general structure of RS configurations: 

.. code::

    [SUBSTANCES]
    A isotropic 1e10
    Cl_1 discrete 10 1 3e-4 0.0
    Cl_2 discrete 20 1 2e-4 0.0

    [REACTIONS]
    Cl_1 + A => Cl_2  | static ; 1.007190e-28 #cl1_forward

It describes a simple reaction system, with three chemical species, named ``A``, ``Cl_1`` and ``Cl_2`` and a single reaction in which ``Cl_1`` reacts with ``A`` to ``Cl_2``. 


Chemical species: ``[SUBSTANCES]`` section
==========================================

The ``[SUBSTANCES]`` section specifies the chemical substances in the RS simulation. Every substance is specified by one line in the ``[SUBSTANCES]`` section. The general format of a substance line is

.. code:: none

    <Substance Name> <Substance Type> <List of substance type specific parameters>

``<Substance Name>`` is a name of the substance, which is also used in the definition of reactions as identifier for the substance. ``<Substance Type>`` is the identifier / name of one of the substance types described below. Note that the name must not contain spaces or other special characters but it is allowed to contain underscores and numbers. After the substance type identifier, a type specific list of parameters defining the characteristics of the substance follows. 

---------------
Substance Types
---------------

``discrete`` : Reactive substance, modeled as particle
    Reactive substance which is modeled as a simulated particle. 

    ``discrete`` substances are defined by: 

    .. code:: none

        <Substance Name> discrete <molecular mass> <charge> <mobility> <collison diameter>

    The parameters are: 

    * ``<molecular mass>``: Mass of a particle of this substance in amu. 
    * ``<mobility>``: Electric mobility of the particles of this substance at standard conditions (:math:`K_0`) in m^2/(V s).
    * ``<collision diameter>`` Effective collision diameter of a particle of this substance in m.
    

``isotropic`` :  Isotropic background substance
    An isotropic substance in the background, not modeled as particle. This substance type can serve as reaction partner to ``discrete`` substances. 

    ``isotropic`` substances are defined by: 

    .. code:: none

        <Substance Name> isotropic <concentration>

    with 

    * ``<concentration>``: Isotropic concentration of this substance. The concentration unit has to be compatible with the reaction rate constants defined in ``[REACTIONS]`` section. 



Reactions: ``[REACTIONS]`` section
==================================

The general format to describe a reaction in the ``[REACTIONS]`` section of an RS configuration is: 

.. code:: none

    <Reaction Stochiometry> | <Reaction Type> <Reaction Parameter List> #<Reaction Name>

the parts of the reaction definition are: 

``<Reaction Stochiometry>`` : String defining the chemical stochiometry of the reaction
    Every reaction in RS is considered an elementary reaction. Therefore, there are only directional reactions but no equilibrium reactions. Equilibria are modeled implicitly by a forward and a backward elementary reaction. 
    
    The stochiometry of the reaction is defined by a string of the format 

    .. code:: none 

        <Educt List> => <Product List>

    Educt and product lists are lists of chemical substances, identified by their names in the ``[SUBSTANCES]`` section. 
    
    They have the form 

    .. code:: none 

        <Quantifier> <Substance Name> + <Quantifier> <Substance Name> + ..

    * ``Quantifier`` is a integer number, quantifying the number of molecules of the following substance for the stochiometry of the reaction. The quantifier can be omitted, it is implicitly replaced by 1 then. 
    * ``Substance Name`` is the name of a chemical substance defined in the ``[SUBSTANCES]`` section. 

    Specifiying the same substance multiple times in a educt or product list is legal, thus

    .. code:: none 

        A + A + 2 A => 2 B 

    is equivalent to 

    .. code:: none

        4 A => B + B

``<Reaction Type>`` : Reaction type identifier
    The name of one of the reaction types described below 

``<Reaction Parameter List>`` : List parameters of the reaction
    List of numeric values, specifying the parameters of the reaction. The list is semicolon separated and has the form 

    .. code:: none 

        ; <Value> ; <Value> ...

``<Reaction Name>``: Name of the reaction 
    Name of the reaction to identify it. Similarly to the substance name, this name must not contain spaces or similar special characters. 

    
For example 

.. code::

    subst_1 + 2 A => subst_2  | static ; 1.0e-5 #subst1_forward

specifies a reaction of type ``static`` named ``subst1_forward`` in which one ``subst_1`` reacts with two ``A`` to one ``subst_2``. The sole parameter to the ``static`` reaction is the static reaction rate constant and has a value of :math:`1.0 \times 10^{-5}`.

--------------
Reaction Types
--------------

There are currently 5 different reaction types, described in detail below. 

.. note::

    As described in the introduction, RS is currently not able to simualte reactions between simulated particles. Therefore, all reaction types have to have only one substance of type ``discrete`` in their educts and products. 

.. note::

    Not all simulation applications support all reaction types. It depends on the individual simulation application which reaction types are applicable. 


``static`` : Reaction with static reaction rate constant
--------------------------------------------------------

    Reaction of a ``discrete`` substance with background reaction partners, with a static rate constant :math:`k`. The reaction probability :math:`p` of this static reaction is calcuated by a simple linearized approach: 

    .. math::

        p = k \cdot \prod_i c_i \cdot \text{d}t

    with the reaction rate constant :math:`k`, the concentrations of the non discrete (``isotropic``) edcucts :math:`c_i` and the time step length :math:`\text{d}t`.

    **Parameter list:** ``; <rate constant>``

    * ``<rate constant>``: Static rate constant :math:`k`



``static_thermalizing`` : Thermalizing reaction with static reaction rate constant
----------------------------------------------------------------------------------

    Reaction of a ``discrete`` substance with background reaction partners, with a static rate constant :math:`k` and thermalization of the ``discrete`` simulated particle. This reaction is basically the same as ``static``, but the simulated particle of the ``discrete`` product is thermalized by the reaction. This means that the velocity vector of the product particle is reinitialized with a random velocity drawn from the Maxwell-Boltzmann distribution during the reaction.

    This reaction type is primarily intendet to model *resonant charge transfer* reactions of the type: 

    .. math:: 

        \text{A}^+ + \text{A} \rightarrow \text{A} + \text{A}^+

    which transfers charge from a potentially electrically accelerated particle to a thermal particle. 

    The reaction probability :math:`p` of this static reaction is calcuated by the same simple linearized approach as in ``static``: 

    .. math::

        p = k \cdot \prod_i c_i \cdot \text{d}t

    with the reaction rate constant :math:`k`, the concentrations of the non discrete (``isotropic``) edcucts :math:`c_i` and the time step length :math:`\text{d}t`.

    **Parameter list:** ``; <rate constant>``

    * ``<rate constant>``: Static rate constant :math:`k`


``vanthoff`` : Reaction with reaction rate given by van't Hoff equation 
-----------------------------------------------------------------------

    Reaction of a ``discrete`` substance with background reaction partners, with a temperature dependent rate constant :math:`k` calculated from a chemical equilibrium with the Van't Hoff equation [wissdorf2013]_.

    We assume that one is interested in the reaction rate constant :math:`k_{\text{A}\rightarrow\text{B}}` of a temperature dependent reaction 

    .. math::

            \text{A} \rightarrow \text{B}

    which is part of a simple chemical equilibrium 

    .. math::

        \text{B} \rightleftharpoons \text{A}. 

    The equilibrium constant :math:`K` of that equilibrium is the ratio of the forward rate constant :math:`k_{\text{B}\rightarrow\text{A}}` and backward rate constant :math:`k_{\text{A}\rightarrow\text{B}}` of the reactions forming the equilibrium: 

    .. math::

        K = \frac{k_{\text{B}\rightarrow\text{A}}}{k_{\text{A}\rightarrow\text{B}}}
    
    The Van't Hoff equation / isobar gives the connection between chemical equilibrium and the temperature: 

    .. math::
        \frac{\text{d}\ln{K}}{\text{d}T} = \frac{\Delta_r H^0(T)}{RT^2}

    with the molar standard reaction enthalpy as function of :math:`T`, :math:`\Delta_r H^0(T)`, and the universal gas constant :math:`R`.

    Assumption of :math:`\Delta_r H^0(T) = \text{const.}` in a temperature interval allows simple integration to get :math:`K_2` at a temperature :math:`T_2` from a known :math:`K_s` at temperature :math:`T_s`:

    .. math::

        \ln\left(\frac{K_2}{K_s}\right) = \frac{-\Delta_r H^0}{R}\left(\frac{1}{T_2} - \frac{1}{T_s}\right)

    it follows for :math:`K_2`: 

    .. math::

        K_2 = \exp\left(\frac{-\Delta_r H^0}{R}\left(\frac{1}{T_2} - \frac{1}{T_s}\right)\right) K_s

    The kinetic definition of :math:`K` therefore allows to calculate a the backward reaction rate constant, :math:`k_{\text{A}\rightarrow\text{B}}`, which is the temperature dependent reaction rate constant we were originally interested in, from a known forward reaction rate and the van't Hoff isobar:

    .. math::

        k_{\text{A}\rightarrow\text{B}} = \frac{k_{\text{B}\rightarrow\text{A}}}{K_2} = \frac{k_{\text{B}\rightarrow\text{A}}}{\exp{\left(\frac{-\Delta_r H^0}{R}\left(\frac{1}{T_2}-\frac{1}{T_s}\right)\right) K_s}}


    The reaction probability :math:`p` of the reaction is calcuated by the usual linearized approach: 

    .. math::

        p = k(T) \cdot \prod_i c_i \cdot \text{d}t

    with the temperature reaction rate constant :math:`k(T)`, the concentrations of the non discrete (``isotropic``) edcucts :math:`c_i` and the time step length :math:`\text{d}t`.

    **Parameter list:** ``; <H_r> ; <K_s> ; <k_backward>``

    * ``<H_r>``: Reaction enthalphy :math:`\Delta_r H^0`
    * ``<K_s>``: Equilibrium constant :math:`K_s` at standard temperature :math:`T_s = 298.15 \text{K}`
    * ``<k_backward>``: Known rate constant :math:`k_{\text{B}\rightarrow\text{A}}` of the equilibrium (which is the backward rate constant with respect to reaction we are interested in)

    References: 

        .. [wissdorf2013] Wissdorf, W., Seifert, L., Derpmann, V., Klee, S., Vautz, W., Benter, T.: Monte Carlo Simulation of Ion Trajectories of Reacting Chemical Systems: Mobility of Small Water Clusters in Ion Mobility Spectrometry. Journal of the American Society for Mass Spectrometry. 24, 632–641 (2013). https://doi.org/10.1007/s13361-012-0553-1


``vanthoff_field`` : Van't Hoff reaction with effective temperature calculated from electric field
--------------------------------------------------------------------------------------------------

    Reaction of a ``discrete`` substance with background reaction partners, with a field dependent rate constant :math:`k` calculated from a chemical equilibrium with the Van't Hoff equation with an effective temperature induced by an electric field on a charged particle. 

    This reaction type is similar to ``vanthoff``. The calculation of the reaction rate from a given temperature is the same as in ``vanthoff``, see the description above for details. In contrast to ``vanthoff``, where the reaction temperature :math:`T_2` is is a free parameter of a reaction event and is usually set by a simulation as the background temperature, this reaction type calculates an effective reaction temperature from the acceleration of a charged particle in an electric field. 

    If charged particles drifting in a gas filled region due to an moderate electric field with the field strengh :math:`E`, the particle velocities can be approximated by a Maxwell-Boltzmann distribution at an elevated temperature :math:`T_\text{eff}`. :math:`T_\text{eff}`, is given in first approximation by [mason1988]_ [viehland2012]_ [shvartsburg2008]_

    .. math::

        T_\text{eff} = T + \frac{M \left(KE\right)^2}{3 k}

    with the mass of the background gas particles :math:`M`, the electric mobility of the charged particles :math:`K`, the background gas temperature :math:`T` and the Boltzmann constant :math:`k`. 

    This reaction type uses this equation to calculate the effective temperature which is used for the estimation of the field dependent rate rate constant with the Van't Hoff reaction isobar as in ``vanthoff``. 

    The mobility :math:`K` in the equation above is calculated from the mobility at standard conditions :math:`K_0` of the ``discrete`` substance in the educts by

    .. math::

        K = K_0 \frac{p_0}{p} \frac{T}{T_0}

    with the local pressure :math:`p` and the temperature :math:`T` and standard pressure :math:`p_0 = 101325 \text{Pa}` and the standard temperature :math:`T_0 = 298.15 \text{K}`. 

    **Parameter list:** ``; <H_r> ; <K_s> ; <k_backward> ; <collison gas mass>``

        * ``<H_r>``: Reaction enthalphy :math:`\Delta_r H^0`
        * ``<K_s>``: Equilibrium constant :math:`K_s` at standard temperature :math:`T_s = 298.15 \text{K}`
        * ``<k_backward>``: Known rate constant :math:`k_{\text{B}\rightarrow\text{A}}` of the equilibrium (which is the backward rate constant with respect to reaction we are interested in)
        * ``<collison gas mass>``: Mass of the background gas particles in amu. 

    References: 

        .. [mason1988] Mason, E.A., McDaniel, E.W.: Transport Properties of Ions in Gases. Wiley-Interscience, New York (1988)
        .. [viehland2012] Viehland, L.A., Siems, W.F.: Uniform moment theory for charged particle motion in gases. Journal of the American Society for Mass Spectrometry. 23, 1841–54 (2012). https://doi.org/10.1007/s13361-012-0450-7
        .. [shvartsburg2008] Shvartsburg, A.A.: Differential Ion Mobility Spectrometry: Nonlinear Ion Transport and Fundamentals of FAIMS. CRC Press, Boca Raton (2008)
        

``simple_step`` : Collision based reaction with step like activation
--------------------------------------------------------------------

    A collision based reaction, which uses a step function at an activation energy :math:`E_A` to determine the reaction probability. 

    The total collison energy :math:`E_k` in a hard sphere collison of two particles with masses :math:`m_1` and :math:`m_2` colliding with the relative velocity :math:`v_{\text{rel}}` is given by

    .. math::

        E_k = \frac{1}{2} \frac{m_1 m_2}{m_1 + m_2} v_{\text{rel}}^2


    This reaction type assumes that if :math:`E_k` is higher than an activation energy :math:`E_A` (:math:`E_k > E_A`) the reaction takes place, if it is below the activation energy (:math:`E_k \leq E_A`) it never takes place. 

    **Parameter list:** ``; <activation_energy>``
        * ``<activation_energy>``: Activation energy :math:`E_A` in eV. 