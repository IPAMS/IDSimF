
========================================================
Welcome to IDSimF, the Ion Dynamics Simulation Framework
========================================================

IDSimF, the *Ion Dynamics Simulation Framework*, is an open source framework for the simulation of non relativistic dynamics of molecular ions. The primary application of IDSimF is in the domain of mass spectrometry and ion mobility spectrometry.


Introduction and Overview
=========================

IDSimF is a collection of modules and programs which allows to simulate the trajectories of molecular ions at non relativistic conditions. It provides individual solver programs for individual simulation tasks and a modular framework for the quick development of new solvers. IDSimF was inspired by other simulation frameworks of a similar structure, particularly the open source fluid dynamics framework `OpenFOAM <https://openfoam.org/>`_. IDSimF is developed primarily in C++. Currently, the main contributor to IDSimF is the `Institute for Pure and Applied Mass Spectrometry (IPAMS) <https://www.ipams.uni-wuppertal.de/en/>`_ and the `Physical and Theroetical Chemistry <https://www.ptc.uni-wuppertal.de/>`_ of the [University of Wuppertal](https://www.uni-wuppertal.de/) . 

The primary application of IDSimF is the simulation of the dynamics of molecular ions in the domain of Mass Spectrometry (MS), Ion Mobility Spectrometry (IMS) and similar methods.  

--------
Features
--------

IDSimF currently provides as main features: 

* Charged particle trajectory integration in electric fields
* Import of electrostatics solutions and electric fields from other codes (e.g. `Comsol Multiphysics <https://www.comsol.com/>`_ and `SIMION <https://simion.com/>`_)
* Electric particle - particle interaction (space charge)
* Interaction of ions with neutral background gas
* Particle based chemical kinetics / simulation of reactive ions

------------------
What IDSimF is not
------------------

IDSimF is developed according to the requirements in molecular mass spectrometry. It is therefore *not at all* a full featured plasma code: 

* Currently only electrostatics are considered, no magnetic fields or magnetic couplings between particles are calculated. However, simple magnetic interactions (e.g. for ion cyclotron resonance simulations) can be added by explicit force terms in the trajectory integration.
* Oscillating / RF-Fields are simulated with an electrostatic approximation. Wave phenomena / RF radiation effects are not simulated. 

--------------------
Application Examples
--------------------

IDSimF is used for many different simulation problems by IPAMS. Examples of the application of IDSimF are: 

HiKE-IMS
--------

Conference contribution: `Simulation of Cluster Dynamics in High Kinetic Energy IMS (HiKE-­IMS) <https://www.ptc.uni-wuppertal.de/fileadmin/Chemie/Physikalische_Chemie/documents/published_docs/conferences/ASMS2019/ASMS2019_ThP299.pdf>`_ 

IDSimF was used to simulate the transport and the chemical kinetics of highly reactive proton bound water clusters through the drift tube of a High Kinetic Energy Ion Mobility Spectrometer (HiKE-IMS). The peak positions and the peak widths of the water cluster system depend largely on the reaction dynamics of the cluster species and are reproduced by the IDSimF simulation. 

----

Differential Mobility Spectrometry
----------------------------------

Conference contributions: 
`Chemical kinetic and ion transport simulations: Temperature dependence of ion mobility and its impact on cluster equilibria <https://www.ptc.uni-wuppertal.de/fileadmin/Chemie/Physikalische_Chemie/documents/published_docs/conferences/ASMS2018/ASMS2018_WP467_ChemicalKineticAndIon_TransportSimulations.pdf>`_ 

`Chemical Kinetics and Ion Transport Simulations: Cluster Dynamics in
Differential Ion Mobility Spectrometry <https://www.ptc.uni-wuppertal.de/fileadmin/Chemie/Physikalische_Chemie/documents/published_docs/conferences/ASMS2019/ASMS2019_TP521.pdf>`_ 

IDSimF was used to simulate the chemical reactions and the ion transport of proton bound water cluster ions in an Differential Ion Mobilty Spectrometry (DMS) cell of a commercial mass spectrometer. The IDSimF model consideres the interactions of the ions with the background gas in the DMS separation region at atmospheric pressure and the chemical reaction dynamics of the ion ensemble in the oscillating sparation field of the DMS. 

----

Space Charge Effects in Ion Trap 
--------------------------------

Conference contribution: `Evaluation of Space Charge Effects in Scanning- vs. Fourier Transform (FT)-Quadrupole Ion Traps (QITs) <https://www.ptc.uni-wuppertal.de/fileadmin/Chemie/Physikalische_Chemie/documents/published_docs/conferences/ASMS2018/ASMS2018_ThP520_EvaluationOfSpaceChargeEffectsInQITs.pdf>`_ 

IDSimF was used to investigate space charge effects in scanning and FT-Quadrupole Ion Trap (QIT) systems. The numerical modeling of the ion motion in the ion trap field and the ion detection clearly shows the space charge limits of the actual ion trap devices and the robustness of the different detection techniques. 


Documentation
=============

.. toctree::
   :maxdepth: 1
   :caption: Installation and Compilation

   installation/prerequisites_compilation
   installation/testing_installation

.. toctree::
    :maxdepth: 1
    :caption: User Guide    

    usersguide/general_structure
    usersguide/applications
    usersguide/rs
    usersguide/benchmarks

.. toctree::
   :maxdepth: 1
   :caption: Modules

   modules/module_overview
   modules/core
   modules/btree
   modules/particlesimulation
   modules/collisionmodels
   modules/rs


==================
Indices and tables
==================

   :ref:`genindex`
