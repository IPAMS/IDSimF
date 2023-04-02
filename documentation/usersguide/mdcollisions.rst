.. _usersguide-mdc:

==================================
Molecular Dynamics Collision Model 
==================================

The molecular dynamics (MD) model is an additional collision model, which allows to describe the interactions between ions and neutral background gas particles.
Fundamentally, the idea of this model is to solve Newton's equations of motion in a separate trajectory integration for the particle system defined by the underlying force field and to update the velocities accordingly.

.. _usersguide-mdc-principle:

MD Collision Principle
======================

There are two main principles underlying the basic molecular dynamics based collision approach:

* The collision partners are modeled as rigid bodys, meaning individual atoms are explicitly modeled by their relative position inside the molecule and by additional parameters, e.g. mass and charge. In the here discussed model, no vibrations or propagation of rotations are considered. All molecules are rotated randomly in the initalization step of the collision, but this rotation stays fixed over the whole integration process. 
* Both particle's trajectories are fully integrated by solving Newton's equations of motion.  

Newton's equations of motion are given by 

.. math::

    m \cdot \frac{d^2r}{dt^2} = \sum_i F_i \quad, 

where :math:`m` is the mass of the particle, :math:`r` the position of the particle and :math:`F_i` the individual force components acting on the particle. 
To be able to solve the equations of motions an appropiate force field has to be supplied. 
The force field used in IDSimF follows the ideas of MOBCAL [shvartsburg1996]_ and IMoS [larriba2013]_ by representing short-range repulsion through the commonly implemented 12-6 Lennard-Jones
potential and longe-range attraction by the ion-induced dipole moment [ieritano2019]_:

.. math::

    \begin{split}
    				V_{ij}(r_{ij}) &= V_{LJ, ij}(r_{ij}) + V_{ind., ij}(r_{ij}) \\
					     &= 4\epsilon_{ij} \left[ \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6} \right] + \\
					     &\frac{\alpha}{8\pi\epsilon_0}\left[ \left( z_ie\frac{r_{x,ij}}{r_{ij}^3}\right)^{2}  + \left( z_ie\frac{r_{y,ij}}{r_{ij}^3}\right)^{2} + 
					     \left( z_ie\frac{r_{z,ij}}{r_{ij}^3}\right)^{2}\right] \quad ,
	\end{split} 

where :math:`\epsilon_{ij}` and :math:`\sigma_{ij}` are the Lennard-Jones parameter for the well-depth and point of zero potential for the interaction between atom :math:`j` and atom :math:`i`, 
:math:`\epsilon_0` is the permittivity of free space, :math:`\alpha` is the polarizability volume of the background gas, :math:`e` is the elementary charge, :math:`z_i` is the charge of atom :math:`i` and 
:math:`\mathbf{r}_{ij} = (r_{x,ij}, r_{y,ij}, r_{z,ij})^\mathsf{T}` is the distance vector between atom :math:`i` and atom :math:`j`, with a magnitude of :math:`r_{ij}`.
Additionally, if the background gas is diatomic nitrogen the quadrupole potential based on the five charge model is also considered [ieritano2019]_

.. math::

    V_{QP, ij}(r_{ij}) = \frac{1}{8\pi\epsilon_0} \sum_{j=1}^3 \frac{z_i z_j e^2}{r_{ij}}\quad .


Before the particle's trajectories can be integrated, the collision itself has to be initalized. 
This process requires sampling of background gas particle position and velocity and is shortly outlined in three steps:

#. Construct the rigid body representation and place atoms based on their relative coordinates. The background gas velocity is estimated by the 1D Maxwell Boltzmann distribution. The system is shifted to the stationary ion frame and the ions center-of-mass is placed at the coordinate system origin.
#. The background gas particle is placed on a randomly drawn point :math:`p` on a sphere (uniformly distributed), with radius :math:`r_{spawn}` around the coordinate origin. 
#. The positions have to be pre-selected to save on compute time. This is done by only selecting position and velocity combinations, which would enter the collision sphere based on geometric considerations. As such, the angle :math:`\alpha` between the velocity vector :math:`v` and the vector from :math:`p` to the midpoint of the inner sphere as well as angle :math:`\theta` between the collision sphere and midpoint vector have to be compared. A position is sampled until :math:`\alpha \leq \theta`. 

After the initalization, the trajectory integration starts. There are three integration methods implemented: Leapfrog, Runge-Kutta 4 and Runge-Kutta-Fehlberg 45. Both the Leapfrog and Runge-Kutta 4 method are a constant step size integration method, but are generally not used in the program itself. Instead, the adpative stepsize method 
Runge-Kutta-Fehlberg is used, which allows for orders of magnitude faster calculation of the trajectories. 

.. _usersguide-mdc-configurations:

MD Configurations
=================

To be able to use the molecular dynamics collision model an additional input file, given in the .json file format, is necessary in which the molecular structures are defined. 
MD configuration files include different molecular structures, uniquely identifiable by their string identifier ``name``. For each molecular structure a diameter (``diameter``) is given and all its atoms are listed with the following parameters:

* ``type``: atome type 
* ``posx``: relative atom position in x (in Å) 
* ``posy``: relative atom position in y (in Å) 
* ``posz``: relative atom position in z (in Å)
* ``mass``: atom mass (in amu) 
* ``charge``: atom charge (in elementary charges) 
* ``partCharge``: atom partial charge (in elementary charges) 
* ``LJsigma``: Lennard-Jones potential :math:`\sigma` parameter (in Å) 
* ``LJeps``: Lennard-Jones potential :math:`\epsilon` parameter (in :math:`\frac{kJ}{mol}`) 

An example of a molecular configuration file is given below: 

.. code-block:: JSON

    [
        {
            "name":"Ar2",
            "diameter":4,
            "atoms":[
                {
                    "type":"Ar", 
                    "posx":0, 
                    "posy":1.85, 
                    "posz":0, 
                    "mass":39.948, 
                    "charge":0, 
                    "partCharge":0, 
                    "LJsigma":3.401, 
                    "LJeps":0.9777
                }, 
                {
                    "type":"Ar", 
                    "posx":0, 
                    "posy":-1.85, 
                    "posz":0, 
                    "mass":39.948, 
                    "charge":0, 
                    "partCharge":0, 
                    "LJsigma":3.401, 
                    "LJeps":0.9777
                }

            ]

        },
        {
            "name":"He",
            "diameter":2.89,
            "atoms":[
                {
                    "type":"He", 
                    "posx":0, 
                    "posy":0, 
                    "posz":0, 
                    "mass":4.003, 
                    "charge":0, 
                    "partCharge":0, 
                    "LJsigma":2.556, 
                    "LJeps":0.0837
                }
            ]

        }
    ]

In addition to the molecular structure file, input files from other applications in IDSimF (see :doc:`applications users guide <applications>`) have to be modified in order to call the
MD collision model. The user has to provide the following information (in accordance to the MD collision model interface :cpp:class:`CollisionModel::MDInteractionsModel`):

* the collision gas polarizability in m³  

* the collision gas identifier for the molecular structure

* the particle identifier for the molecular structure of the ion

* the path to the molecular structure .csv file

* the maximum integration time for the sub-integrator of the MD collision model in seconds

* the step size for the sub-integrator in seconds

* the collision radius scaling

* the scaling for angle :math:`\theta`

* the radius of the spawn sphere in m 

* a flag to possibly save the simulated trajectories of the MD collision model

* the distance in m after which the trajectory data should be saved  

* the time step after which the trajectory data should be saved.

.. note::
   Special care should be taken when setting the collision radius scaling and spawn sphere radius. It should be ensured that 
   the background gas particle is placed in an area of only marginal potential field, as otherwise the collision cannot 
   perserve energy conservation. 
   Also an appropiate colision radius scaling has to bet choosen, which allows to capture all long-range collision interactions. 
   The scaling is typically at least set to 2 if Helium is used as the background gas. 
   If diatomic nitrogen is used it is recommended to chose a scaling of at least 3 or higher.

References:
-----------

    .. [larriba2013] C. Larriba and C. J. Hogan, “Ion mobilities in diatomic gases: Measurement versus prediction with non-specular scattering models,” 
        The Journal of Physical Chemistry A, vol. 117, no. 19, pp. 3887–3901, 2013. PMID: 23488939.
    .. [ieritano2019] C. Ieritano et al., “A parallelized molecular collision cross section package with optimized accuracy and efficiency,” 
        Analyst, vol. 144, no. 5, pp. 1660-1670, 2019. doi: 10.1039/C8AN02150C.
    .. [shvartsburg1996] A. A. Shvartsburg and M. F. Jarrold, “An exact hard-spheres scattering model for the mobilities of polyatomic ions,” 
        Chemical Physics Letters, vol. 261, no. 1, pp. 86–91, 1996.
