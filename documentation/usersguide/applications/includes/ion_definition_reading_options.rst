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

``n_ions`` : Vector of integers
    Number of ions of the ionic species defined by the the masses defined in ``ion_masses``. 

``ion_masses`` : Vector of float 
    Ion masses in amu. 

``ion_charges`` : Vector of float
    Ion charges in elementary charges.     

``ion_collision_gas_diameters_angstrom`` : Vector of float
    Effective hard sphere collision diameters of the ionic species in angström. 

``ion_time_of_birth_range_s`` : float
    Time range in which ions are generated, in seconds. The specified number of ions are generated uniformly in this time range.

Ion start configuration
.......................

The initial positions of the simulated ions can be a cubic box with faces parallel to the spatial axes or a rotatable and translatable cylinder. The position of the ion start zone is defined by ``ion_start_base_position_m``. This partameter defines the center of a box start zone or the position of the center of the bottom face of a cylinder start zone. 

``ion_start_geometry`` : Keyword:[``box``, ``cylinder``]
    Sets the ion start geometry.

    ``box`` : Ion start zone is a box
        The ion start zone is a cubic box with faces parallel to the spatial axes, randomly filled with particles. It is defined by

        ``ion_start_box_size_m``: Vector of three floats
            Size of the box in x, y, z direction

        ``ion_start_base_position_m``: Vector of three floats
            Position of the box center

    ``cylinder`` : Ion start zone is a cylinder
        The ion start zone is a rotatable and translatable cylinder. The cylinder is defined by

        ``ion_start_radius_m``: float
            Radius of the cylinder in m. 

        ``ion_start_length_m`` : float
            Distance from bottom face of the cylinder to the top face of the cylinder (cylinder length). 

        ``ion_start_cylinder_normal_vector``: Vector of three floats
            Normal vector of the bottom face of the cylinder, defining the rotation of the cylinder. The normal vector is also the directon of the axis of rotational symmetry of the cylinder (e.g. [1,0,0] means the cylinder axis is parallel to the x-axis).

        ``ion_start_base_position_m``: Vector of three floats
            Position of the center of the bottom face of the cylinder




