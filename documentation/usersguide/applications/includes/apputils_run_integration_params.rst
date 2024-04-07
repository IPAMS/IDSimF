``integrator_mode`` : Keyword:[``verlet``, ``parallel_verlet``, ``full_sum_verlet``, ``RK4``, ``full_sum_RK4``], optional: [``FMM3D_verlet``, ``ExaFMM_verlet``]
    Selects the trajectory integrator

    * ``verlet``: Serial (non parallelized) Velocity Verlet integrator with Barnes Hut Tree for space charge calculation
    * ``parallel_verlet``: Parallelized Velocity Verlet integrator with Barnes Hut Tree for space charge calculation
    * ``full_sum_verlet``: Parallelized Velocity Verlet integrator with full sum between all particles for space charge calculation 
    * ``RK4``: Parallelized Runge-Kutta 4 integrator with Barnes Hut Tree for space charge calculation
    * ``full_sum_RK4``: Parallelized Runge-Kutta 4 integrator full sum between all particles for space charge calculation 

With activated FMM libraries (Exa-FMM and or FMM3D) two FMM integrators are also available: 

    * ``FMM3D_verlet``: Velocity Verlet with Fast Multipole space charge calculation with FMM3D library 
    * ``ExaFMM_verlet``: Velocity Verlet with Fast Multipole space charge calculation with Exa-FMM library 
