.....
FMM3D
.....

Clone / download `FMM3D <https://fmm3d.readthedocs.io/en/latest/index.html>`_ into a local directory and follow the install instructions `install instructions <https://fmm3d.readthedocs.io/en/latest/install.html#obtaining-fmm3d>`_. If all dependencies are available, configuration is done by running make in the FMM3D directory: 

.. code-block:: console
    
    cd FMM3D
    make lib


After building, FMM3D can be used in an IDSimF build by providing the path into the FMM3D directory as ``FMM_3D_PATH`` parameter. For example, if exafmm-t was cloned to ``/usr/libs/FMM3D``, the additional cmake option is:

.. code-block:: console

    -DFMM_3D_PATH=/usr/libs/FMM3D

