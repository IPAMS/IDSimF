........
Exafmm-t
........

Clone / download `exafmm-t <https://exafmm.github.io/exafmm-t>`_ into a local directory and follow the `installation guide <https://exafmm.github.io/exafmm-t/compile.html#install-exafmm-t>`_. If all dependencies are available, configuration is done by changing into the exafmm-t directory and using autotools: 

.. code-block:: console
    
    cd exafmm-t
    ./configure

Then building and testing the build is done with: 

.. code-block:: console
    
    make check

After building, exafmm-t can be used in an IDSimF build by providing the path into the exafmm-t directory as ``EXAFMMT_PATH`` parameter. For example, if exafmm-t was cloned to ``/usr/libs/exafmm-t``, the additional cmake option is:

.. code-block:: console

    -DEXAFMMT_PATH=/usr/libs/exafmm-t