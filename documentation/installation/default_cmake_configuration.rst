IDSimF uses ``cmake`` as helping tool for configuration and compilation. ``cmake`` allows a so called "out of source build". This creates a separated "binary tree" in a separated build folder, in which the compilation of executable binaries takes place without interfering with the cloned sources. 

To do an out of source build, change into the cloned IDSimF folder and create a build folder, for example ``build`` in it and change into it: 

.. code-block::
    
    cd IDSimF
    mkdir build
    cd build

Basically ``cmake`` prepares a build tree in the current folder if it is called with an source folder as argument. Optional parameters for the build configuration are set with optional arguments the form ``-D<OPTION NAME>=<VALUE>``.

In the build folder, prepare build / binary tree with 

.. code-block::

    cmake .. -DCMAKE_BUILD_TYPE=Release

The build option mean the following: 

    + ``-DCMAKE_BUILD_TYPE=Release`` sets the build type to `Release` which means, that optimizations are switched on and debugging information is removed from the compiled binary. This results in a significantly faster binary than when building with debugging information switched on (``-DCMAKE_BUILD_TYPE=Debug``) which is required to analyze the compiled binary with debugging tools. 