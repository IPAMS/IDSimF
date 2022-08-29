.. _modules-fileio:

===============
Module: File IO 
===============


File Writers
============

File writer classes are used to export data to persistent files. 


------------------------------------------
Primary Simulation Result Data File Writer
------------------------------------------

:cpp:class:`FileIO::TrajectoryHDF5Writer` writes trajectory HDF5 files, which is the current primary trajectory data export format of IDSimF.

.. doxygenclass:: FileIO::TrajectoryHDF5Writer
    :members:
    :undoc-members:


:cpp:class:`FileIO::Scalar_writer` writes tables of scalar values from simulations:

.. doxygenclass:: FileIO::Scalar_writer
    :members:
    :undoc-members:


-----------------------------
Additional Result File Writer
-----------------------------

Additional file writer provide additional export file formats. 

.. note:: 

    The additional file writer are currently not well maintained. 

.. doxygenclass:: FileIO::TrajectoryExplorerJSONwriter
    :members:
    :undoc-members:

.. doxygenclass:: FileIO::SimpleVTKwriter
    :members:
    :undoc-members:


------------------------------
Special Simulation File Writer
------------------------------

There are some file writers for special simulation requirements: 

.. doxygenclass:: FileIO::InductionCurrentWriter
    :members:
    :undoc-members:

.. doxygenclass:: FileIO::IdealizedQitFFTWriter
    :members:
    :undoc-members:

.. doxygenclass:: FileIO::AverageChargePositionWriter
    :members:
    :undoc-members:        


File Readers
============

File readers import data from persistent files 

:cpp:class:`FileIO::HDF5Reader` is a general reader for HDF5 files.

.. doxygenclass:: FileIO::HDF5Reader
    :members:
    :undoc-members:

.. doxygenclass:: FileIO::IonCloudReader
    :members:
    :undoc-members:

.. doxygenclass:: FileIO::MolecularStructureReader
    :members:
    :undoc-members:
