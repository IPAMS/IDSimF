project(IDSimF_fileIO)

set(SOURCE_FILES
        FileIO_averageChargePositionWriter.hpp
        FileIO_averageChargePositionWriter.cpp
        FileIO_simpleVTKwriter.hpp
        FileIO_simpleVTKwriter.cpp
        FileIO_trajectoryExplorerJSONwriter.hpp
        FileIO_trajectoryExplorerJSONwriter.cpp
        FileIO_ionCloudReader.cpp
        FileIO_ionCloudReader.hpp
        FileIO_scalar_writer.cpp
        FileIO_scalar_writer.hpp
        FileIO_idealizedQitFFTWriter.cpp
        FileIO_idealizedQitFFTWriter.hpp
        FileIO_inductionCurrentWriter.cpp
        FileIO_inductionCurrentWriter.hpp
        FileIO_HDF5Reader.hpp
        FileIO_HDF5Reader.tpp
        FileIO_HDF5Reader.cpp
        FileIO_trajectoryHDF5Writer.hpp
        FileIO_trajectoryHDF5Writer.cpp
        FileIO_trajectoryWriterDefs.hpp
        FileIO_MolecularStructureReader.hpp
        FileIO_MolecularStructureReader.cpp)

add_library(file_io STATIC ${SOURCE_FILES})
target_include_directories(file_io PUBLIC .)

target_link_libraries(file_io core spacecharge particlesimulation)
target_link_libraries(file_io ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES} spdlog::spdlog)