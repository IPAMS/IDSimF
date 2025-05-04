project(IDSimF_fileIO)

set(SOURCE_FILES
        FileIO_averageChargePositionWriter.cpp
        FileIO_simpleVTKwriter.cpp
        FileIO_ionCloudReader.cpp
        FileIO_scalar_writer.cpp
        FileIO_idealizedQitFFTWriter.cpp
        FileIO_inductionCurrentWriter.cpp
        FileIO_HDF5Reader.cpp
        FileIO_trajectoryHDF5Writer.cpp
        FileIO_MolecularStructureReader.cpp
        FileIO_CSVReader.cpp)

add_library(file_io STATIC ${SOURCE_FILES})
target_include_directories(file_io PUBLIC .)

target_link_libraries(file_io PRIVATE core particlesimulation spacecharge collisionmodels
        ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES} jsoncpp spdlog::spdlog)