project(IDSimF_fileIO)

set(SOURCE_FILES
        PSim_averageChargePositionWriter.hpp
        PSim_averageChargePositionWriter.cpp
        PSim_simpleVTKwriter.hpp
        PSim_simpleVTKwriter.cpp
        PSim_trajectoryExplorerJSONwriter.hpp
        PSim_trajectoryExplorerJSONwriter.cpp
        PSim_ionCloudReader.cpp
        PSim_ionCloudReader.hpp
        PSim_scalar_writer.cpp
        PSim_scalar_writer.hpp
        PSim_idealizedQitFFTWriter.cpp
        PSim_idealizedQitFFTWriter.hpp
        PSim_inductionCurrentWriter.cpp
        PSim_inductionCurrentWriter.hpp
        PSim_HDF5Reader.hpp
        PSim_HDF5Reader.tpp
        PSim_HDF5Reader.cpp
        PSim_trajectoryHDF5Writer.hpp
        PSim_trajectoryHDF5Writer.cpp
        PSim_trajectoryWriterDefs.hpp)

add_library(file_io STATIC ${SOURCE_FILES})
target_include_directories(file_io PUBLIC .)

target_link_libraries(file_io core btree particlesimulation)
target_link_libraries(file_io ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES} spdlog::spdlog)