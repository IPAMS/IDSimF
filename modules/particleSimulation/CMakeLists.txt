project(IDSimF_particleSimulation)

set(SOURCE_FILES
        PSim_interpolatedField.hpp
        PSim_interpolatedField.cpp
        PSim_util.hpp
        PSim_util.cpp
        PSim_averageChargePositionWriter.hpp
        PSim_averageChargePositionWriter.cpp
        PSim_simpleVTKwriter.hpp
        PSim_simpleVTKwriter.cpp
        PSim_sampledWaveform.hpp
        PSim_sampledWaveform.cpp
        PSim_trajectoryExplorerJSONwriter.hpp
        PSim_trajectoryExplorerJSONwriter.cpp
        PSim_abstractTimeIntegrator.hpp
        PSim_abstractTimeIntegrator.cpp
        PSim_velocityIntegrator.hpp
        PSim_velocityIntegrator.cpp
        PSim_verletIntegrator.hpp
        PSim_verletIntegrator.cpp
        PSim_parallelVerletIntegrator.hpp
        PSim_parallelVerletIntegrator.cpp
        PSim_math.cpp
        PSim_math.hpp
        PSim_ionCloudReader.cpp
        PSim_ionCloudReader.hpp
        PSim_scalar_writer.cpp
        PSim_scalar_writer.hpp
        PSim_idealizedQitFFTWriter.cpp
        PSim_idealizedQitFFTWriter.hpp
        PSim_inductionCurrentWriter.cpp
        PSim_inductionCurrentWriter.hpp
        PSim_trajectoryHDF5Writer.hpp
        PSim_trajectoryHDF5Writer.cpp
        PSim_trajectoryWriterDefs.hpp
        PSim_simionPotentialArray.cpp
        PSim_simionPotentialArray.hpp
        PSim_HDF5Reader.cpp
        PSim_HDF5Reader.hpp
        PSim_particleStartZone.hpp
        PSim_cylinderStartZone.cpp
        PSim_cylinderStartZone.hpp
        PSim_boxStartZone.cpp
        PSim_boxStartZone.hpp PSim_particleStartZone.cpp PSim_particleStartSplatTracker.cpp PSim_particleStartSplatTracker.hpp PSim_constants.hpp)

add_library(particlesimulation STATIC ${SOURCE_FILES})
target_include_directories(particlesimulation PUBLIC .)

target_link_libraries(particlesimulation core btree collisionmodels)
target_link_libraries(particlesimulation ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES} spdlog::spdlog)