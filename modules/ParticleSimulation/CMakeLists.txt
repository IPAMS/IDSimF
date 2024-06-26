project(IDSimF_particleSimulation)

set(SOURCE_FILES
        PSim_interpolatedField.hpp
        PSim_interpolatedField.cpp
        PSim_util.hpp
        PSim_util.cpp
        PSim_sampledWaveform.hpp
        PSim_sampledWaveform.cpp
        PSim_math.cpp
        PSim_math.hpp
        PSim_simionPotentialArray.cpp
        PSim_simionPotentialArray.hpp
        PSim_particleStartZone.hpp
        PSim_cylinderStartZone.cpp
        PSim_cylinderStartZone.hpp
        PSim_boxStartZone.cpp
        PSim_boxStartZone.hpp
        PSim_particleStartZone.cpp
        PSim_particleStartSplatTracker.cpp
        PSim_particleStartSplatTracker.hpp
        PSim_constants.hpp
        PSim_sampledFunction.cpp
        PSim_sampledFunction.hpp PSim_sphereStartZone.cpp PSim_sphereStartZone.cpp)

add_library(particlesimulation STATIC ${SOURCE_FILES})
target_include_directories(particlesimulation PUBLIC .)

target_link_libraries(particlesimulation core spacecharge file_io collisionmodels)
target_link_libraries(particlesimulation spdlog::spdlog)