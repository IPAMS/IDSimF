project(IDSimF_particleSimulation)

set(SOURCE_FILES
        PSim_interpolatedField.cpp
        PSim_util.cpp
        PSim_sampledWaveform.cpp
        PSim_math.cpp
        PSim_simionPotentialArray.cpp
        PSim_cylinderStartZone.cpp
        PSim_boxStartZone.cpp
        PSim_particleStartZone.cpp
        PSim_particleStartSplatTracker.cpp
        PSim_sampledFunction.cpp
        PSim_sphereStartZone.cpp
        PSim_sphereStartZone.cpp)

add_library(particlesimulation STATIC ${SOURCE_FILES})
target_include_directories(particlesimulation PUBLIC .)

target_link_libraries(particlesimulation PRIVATE core spacecharge file_io collisionmodels spdlog::spdlog)