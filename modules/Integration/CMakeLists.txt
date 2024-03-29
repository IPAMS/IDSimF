project(IDSimF_integration)

set(SOURCE_FILES
        Integration_generic.hpp
        Integration_abstractTimeIntegrator.hpp
        Integration_abstractTimeIntegrator.cpp
        Integration_velocityIntegrator.hpp
        Integration_velocityIntegrator.cpp
        Integration_verletIntegrator.hpp
        Integration_verletIntegrator.cpp
        Integration_parallelVerletIntegrator.hpp
        Integration_parallelVerletIntegrator.cpp
        Integration_parallelRK4Integrator.hpp
        Integration_parallelRK4Integrator.cpp
        Integration_fullSumRK4Integrator.hpp
        Integration_fullSumRK4Integrator.cpp)

if(USE_FMM_3D)
    set (SOURCE_FILES ${SOURCE_FILES}
            Integration_fmmIntegrator.cpp
            Integration_fmmIntegrator.hpp)
endif()

add_library(integration STATIC ${SOURCE_FILES})
target_include_directories(integration PUBLIC .)

target_link_libraries(integration core spacecharge collisionmodels file_io)
target_link_libraries(integration spdlog::spdlog)
if(OpenMP_CXX_FOUND)
    target_link_libraries(integration OpenMP::OpenMP_CXX)
endif()