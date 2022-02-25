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
        Integration_parallelVerletIntegrator.cpp)

if(USE_FMM_3D)
    set (SOURCE_FILES ${SOURCE_FILES}
            Integration_fmm3dIntegrator.cpp
            Integration_fmm3dIntegrator.hpp)
endif()

add_library(integration STATIC ${SOURCE_FILES})
target_include_directories(integration PUBLIC .)

target_link_libraries(integration core spacecharge collisionmodels file_io)
target_link_libraries(integration spdlog::spdlog)