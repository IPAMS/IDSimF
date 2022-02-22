project(IDSimF_integration)

set(SOURCE_FILES
        PSim_abstractTimeIntegrator.hpp
        PSim_abstractTimeIntegrator.cpp
        PSim_velocityIntegrator.hpp
        PSim_velocityIntegrator.cpp
        PSim_verletIntegrator.hpp
        PSim_verletIntegrator.cpp
        PSim_parallelVerletIntegrator.hpp
        PSim_parallelVerletIntegrator.cpp)

add_library(integration STATIC ${SOURCE_FILES})
target_include_directories(integration PUBLIC .)

target_link_libraries(integration core btree collisionmodels file_io)
target_link_libraries(integration spdlog::spdlog)