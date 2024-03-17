project(test_Integration)

set(SOURCE_FILES
        test_main.cpp
        test_verletIntegrator.cpp
        test_parallelVerletIntegrator.cpp
        test_parallelRK4Integrator.cpp
        test_velocityIntegrator.cpp)

if(USE_FMM_3D)
    set (SOURCE_FILES ${SOURCE_FILES} test_fmmVerletIntegrator.cpp)
endif()

set(TEST_FILE_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields)
set(TEST_FILES
        ${TEST_FILE_FOLDER}/RS_verlet_test.conf
        )
file(COPY ${TEST_FILES} DESTINATION .)

add_executable(test_Integration ${SOURCE_FILES})

target_include_directories(test_Integration PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_Integration core integration file_io spacecharge particlesimulation collisionmodels ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

add_test(NAME test_Integration COMMAND test_Integration)