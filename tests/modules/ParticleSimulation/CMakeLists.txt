project(test_ParticleSimulation)

set(SOURCE_FILES
        test_main.cpp
        test_interpolatedField.cpp
        test_sampledWaveform.cpp
        test_math.cpp
        test_util.cpp
        test_simionPotentialArray.cpp
        test_startZones.cpp
        test_particleStartSplatTracker.cpp
        test_sampledFunction.cpp)

set(TEST_FILE_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields)
set(TEST_FILES
        ${TEST_FILE_FOLDER}/swift_test_sin.csv
        ${TEST_FILE_FOLDER}/swift_test_transient.csv
        ${TEST_FILE_FOLDER}/low_sample_waveform.csv
        ${TEST_FILE_FOLDER}/triangle_waveform.csv
        ${TEST_FILE_FOLDER}/function_linear.csv
        ${TEST_FILE_FOLDER}/function_quadratic.csv
        ${TEST_FILE_FOLDER}/simion_test_cylindrical.pa
        ${TEST_FILE_FOLDER}/simion_test_cylindrical_mirrored.pa
        ${TEST_FILE_FOLDER}/simion_test_planar_2d.pa
        ${TEST_FILE_FOLDER}/simion_test_planar_2d_scaled.pa
        ${TEST_FILE_FOLDER}/simion_test_planar_2d_nonUniformlyScaled.pa
        ${TEST_FILE_FOLDER}/simion_test_planar_3d.pa
        ${TEST_FILE_FOLDER}/simion_test_planar_3d_mirrored.pa
        ${TEST_FILE_FOLDER}/test_linear_vector_field_01.h5
        ${TEST_FILE_FOLDER}/test_linear_scalar_field_01.h5
        ${TEST_FILE_FOLDER}/quad_dev_flow_3d.h5
        )
file(COPY ${TEST_FILES} DESTINATION .)

add_executable(test_ParticleSimulation ${SOURCE_FILES})
target_include_directories(test_ParticleSimulation PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_ParticleSimulation core particlesimulation spacecharge integration file_io collisionmodels rs ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

add_test(NAME test_ParticleSimulation COMMAND test_ParticleSimulation)