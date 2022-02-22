cmake_minimum_required(VERSION 3.5)
project(ionSimulation_cpp)

set(SOURCE_FILES
        interpolatedField_runtime_benchmark.cpp)

set(TEST_FILE_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields)
set(TEST_FILES
        ${TEST_FILE_FOLDER}/test_linear_scalar_field_01.h5
        ${TEST_FILE_FOLDER}/test_linear_vector_field_01.h5
)
file(COPY ${TEST_FILES} DESTINATION .)

add_executable(benchmark_interpolated_field_runtime ${SOURCE_FILES})
target_link_libraries(benchmark_interpolated_field_runtime core file_io particlesimulation)