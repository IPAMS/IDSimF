cmake_minimum_required(VERSION 3.5)
project(ionSimulation_cpp)

set(SOURCE_FILES
        simionPA_runtime_benchmark.cpp)

set(TEST_FILE_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields)
set(TEST_FILES
        ${TEST_FILE_FOLDER}/simion_test_planar_2d.pa
        ${TEST_FILE_FOLDER}/simion_test_planar_3d_mirrored.pa
        ${TEST_FILE_FOLDER}/simion_test_cylindrical_mirrored.pa
)
file(COPY ${TEST_FILES} DESTINATION .)

include_directories(${CMAKE_SOURCE_DIR}/modules/BTree)
include_directories(${CMAKE_SOURCE_DIR}/modules/particleSimulation)

add_executable(benchmark_simionpa_runtime ${SOURCE_FILES})
target_link_libraries(benchmark_simionpa_runtime core particlesimulation apputils)