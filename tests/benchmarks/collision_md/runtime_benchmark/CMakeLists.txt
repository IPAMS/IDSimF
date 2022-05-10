cmake_minimum_required(VERSION 3.5)
project(ionSimulation_cpp)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES
        md_runtime_benchmark.cpp)

include_directories(${CMAKE_SOURCE_DIR}/modules/BTree)
include_directories(${CMAKE_SOURCE_DIR}/modules/backgroundInteraction/collisionModels)
include_directories(${CMAKE_SOURCE_DIR}/modules/RS)
include_directories(${CMAKE_SOURCE_DIR}/modules/FileIO)

set(TEST_FILE_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields)
set(TEST_FILES
        ${TEST_FILE_FOLDER}/test_molecularstructure_reader.csv )
file(COPY ${TEST_FILES} DESTINATION .)

add_executable(benchmark_collision_md_runtime ${SOURCE_FILES})
target_link_libraries(benchmark_collision_md_runtime core collisionmodels apputils)