cmake_minimum_required(VERSION 3.5)
project(ionSimulation_cpp)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES
        hs_runtime_benchmark.cpp)

include_directories(${CMAKE_SOURCE_DIR}/modules/BTree)
include_directories(${CMAKE_SOURCE_DIR}/modules/backgroundInteraction/collisionModels)
include_directories(${CMAKE_SOURCE_DIR}/modules/RS)

add_executable(benchmark_collision_hardsphere_runtime ${SOURCE_FILES})
target_link_libraries(benchmark_collision_hardsphere_runtime core collisionmodels apputils)