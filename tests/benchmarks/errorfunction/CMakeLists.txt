cmake_minimum_required(VERSION 3.5)
project(ionSimulation_cpp)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES
        main.cpp)

include_directories(${CMAKE_SOURCE_DIR}/modules/BTree)
include_directories(${CMAKE_SOURCE_DIR}/modules/backgroundInteraction/collisionModels)

add_executable(benchmark_errorfunction ${SOURCE_FILES})
target_link_libraries(benchmark_errorfunction core collisionmodels apputils)