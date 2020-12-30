cmake_minimum_required(VERSION 3.5)
project(ionSimulation_cpp)

set(SOURCE_FILES
        rs_runtime_benchmark.cpp)

set(TEST_FILES
        RS_waterCl_static.conf
)
file(COPY ${TEST_FILES} DESTINATION .)

add_executable(benchmark_rs_runtime ${SOURCE_FILES})
target_link_libraries(benchmark_rs_runtime core particlesimulation rs)