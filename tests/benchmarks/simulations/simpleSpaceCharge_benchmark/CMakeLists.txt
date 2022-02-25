project(benchmark_simpleSpaceCharge)

set(SOURCE_FILES
        simpleSpaceCharge_benchmark.cpp
)

add_executable(benchmark_simpleSpaceCharge ${SOURCE_FILES})
target_include_directories(benchmark_simpleSpaceCharge PRIVATE
        ${CMAKE_SOURCE_DIR}/libs/CLI11)
target_link_libraries(benchmark_simpleSpaceCharge core spacecharge particlesimulation collisionmodels apputils)