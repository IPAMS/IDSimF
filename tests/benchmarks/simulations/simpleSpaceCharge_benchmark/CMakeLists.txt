project(benchmark_simpleSpaceCharge)

set(SOURCE_FILES
        simpleSpaceCharge_benchmark.cpp
)

add_executable(benchmark_simpleSpaceCharge ${SOURCE_FILES})

target_link_libraries(benchmark_simpleSpaceCharge core btree particlesimulation collisionmodels apputils)