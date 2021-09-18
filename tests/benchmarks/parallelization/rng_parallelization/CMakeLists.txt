project(benchmark_rng)

set(SOURCE_FILES
        rng_runtime_benchmark.cpp)

add_executable(benchmark_rng ${SOURCE_FILES})
target_include_directories(benchmark_rng PRIVATE ${CMAKE_SOURCE_DIR}/libs/CLI11)
target_link_libraries(benchmark_rng core apputils)