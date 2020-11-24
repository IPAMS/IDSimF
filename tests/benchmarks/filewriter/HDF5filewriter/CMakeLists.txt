project(benchmark_hdf5Writer)

set(SOURCE_FILES
        HDF5filewriter_benchmark.cpp
)

add_executable(benchmark_hdf5Writer ${SOURCE_FILES})

target_link_libraries(benchmark_hdf5Writer core particlesimulation)