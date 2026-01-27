project(benchmark_matrix3Performance)

set(SOURCE_FILES
        matrix3_performanceBenchmark.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

#find_package (Eigen3 3.3 REQUIRED NO_MODULE)

target_link_libraries(benchmark_matrix3Performance core apputils)