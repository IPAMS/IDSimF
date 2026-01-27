project(benchmark_vectorCompareImplementations)

set(SOURCE_FILES
        vectorImplementationsCompare.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

target_link_libraries(benchmark_vectorCompareImplementations core apputils Eigen3::Eigen)