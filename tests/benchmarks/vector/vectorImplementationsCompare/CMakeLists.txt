project(benchmark_vectorImplementationCompare)

set(SOURCE_FILES
        vectorImplementationsCompare.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

target_link_libraries(benchmark_vectorImplementationCompare core apputils collisionmodels particlesimulation Eigen3::Eigen)