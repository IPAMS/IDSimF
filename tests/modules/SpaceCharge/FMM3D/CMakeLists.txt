project(test_FMM3D)

set(SOURCE_FILES
        test_main.cpp
        test_fmm3DSolver.cpp)

add_executable(test_FMM3D ${SOURCE_FILES})
target_include_directories(test_FMM3D PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_FMM3D core spacecharge particlesimulation)

add_test(NAME test_FMM3D COMMAND test_FMM3D)