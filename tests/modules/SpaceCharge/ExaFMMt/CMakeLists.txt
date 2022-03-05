project(test_ExaFMMt)

set(SOURCE_FILES
        test_main.cpp
        test_exafmmtSolver.cpp)

add_executable(test_ExaFMMt ${SOURCE_FILES})
target_include_directories(test_ExaFMMt PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_ExaFMMt core spacecharge particlesimulation)

add_test(NAME test_ExaFMMt COMMAND test_ExaFMMt)