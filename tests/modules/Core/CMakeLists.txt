project(test_core)

set(SOURCE_FILES
        test_main.cpp
        test_vector.cpp
        test_particle.cpp
        test_utils.cpp
        test_math.cpp)

add_executable(test_core ${SOURCE_FILES})
target_include_directories(test_core PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_core core)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_core OpenMP::OpenMP_CXX)
endif()
add_test(NAME test_core COMMAND test_core)