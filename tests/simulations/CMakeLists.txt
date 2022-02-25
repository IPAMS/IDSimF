project(test_simulations)

set(SOURCE_FILES
        test_main.cpp
        test_integrators.cpp
)

add_executable(test_simulations ${SOURCE_FILES})
target_include_directories(test_simulations PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_simulations core spacecharge integration particlesimulation collisionmodels file_io)

add_test(NAME test_simulations COMMAND test_simulations)