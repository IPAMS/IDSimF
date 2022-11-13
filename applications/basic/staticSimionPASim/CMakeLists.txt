project(staticSimionPASim)

set(SOURCE_FILES
        staticSimionPASim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils spacecharge collisionmodels particlesimulation)

# register example runs for CTest
add_test(NAME app_basic_staticSimionPASim_benchmarkRunShort COMMAND ${PROJECT_NAME} "example/benchmarkRun_short.json" "run_app_basic_staticSimionPASim_benchmarkRunShort")

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)