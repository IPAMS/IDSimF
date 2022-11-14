project(idealIsothermReactorSim)

set(SOURCE_FILES
        idealIsothermReactorSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils spacecharge rs particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})

# register example runs for CTest
add_test(NAME app_chemistry_idealIsothermReactorSim_waterCluster_benchmarkRun COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_benchmarkRun.json" "run_app_chemisty_idealIsothermReactorSim_waterCluster_benchmarkRun"
        -n ${N_THREADS})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)