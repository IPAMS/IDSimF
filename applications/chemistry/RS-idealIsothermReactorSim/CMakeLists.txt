project(RS-idealIsothermReactorSim)

set(SOURCE_FILES
        RS-idealIsothermReactorSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils btree rs particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})

# register example runs for CTest
add_test(NAME RS-idealIsothermReactorSim_waterCluster_benchmarkRun COMMAND ${PROJECT_NAME} "example/RS_waterCluster_benchmarkRun.json" "run_RS-idealIsothermReactorSim_waterCluster_benchmarkRun")

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)