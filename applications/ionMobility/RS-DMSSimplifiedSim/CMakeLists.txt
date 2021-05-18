project(RS-DMSSimplifiedSim)

set(SOURCE_FILES
        RS-DMSSimplifiedSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils rs collisionmodels particlesimulation)

add_test(NAME BT-RS-DMSSimplifiedSim_waterCluster_staticRun COMMAND ${PROJECT_NAME}
        "example/DMS_waterCluster_static.json" "run_DMS_waterCluster_static")

add_test(NAME BT-RS-DMSSimplifiedSim_waterCluster_tempDependentRun COMMAND ${PROJECT_NAME}
        "example/DMS_waterCluster_tempDependent.json" "run_DMS_waterCluster_tempDependent")

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)