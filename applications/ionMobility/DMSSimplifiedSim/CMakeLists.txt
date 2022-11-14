project(DMSSimplifiedSim)

set(SOURCE_FILES
        DMSSimplifiedSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils rs collisionmodels particlesimulation apputils_ionMobility)

add_test(NAME app_ionMobility_DMSSimplifiedSim_waterCluster_staticRun COMMAND ${PROJECT_NAME}
        "example/DMS_waterCluster_static.json" "run_app_ionMobility_waterCluster_static" -n ${N_THREADS})

add_test(NAME app_ionMobility_DMSSimplifiedSim_waterCluster_tempDependentRun COMMAND ${PROJECT_NAME}
        "example/DMS_waterCluster_tempDependent.json" "run_app_ionMobility_waterCluster_tempDependent" -n ${N_THREADS})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)