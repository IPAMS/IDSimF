project(TIMSSim)

set(SOURCE_FILES
        TIMSSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils rs collisionmodels particlesimulation apputils_ionMobility)

add_test(NAME app_ionMobility_TWIMSSim_waterCluster_HS COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_HS.json" "run_BT-RS-TWIMSSim_waterCluster_HS" -n ${N_THREADS})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)