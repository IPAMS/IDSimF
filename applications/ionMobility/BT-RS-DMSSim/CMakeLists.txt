project(BT-RS-DMSSim)

set(SOURCE_FILES
        BT-RS-DMSSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils rs collisionmodels particlesimulation)

add_test(NAME BT-RS-DMSSim_waterCluster_benchmarkRun COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_benchmarkRun.json" "run_BT-RS-DMSSim_waterCluster_benchmarkRun")

add_test(NAME BT-RS-DMSSim_waterCluster_benchmarkRun_noCollison COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_benchmarkRun_noCollision.json" "run_BT-RS-DMSSim_waterCluster_benchmarkRun_noCollision")

add_test(NAME BT-RS-DMSSim_waterCluster_temperatureGradient COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_temperatureGradient.json" "run_BT-RS-DMSSim_waterCluster_temperatureGradient")

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)