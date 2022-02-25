project(BT-RS-IMSSim)

set(SOURCE_FILES
        BT-RS-IMSSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils spacecharge rs collisionmodels particlesimulation)

add_test(NAME BT-RS-IMSSim_waterCluster_HS COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_HS.json" "run_BT-RS-IMSSim_waterCluster_HS")

add_test(NAME BT-RS-IMSSim_waterCluster_SDS COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_SDS.json" "run_BT-RS-IMSSim_waterCluster_SDS")

add_test(NAME BT-RS-IMSSim_no_reaction COMMAND ${PROJECT_NAME}
        "example/RS_no_reaction.json" "run_BT-RS-IMSSim_no_reaction")

add_test(NAME BT-RS-IMSSim_waterCluster_no_transport COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_no_transport.json" "run_BT-RS-IMSSim_waterCluster_no_transport")

add_test(NAME BT-RS-IMSSim_waterCluster_simple COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_simple.json" "run_BT-RS-IMSSim_waterCluster_simple")

add_test(NAME BT-RS-IMSSim_waterCluster_thermalizing COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_thermalizing.json" "run_BT-RS-IMSSim_waterCluster_thermalizing")

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)