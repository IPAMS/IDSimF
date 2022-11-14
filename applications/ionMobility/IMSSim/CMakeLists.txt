project(IMSSim)

set(SOURCE_FILES
        IMSSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils spacecharge rs collisionmodels particlesimulation)

add_test(NAME app_ionMobility_IMSSim_waterCluster_HS COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_HS.json" "run_BT-RS-IMSSim_waterCluster_HS" -n ${N_THREADS})

add_test(NAME app_ionMobility_IMSSim_waterCluster_SDS COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_SDS.json" "run_app_ionMobility_IMSSim_waterCluster_SDS" -n ${N_THREADS})

add_test(NAME app_ionMobility_IMSSim_no_reaction COMMAND ${PROJECT_NAME}
        "example/RS_no_reaction.json" "run_app_ionMobility_IMSSim_no_reaction" -n ${N_THREADS})

add_test(NAME app_ionMobility_IMSSim_waterCluster_no_transport COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_no_transport.json" "run_app_ionMobility_IMSSim_waterCluster_no_transport" -n ${N_THREADS})

add_test(NAME app_ionMobility_IMSSim_waterCluster_simple COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_simple.json" "run_app_ionMobility_IMSSim_waterCluster_simple" -n ${N_THREADS})

add_test(NAME app_ionMobility_IMSSim_waterCluster_thermalizing COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_thermalizing.json" "run_app_ionMobility_IMSSim_waterCluster_thermalizing" -n ${N_THREADS})

add_test(NAME app_ionMobility_IMSSim_MD_short_no_recording COMMAND ${PROJECT_NAME}
        "example/RS_MDCollisions_O2-He_short_no_recording.json" "run_app_ionMobility_IMSSim_RS_MDCollisions_O2-He_short_no_recording" -n ${N_THREADS})

add_test(NAME app_ionMobility_IMSSim_MD_short_recording COMMAND ${PROJECT_NAME}
        "example/RS_MDCollisions_O2-He_short_recording.json" "run_app_ionMobility_IMSSim_RS_MDCollisions_O2-He_short_recording" -n ${N_THREADS})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)