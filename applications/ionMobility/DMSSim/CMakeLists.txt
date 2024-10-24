project(DMSSim)

set(SOURCE_FILES
        DMSSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils rs collisionmodels particlesimulation apputils_ionMobility)

add_test(NAME app_ionMobility_DMSSim_waterCluster_HS COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_HS.json" "run_BT-RS-DMSSim_waterCluster_HS" -n ${N_THREADS})

add_test(NAME app_ionMobility_DMSSim_waterCluster_SDS COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_SDS.json" "run_app_ionMobility_DMSSim_waterCluster_SDS" -n ${N_THREADS})

add_test(NAME app_ionMobility_DMSSim_waterCluster_noCollison COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_noCollision.json" "run_app_ionMobility_DMSSim_waterCluster_noCollision" -n ${N_THREADS})

add_test(NAME app_ionMobility_DMSSim_waterCluster_temperatureGradient COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_temperatureGradient.json" "run_app_ionMobility_DMSSim_waterCluster_temperatureGradient" -n ${N_THREADS})

add_test(NAME app_ionMobility_DMSSim_waterCluster_modulatedCV COMMAND ${PROJECT_NAME}
        "example/RS_waterCluster_modulatedCV.json" "run_app_ionMobility_DMSSim_waterCluster_modulatedCV" -n ${N_THREADS})

add_test(NAME app_ionMobility_DMSSim_MD_colllisions COMMAND ${PROJECT_NAME}
        "example/RS_MDCollisions_O2-He_short.json" "run_app_ionMobility_MDCollisions_O2-He_short" -n ${N_THREADS})


add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)