project(TIMSSim)

set(SOURCE_FILES
        TIMSSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils rs collisionmodels particlesimulation apputils_ionMobility)

add_test(NAME app_ionMobility_TIMSSim_HS_uniformFlow COMMAND ${PROJECT_NAME}
        "example/TIMS_HS_uniformFlow.json" "run_IMS_HS_uniformFlow" -n ${N_THREADS})

add_test(NAME app_ionMobility_TIMSSim_HS_parabolicFlow COMMAND ${PROJECT_NAME}
        "example/TIMS_HS_parabolicFlow.json" "run_IMS_HS_parabolicFlow" -n ${N_THREADS})

add_test(NAME app_ionMobility_TIMSSim_HS_flowfield_2d_axial_symmetric COMMAND ${PROJECT_NAME}
        "example/TIMS_HS_flowfield_2dAxi.json" "TIMS_HS_flowfield_2dAxi" -n ${N_THREADS})

add_test(NAME app_ionMobility_TIMSSim_HS_flowfield_3d COMMAND ${PROJECT_NAME}
        "example/TIMS_HS_flowfield_3d.json" "run_TIMS_HS_flowfield_3d" -n ${N_THREADS})


add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)