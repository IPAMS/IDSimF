project(BT-QITSim)

set(SOURCE_FILES
        BT-QITSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

add_test(NAME BT-QITSim_fixedRF COMMAND ${PROJECT_NAME} "example/QIT_fixedRF.json" "run_BT-QITSim_fixedRF" -n ${N_THREADS})
add_test(NAME BT-QITSim_fixedRF_ionCloud COMMAND ${PROJECT_NAME} "example/QIT_fixedRF_ionCloud.json" "run_BT-QITSim_fixedRF_ionCloud" -n ${N_THREADS})
add_test(NAME BT-QITSim_particlesScaled COMMAND ${PROJECT_NAME} "example/QIT_particlesScaled.json" "run_BT-QITSim_particlesScaled" -n ${N_THREADS})
add_test(NAME BT-QITSim_variableRF COMMAND ${PROJECT_NAME} "example/QIT_variableRF.json" "run_BT-QITSim_variableRF" -n ${N_THREADS})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)