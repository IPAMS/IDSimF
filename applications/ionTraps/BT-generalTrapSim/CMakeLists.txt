project(BT-generalTrapSim)

set(SOURCE_FILES
        BT-generalTrapSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

add_test(NAME BT-generalTrapSim_simpleCapacitor COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor.json" "run_BT-generalTrapSim_simpleCapacitor")

add_test(NAME BT-generalTrapSim_simpleCapacitor_allSplat COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor_allSplat.json" "run_BT-generalTrapSim_simpleCapacitor_allSplat")

add_test(NAME BT-generalTrapSim_simpleCapacitor_ionCloud COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor_ionCloud.json" "run_BT-generalTrapSim_simpleCapacitor_ionCloud")

add_test(NAME BT-generalTrapSim_simpleCapacitor_serial COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor_serial.json" "run_BT-generalTrapSim_simpleCapacitor_serial")

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)