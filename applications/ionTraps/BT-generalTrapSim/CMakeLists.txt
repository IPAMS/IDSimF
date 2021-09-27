project(BT-generalTrapSim)

set(SOURCE_FILES
        BT-generalTrapSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(${PROJECT_NAME} spdlog::spdlog)
#target_include_directories(${PROJECT_NAME} PRIVATE
#        ${CMAKE_SOURCE_DIR}/libs)

add_test(NAME BT-generalTrapSim_simpleCapacitor COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor.json" "run_BT-generalTrapSim_simpleCapacitor" -n ${N_THREADS})

add_test(NAME BT-generalTrapSim_simpleCapacitor_allSplat COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor_allSplat.json" "run_BT-generalTrapSim_simpleCapacitor_allSplat" -n ${N_THREADS})

add_test(NAME BT-generalTrapSim_simpleCapacitor_ionCloud COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor_ionCloud.json" "run_BT-generalTrapSim_simpleCapacitor_ionCloud" -n ${N_THREADS})

add_test(NAME BT-generalTrapSim_simpleCapacitor_serial COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor_serial.json" "run_BT-generalTrapSim_simpleCapacitor_serial" -n ${N_THREADS})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)