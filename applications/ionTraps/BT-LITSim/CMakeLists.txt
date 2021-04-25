project(BT-LITSim)

set(SOURCE_FILES
        BT-LITSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

add_test(NAME BT-LITSim_simpleCapacitor COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor.json" "run_BT-LITSim_simpleCapacitor")

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)