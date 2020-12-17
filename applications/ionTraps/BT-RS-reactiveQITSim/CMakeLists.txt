project(BT-RS-reactiveQITSim)

set(SOURCE_FILES
        BT-RS-reactiveQITSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils rs collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

add_test(NAME BT-RS-reactiveQITSim_benchmarkRun COMMAND ${PROJECT_NAME}
        "example/QIT_benchmarkRun.json" "run_BT-RS-reactiveQITSim_benchmarkRun")

add_test(NAME BT-RS-reactiveQITSim_benchmarkRun_scanningRF COMMAND ${PROJECT_NAME}
        "example/QIT_benchmarkRun_scanningRf.json" "run_BT-RS-reactiveQITSim_benchmarkRun_scanningRF")

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)