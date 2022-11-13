project(generalQuadSim)

set(SOURCE_FILES
        generalQuadSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})

add_test(NAME app_ionTransfer_generalQuadSim_benchmarkRun COMMAND ${PROJECT_NAME}
        "example/BTree_quad_benchmarkRun.json" "run_app_ionTransfer_generalQuadSim_benchmarkRun")

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)