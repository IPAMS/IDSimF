project(LITSim)

set(SOURCE_FILES
        LITSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils particlesimulation file_io integration spacecharge collisionmodels rs)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

add_test(NAME app_ionTraps_LITSim_simpleQuadLIT COMMAND ${PROJECT_NAME}
        "example/simpleQuadLIT.json" "run_app_ionTraps_LITSim_simpleQuadLIT" -n ${N_THREADS})

add_test(NAME app_ionTraps_LITSim_simpleQuadLIT_RungeKutta COMMAND ${PROJECT_NAME}
        "example/simpleQuadLIT_RK4.json" "run_app_ionTraps_LITSim_simpleQuadLIT_RK4" -n ${N_THREADS})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)