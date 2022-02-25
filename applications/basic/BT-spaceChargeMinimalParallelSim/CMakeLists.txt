project(BT-spaceChargeMinimalParallelSim)

set(SOURCE_FILES
        BT-spaceChargeMinimalParallelSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils spacecharge collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

# register example runs for CTest
add_test(NAME BT-BT-spaceChargeMinimalParallelSim_benchmarkRunShort COMMAND ${PROJECT_NAME}
        "example/benchmarkRun_short.json" "run_BT-spaceChargeMinimalParallelSim_benchmarkRunShort"
        -n ${N_THREADS})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)