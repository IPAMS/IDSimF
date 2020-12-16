project(BT-spaceChargeMinimalSim)

set(SOURCE_FILES
        BT-spaceChargeMinimalSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils btree collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

# register example run for CTest
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} "example/benchmarkRun_short.json" "benchmarkRun_short")

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)