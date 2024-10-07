project(idealizedMultipleCollisions)

set(SOURCE_FILES
        idealizedMultipleCollisions.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils rs collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)