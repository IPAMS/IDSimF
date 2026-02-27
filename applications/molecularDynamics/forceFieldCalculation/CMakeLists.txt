project(forceFieldCalculation)

set(SOURCE_FILES
        forceFieldCalculation.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils particlesimulation integration file_io spacecharge collisionmodels rs)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

add_test(NAME app_MDinteractions_forceFieldCalculation COMMAND ${PROJECT_NAME} "example/MD_Collisions_ForceField.json" "run_app_mdForceFields")

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)