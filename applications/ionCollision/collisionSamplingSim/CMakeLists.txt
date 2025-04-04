project(collisionSamplingSim)

set(SOURCE_FILES
        collisionSamplingSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils collisionmodels)

# register example runs for CTest
add_test(NAME app_ionCollision_collisionSamplingSim_Li+ COMMAND ${PROJECT_NAME}
        "example/collisionSampling_Li.json.json" "run_app_ionCollision_collisionSamplingSim" -n ${N_THREADS})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)