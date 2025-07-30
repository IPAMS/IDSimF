project(collisionSamplingSim)

set(SOURCE_FILES
        collisionSamplingSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils particlesimulation file_io integration spacecharge collisionmodels rs)

# register example runs for CTest
add_test(NAME app_ionCollision_collisionSamplingSim_Li+_N2 COMMAND ${PROJECT_NAME}
        "example/collisionSampling_Li_N2.json" "run_app_ionCollision_collisionSamplingSim_Li_N2" -n ${N_THREADS})

add_test(NAME app_ionCollision_collisionSamplingSim_Li+_Ar COMMAND ${PROJECT_NAME}
        "example/collisionSampling_O2_Ar.json" "run_app_ionCollision_collisionSamplingSim_O2_Ar" -n ${N_THREADS})

add_test(NAME app_ionCollision_collisionSamplingSim_He_He COMMAND ${PROJECT_NAME}
        "example/collisionSampling_He_he.json" "run_app_ionCollision_collisionSamplingSim_He_He" -n ${N_THREADS})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)