project(quadrupoleCollisionCellSim)

set(SOURCE_FILES
        quadrupoleCollisionCellSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})

# register example runs for CTest
add_test(NAME app_ionCollision_quadrupoleCollisionCellSim_simpleCapacitor COMMAND ${PROJECT_NAME}
        "example/BTree_quadCollisionCell_simpleCapacitor.json" "run_app_ionCollision_quadrupoleCollisionCellSim_simpleCapacitor" -n ${N_THREADS})

add_test(NAME app_ionCollision_quadrupoleCollisionCellSim_simpleCapacitor_simpleRec COMMAND ${PROJECT_NAME}
        "example/BTree_quadCollisionCell_simpleCapacitor_simpleRec.json" "run_app_ionCollision_quadrupoleCollisionCellSim_simpleCapacitor_simpleRec"
        -n ${N_THREADS})

add_test(NAME app_ionCollision_quadrupoleCollisionCellSim_simpleCapacitor_cylinder COMMAND ${PROJECT_NAME}
        "example/BTree_quadCollisionCell_simpleCapacitor_cylinder.json" "run_app_ionCollision_quadrupoleCollisionCellSim_simpleCapacitor_cylinder"
        -n ${N_THREADS})

add_test(NAME app_ionCollision_quadrupoleCollisionCellSim_simpleCapacitor_ionCloud COMMAND ${PROJECT_NAME}
        "example/BTree_quadCollisionCell_simpleCapacitor_ionCloud.json" "run_app_ionCollision_quadrupoleCollisionCellSim_simpleCapacitor_ionCloud"
        -n ${N_THREADS})

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)