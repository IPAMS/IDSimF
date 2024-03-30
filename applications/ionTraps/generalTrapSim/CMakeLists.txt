project(generalTrapSim)

set(SOURCE_FILES
        generalTrapSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(${PROJECT_NAME} spdlog::spdlog)
#target_include_directories(${PROJECT_NAME} PRIVATE
#        ${CMAKE_SOURCE_DIR}/libs)

add_test(NAME app_ionTraps_generalTrapSim_simpleCapacitor COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor.json" "run_app_ionTraps_generalTrapSim_simpleCapacitor" -n ${N_THREADS})

add_test(NAME app_ionTraps_generalTrapSim_simpleCapacitor_allSplat COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor_allSplat.json" "run_app_ionTraps_generalTrapSim_simpleCapacitor_allSplat" -n ${N_THREADS})

add_test(NAME app_ionTraps_generalTrapSim_simpleCapacitor_RK4 COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor_RK4.json" "run_app_ionTraps_generalTrapSim_simpleCapacitor_RK4" -n ${N_THREADS})

add_test(NAME app_ionTraps_generalTrapSim_simpleCapacitor_fullSumRK4 COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor_fullSumRK4.json" "run_app_ionTraps_generalTrapSim_simpleCapacitor_fullSumRK4" -n ${N_THREADS})

add_test(NAME app_ionTraps_generalTrapSim_simpleCapacitor_ionCloud COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor_ionCloud.json" "run_app_ionTraps_generalTrapSim_simpleCapacitor_ionCloud" -n ${N_THREADS})

add_test(NAME app_ionTraps_generalTrapSim_simpleCapacitor_serial COMMAND ${PROJECT_NAME}
        "example/simpleCapacitor_serial.json" "run_app_ionTraps_generalTrapSim_simpleCapacitor_serial" -n ${N_THREADS})

if(USE_FMM_3D)
    add_test(NAME app_ionTraps_generalTrapSim_simpleCapacitor_FMM3D COMMAND ${PROJECT_NAME}
            "example/simpleCapacitor_serial.json" "run_app_ionTraps_generalTrapSim_simpleCapacitor_FMM3D" -n ${N_THREADS})
endif()

if(USE_EXAFMMT)
    add_test(NAME app_ionTraps_generalTrapSim_simpleCapacitor_exaFMM COMMAND ${PROJECT_NAME}
            "example/simpleCapacitor_serial.json" "run_app_ionTraps_generalTrapSim_simpleCapacitor_exaFMM" -n ${N_THREADS})
endif()


add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)