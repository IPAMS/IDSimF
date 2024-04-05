project(spaceChargeSimpleSim)

set(SOURCE_FILES
        spaceChargeSimpleSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils spacecharge collisionmodels particlesimulation)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

# register example runs for CTest
add_test(NAME app_basic_spaceChargeSimpleSim_benchmarkRunShort_parallel COMMAND ${PROJECT_NAME} "example/benchmarkRun_short_parallel.json" "run_app_basic_spaceChargeSimpleSim_benchmarkRunShort_parallel")
add_test(NAME app_basic_spaceChargeSimpleSim_benchmarkRunShort_serial COMMAND ${PROJECT_NAME} "example/benchmarkRun_short_serial.json" "run_app_basic_spaceChargeSimpleSim_benchmarkRunShort_serial")

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/example/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/example/)