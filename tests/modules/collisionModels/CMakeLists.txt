project(test_CollisionModels)

set(SOURCE_FILES
        test_HardSphere.cpp
        test_main.cpp
        test_MathFunctions.cpp
        test_CollisionStatistics.cpp
        test_StatisticalDiffusion.cpp
        test_util.cpp
        test_MultiCollisionModel.cpp)

set(TEST_FILE_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields)
set(TEST_FILES
        ${TEST_FILE_FOLDER}/cs_icdf_2020_02_27_001_test.dat
        ${TEST_FILE_FOLDER}/cs_icdf_2020_03_15_001_parameterTest.dat
        )
file(COPY ${TEST_FILES} DESTINATION .)

add_executable(test_CollisionModels ${SOURCE_FILES})
target_include_directories(test_CollisionModels PUBLIC ${CMAKE_SOURCE_DIR}/libs/catch)
target_link_libraries(test_CollisionModels core btree rs collisionmodels particlesimulation)
target_link_libraries(test_CollisionModels ${VTK_LIBRARIES})

add_test(NAME test_CollisionModels COMMAND test_CollisionModels)