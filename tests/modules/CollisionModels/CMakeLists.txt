project(test_CollisionModels)

set(SOURCE_FILES
        test_HardSphere.cpp
        test_main.cpp
        test_MathFunctions.cpp
        test_CollisionStatistics.cpp
        test_energyDiatomicMD.cpp
        test_StatisticalDiffusion.cpp
        test_util.cpp
        test_MultiCollisionModel.cpp 
        test_SoftSphere.cpp
        test_MDInteractions.cpp
        test_Atom.cpp
        test_Molecule.cpp
        test_MDVerlet.cpp
        test_preconstructedMD.cpp
        test_SpatialFunctions.cpp)

set(TEST_FILE_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields)
set(TEST_FILES
        ${TEST_FILE_FOLDER}/cs_icdf_2020_02_27_001_test.dat
        ${TEST_FILE_FOLDER}/cs_icdf_2020_03_15_001_parameterTest.dat
        ${TEST_FILE_FOLDER}/test_molecularstructure_reader.json
        )

file(COPY ${TEST_FILES} DESTINATION .)

add_executable(test_CollisionModels ${SOURCE_FILES})
target_include_directories(test_CollisionModels PUBLIC ${CMAKE_SOURCE_DIR}/libs/catch ${CMAKE_SOURCE_DIR}/tests/util ${CMAKE_SOURCE_DIR}/libs/jsoncpp/json
                                ${CMAKE_SOURCE_DIR}/applications/util)
target_link_libraries(test_CollisionModels core spacecharge integration rs collisionmodels particlesimulation)
target_link_libraries(test_CollisionModels ${VTK_LIBRARIES})

add_test(NAME test_CollisionModels COMMAND test_CollisionModels)

add_executable(test_MDPre test_preconstructedMD.cpp test_main.cpp)
target_include_directories(test_MDPre PUBLIC ${CMAKE_SOURCE_DIR}/libs/catch ${CMAKE_SOURCE_DIR}/tests/util ${CMAKE_SOURCE_DIR}/libs/jsoncpp/json
                                ${CMAKE_SOURCE_DIR}/applications/util)
target_link_libraries(test_MDPre core spacecharge integration rs collisionmodels particlesimulation)
target_link_libraries(test_MDPre ${VTK_LIBRARIES})
add_test(NAME test_MDPre COMMAND test_MDPre)
set_tests_properties(test_MDPre PROPERTIES TIMEOUT 7200) 