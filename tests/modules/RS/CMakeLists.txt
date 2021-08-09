set(SOURCE_FILES
        test_main.cpp
        test_reaction.cpp
        test_simulation.cpp
        test_substance.cpp
        test_configFileParser.cpp
        test_reactiveParticle.cpp
        test_simulationConfiguration.cpp test_util.cpp)

set(TEST_FILE_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields)
set(TEST_FILES
        ${TEST_FILE_FOLDER}/RS_waterCluster_test.conf
        ${TEST_FILE_FOLDER}/RS_waterCluster_test_temperatureDependent.conf
        ${TEST_FILE_FOLDER}/RS_waterCluster_test_fieldDependent.conf
        ${TEST_FILE_FOLDER}/RS_collisionBasedReactions_test.conf
        ${TEST_FILE_FOLDER}/RS_minimal_test.conf
        )
file(COPY ${TEST_FILES} DESTINATION .)

add_executable(test_RS ${SOURCE_FILES})
target_include_directories(test_RS PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_RS core btree rs)

add_test(NAME test_RS COMMAND test_RS)