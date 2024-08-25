project(test_Applications)

set(SOURCE_FILES
        test_main.cpp
        test_ionDefinitionReading.cpp
        test_logging_timing.cpp
        test_simulation_configuration.cpp)

set(TEST_FILE_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields/simulation_configurations)
set(TEST_PA_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields/capacitor_test_pa)

set(TEST_FILES
        ${TEST_FILE_FOLDER}/ionBox.json
        ${TEST_FILE_FOLDER}/ionCylinder_minimal.json
        ${TEST_FILE_FOLDER}/ionCylinder_full.json
        ${TEST_FILE_FOLDER}/ionCloudFile.json
        ${TEST_FILE_FOLDER}/test_ion_cloud_01.csv
        ${TEST_FILE_FOLDER}/ionDefinition_invalid.json
        ${TEST_FILE_FOLDER}/ionDefinition_invalid_2.json
        ${TEST_FILE_FOLDER}/simulationConfiguration_typesTest.json
        ${TEST_FILE_FOLDER}/simulationConfiguration_fileReadingTest.json
        ${TEST_PA_FOLDER}/cylinder_capacitor.pa1
        ${TEST_PA_FOLDER}/cylinder_capacitor.pa2
        ${TEST_PA_FOLDER}/cylinder_capacitor.pa3
        )
file(COPY ${TEST_FILES} DESTINATION .)

add_executable(test_Applications ${SOURCE_FILES})
target_include_directories(test_Applications PRIVATE
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_Applications core spacecharge particlesimulation apputils)

add_test(NAME test_Applications COMMAND test_Applications)