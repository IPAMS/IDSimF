project(test_Applications)

set(SOURCE_FILES
        test_main.cpp
        test_ionDefinitionReading.cpp
)

set(TEST_FILE_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields/simulation_configurations)
set(TEST_FILES
        ${TEST_FILE_FOLDER}/ionBox.json
        ${TEST_FILE_FOLDER}/ionCloudFile.json
        ${TEST_FILE_FOLDER}/test_ion_cloud_01.csv
        )
file(COPY ${TEST_FILES} DESTINATION .)

add_executable(test_Applications ${SOURCE_FILES})
target_include_directories(test_Applications PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_Applications core btree particlesimulation apputils)