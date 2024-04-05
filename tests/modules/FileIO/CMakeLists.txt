project(test_FileIO)

set(SOURCE_FILES
        ${CMAKE_SOURCE_DIR}/libs/jsoncpp/jsoncpp.cpp
        test_main.cpp
        test_vtkFileWriter.cpp
        test_ionCloudReader.cpp
        test_simpleFileWriters.cpp
        test_inductionCurrentWriter.cpp
        test_hdf5FileWriter.cpp
        test_hdf5FileReader.cpp
        test_MolecularStructureReader.cpp)

set(TEST_FILE_FOLDER ${CMAKE_SOURCE_DIR}/tests/testfields)
set(TEST_FILES
        ${TEST_FILE_FOLDER}/swift_test_sin.csv
        ${TEST_FILE_FOLDER}/swift_test_transient.csv
        ${TEST_FILE_FOLDER}/low_sample_waveform.csv
        ${TEST_FILE_FOLDER}/test_linear_vector_field_01.h5
        ${TEST_FILE_FOLDER}/test_linear_scalar_field_01.h5
        ${TEST_FILE_FOLDER}/test_ion_cloud_01.csv
        ${TEST_FILE_FOLDER}/test_particle_trajectory.h5
        ${TEST_FILE_FOLDER}/quad_dev_flow_3d.h5
        ${TEST_FILE_FOLDER}/quad_dev_pressure_3d.h5
        ${TEST_FILE_FOLDER}/quad_dev_pressure_3d.h5
        ${TEST_FILE_FOLDER}/capacitor_test_pa/cylinder_capacitor.pa1
        ${TEST_FILE_FOLDER}/capacitor_test_pa/cylinder_capacitor.pa2
        ${TEST_FILE_FOLDER}/test_molecularstructure_reader.json )
file(COPY ${TEST_FILES} DESTINATION .)

add_executable(test_FileIO ${SOURCE_FILES})
target_include_directories(test_FileIO PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/libs/jsoncpp/json
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_FileIO core file_io spacecharge collisionmodels rs ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

add_test(NAME test_FileIO COMMAND test_FileIO)