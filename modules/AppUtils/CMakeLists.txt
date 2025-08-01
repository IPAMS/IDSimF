project(apputils)

set(SOURCE_FILES
        appUtils_inputFileUtilities.cpp
        appUtils_ionDefinitionReading.cpp
        appUtils_simulationConfiguration.cpp
        appUtils_logging.cpp
        appUtils_stopwatch.cpp
        appUtils_signalHandler.cpp
        appUtils_commandlineParser.cpp
        appUtils_integrationRunning.cpp)

add_library(apputils STATIC ${SOURCE_FILES})
target_include_directories(apputils PUBLIC .)
target_include_directories(apputils PUBLIC ${CMAKE_SOURCE_DIR}/libs/CLI11)

target_link_libraries(apputils PRIVATE core particlesimulation collisionmodels integration spacecharge file_io)
target_link_libraries(apputils PUBLIC spdlog::spdlog jsoncpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(apputils PRIVATE OpenMP::OpenMP_CXX)
endif()
if (USE_CPP_FSLIB)
    target_link_libraries(apputils PRIVATE stdc++fs)
endif()