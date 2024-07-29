project(apputils)

set(SOURCE_FILES
        ${CMAKE_SOURCE_DIR}/libs/jsoncpp/jsoncpp.cpp
        appUtils_simulationConfiguration.hpp
        appUtils_inputFileUtilities.hpp
        appUtils_inputFileUtilities.cpp
        appUtils_ionDefinitionReading.cpp
        appUtils_ionDefinitionReading.hpp
        appUtils_simulationConfiguration.cpp
        appUtils_logging.cpp
        appUtils_logging.hpp
        appUtils_stopwatch.cpp
        appUtils_stopwatch.hpp
        appUtils_signalHandler.cpp
        appUtils_signalHandler.hpp
        appUtils_commandlineParser.cpp
        appUtils_commandlineParser.hpp
        appUtils_integrationRunning.cpp
        appUtils_integrationRunning.hpp)

add_library(apputils STATIC ${SOURCE_FILES})
target_include_directories(apputils PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/jsoncpp/json
        ${CMAKE_SOURCE_DIR}/applications/util
        ${CMAKE_SOURCE_DIR}/libs/CLI11)

target_link_libraries(apputils particlesimulation integration spdlog::spdlog)
if(OpenMP_CXX_FOUND)
    target_link_libraries(apputils OpenMP::OpenMP_CXX)
endif()
if (USE_CPP_FSLIB)
    target_link_libraries(apputils stdc++fs)
endif()