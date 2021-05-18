project(apputils)

set(SOURCE_FILES
        ${CMAKE_SOURCE_DIR}/libs/jsoncpp/jsoncpp.cpp
        appUtils_simulationConfiguration.hpp
        appUtils_inputFileUtilities.hpp
        appUtils_inputFileUtilities.cpp
        appUtils_ionDefinitionReading.cpp
        appUtils_ionDefinitionReading.hpp appUtils_simulationConfiguration.cpp appUtils_logging.cpp appUtils_logging.hpp appUtils_stopwatch.cpp appUtils_stopwatch.hpp appUtils_signalHandler.cpp appUtils_signalHandler.hpp)

add_library(apputils STATIC ${SOURCE_FILES})
target_include_directories(apputils PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/jsoncpp/json
        ${CMAKE_SOURCE_DIR}/applications/util)

target_link_libraries(apputils particlesimulation spdlog::spdlog)
if (USE_CPP_FSLIB)
    target_link_libraries(apputils stdc++fs)
endif()