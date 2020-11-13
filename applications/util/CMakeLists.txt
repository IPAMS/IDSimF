project(apputils)

set(SOURCE_FILES
        ${CMAKE_SOURCE_DIR}/libs/jsoncpp/jsoncpp.cpp
        appUtils_parameterParsing.hpp
        appUtils_inputFileUtilities.hpp
        appUtils_inputFileUtilities.cpp
        appUtils_ionDefinitionReading.cpp
        appUtils_ionDefinitionReading.hpp appUtils_parameterParsing.cpp)

add_library(apputils STATIC ${SOURCE_FILES})
target_include_directories(apputils PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/jsoncpp/json
        ${CMAKE_SOURCE_DIR}/applications/util)

target_link_libraries(apputils particlesimulation)
if (USE_CPP_FSLIB)
    target_link_libraries(apputils stdc++fs)
endif()