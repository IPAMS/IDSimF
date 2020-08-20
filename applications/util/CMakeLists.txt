project(apputils)

set(SOURCE_FILES
        ${CMAKE_SOURCE_DIR}/libs/jsoncpp/jsoncpp.cpp
        parameterParsing.hpp
        inputFileUtilities.hpp
        inputFileUtilities.cpp)

add_library(apputils STATIC ${SOURCE_FILES})
target_include_directories(apputils PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/jsoncpp/json
        ${CMAKE_SOURCE_DIR}/applications/util)

target_link_libraries(apputils particlesimulation)
if (USE_CPP_FSLIB)
    target_link_libraries(apputils stdc++fs)
endif()