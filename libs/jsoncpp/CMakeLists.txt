project(jsoncpp)

set(SOURCE_FILES
    jsoncpp.cpp)

add_library(jsoncpp STATIC ${SOURCE_FILES})
target_include_directories(jsoncpp PUBLIC json)