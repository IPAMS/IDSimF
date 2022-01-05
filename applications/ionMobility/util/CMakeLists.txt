project(apputils_ionMobility)

set(SOURCE_FILES
        dmsSim_dmsFields.hpp
        dmsSim_dmsFields.cpp)

add_library(apputils_ionMobility STATIC ${SOURCE_FILES})
target_include_directories(apputils_ionMobility PUBLIC
        ${CMAKE_SOURCE_DIR}/applications/ionMobility/util)

target_link_libraries(apputils_ionMobility apputils)
