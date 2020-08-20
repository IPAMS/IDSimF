project(IDSimF_core)

set(SOURCE_FILES
        Core_constants.hpp
        Core_vector.cpp
        Core_vector.hpp
        Core_randomGenerators.cpp
        Core_randomGenerators.hpp
        Core_randomTestSamples.hpp
        Core_debug.hpp)

add_library(core STATIC ${SOURCE_FILES})
target_include_directories(core PUBLIC .)