project(IDSimF_core)

set(SOURCE_FILES
        Core_constants.hpp
        Core_vector.cpp
        Core_vector.hpp
        Core_particle.cpp
        Core_particle.hpp
        Core_randomGenerators.cpp
        Core_randomGenerators.hpp
        Core_randomTestSamples.hpp
        Core_debug.hpp
        Core_math.cpp
        Core_math.hpp)

add_library(core STATIC ${SOURCE_FILES})
if(OpenMP_CXX_FOUND)
    target_link_libraries(core PRIVATE OpenMP::OpenMP_CXX)
endif()
target_include_directories(core PUBLIC .)
# target_include_directories(core PUBLIC . ${CMAKE_SOURCE_DIR}/modules/BackgroundInteraction/CollisionModels)