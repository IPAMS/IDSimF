project(IDSimF_RS)

set(SOURCE_FILES
        RS_AbstractReaction.cpp
        RS_Simulation.cpp
        RS_Substance.cpp
        RS_ConfigFileParser.cpp
        RS_StaticReaction.cpp
        RS_VantHoffReaction.cpp
        RS_FieldDependentVantHoffReaction.cpp
        RS_ReactiveParticle.cpp
        RS_SimulationConfiguration.cpp
        RS_ConcentrationFileWriter.cpp
        RS_SimpleCollisionStepReaction.cpp
        RS_StaticThermalizingReaction.cpp
        RS_util.cpp)

add_library(rs STATIC ${SOURCE_FILES})
target_include_directories(rs PUBLIC .)
target_link_libraries(rs PRIVATE core spacecharge spdlog::spdlog)
if(OpenMP_CXX_FOUND)
    target_link_libraries(rs PRIVATE OpenMP::OpenMP_CXX)
endif()