project(IDSimF_RS)

set(SOURCE_FILES
        RS_AbstractReaction.hpp
        RS_AbstractReaction.cpp
        RS_Simulation.hpp
        RS_Simulation.cpp
        RS_Substance.hpp
        RS_Substance.cpp
        RS_ConfigFileParser.hpp
        RS_ConfigFileParser.cpp
        RS_StaticReaction.hpp
        RS_StaticReaction.cpp
        RS_VantHoffReaction.hpp
        RS_VantHoffReaction.cpp
        RS_FieldDependentVantHoffReaction.hpp
        RS_FieldDependentVantHoffReaction.cpp
        RS_constants.hpp
        RS_ReactiveParticle.cpp
        RS_ReactiveParticle.hpp
        RS_SimulationConfiguration.cpp
        RS_SimulationConfiguration.hpp
        RS_ConcentrationFileWriter.cpp
        RS_ConcentrationFileWriter.hpp
        RS_SimpleCollisionStepReaction.cpp
        RS_SimpleCollisionStepReaction.hpp
        RS_StaticThermalizingReaction.cpp
        RS_StaticThermalizingReaction.hpp
        RS_util.cpp
        RS_util.hpp)

add_library(rs STATIC ${SOURCE_FILES})
target_include_directories(rs PUBLIC .)
target_link_libraries(rs core spacecharge spdlog::spdlog)
if(OpenMP_CXX_FOUND)
    target_link_libraries(rs OpenMP::OpenMP_CXX)
endif()