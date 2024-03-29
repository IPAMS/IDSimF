project(IDSimF_collisionModels)

set(SOURCE_FILES
        CollisionModel_HardSphere.hpp
        CollisionModel_HardSphere.cpp
        CollisionModel_MathFunctions.hpp
        CollisionModel_MathFunctions.cpp
        CollisionModel_SpatialFieldFunctions.hpp
        CollisionModel_SpatialFieldFunctions.cpp
        CollisionModel_StatisticalDiffusion.cpp
        CollisionModel_StatisticalDiffusion.hpp
        CollisionModel_AbstractCollisionModel.hpp
        CollisionModel_util.cpp
        CollisionModel_util.hpp
        CollisionModel_MultiCollisionModel.cpp
        CollisionModel_MultiCollisionModel.hpp
        CollisionModel_CollisionStatistics.cpp
        CollisionModel_CollisionStatistics.hpp
        CollisionStatistic_default.hpp 
        CollisionModel_SoftSphere.cpp 
        CollisionModel_SoftSphere.hpp
        CollisionModel_MDInteractions.cpp
        CollisionModel_MDInteractionsPreconstructed.cpp
        CollisionModel_MDInteractionsPreconstructed.hpp
        CollisionModel_MDInteractions.hpp
        CollisionModel_AbstractMDForceField.cpp
        CollisionModel_AbstractMDForceField.hpp
        CollisionModel_Atom.cpp
        CollisionModel_Atom.hpp
        CollisionModel_Molecule.cpp
        CollisionModel_Molecule.hpp
        CollisionModel_MolecularStructure.cpp
        CollisionModel_MolecularStructure.hpp CollisionModel_MDForceField_LJ12_6.cpp CollisionModel_MDForceField_LJ12_6.hpp)

add_library(collisionmodels STATIC ${SOURCE_FILES})
target_include_directories(collisionmodels PUBLIC .)
target_link_libraries(collisionmodels core spacecharge rs apputils)