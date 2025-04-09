project(IDSimF_collisionModels)

set(SOURCE_FILES
        CollisionModel_MathFunctions.hpp
        CollisionModel_MathFunctions.cpp
        CollisionModel_SpatialFieldFunctions.hpp
        CollisionModel_SpatialFieldFunctions.cpp
        CollisionModel_AbstractCollisionModel.hpp
        CollisionModel_util.cpp
        CollisionModel_util.hpp
        CollisionModel_MultiCollisionModel.cpp
        CollisionModel_MultiCollisionModel.hpp
        HardSoftSphere/CollisionModel_HardSphere.hpp
        HardSoftSphere/CollisionModel_HardSphere.cpp
        HardSoftSphere/CollisionModel_SoftSphere.cpp
        HardSoftSphere/CollisionModel_SoftSphere.hpp
        StatisticalDiffusion/CollisionModel_CollisionStatistics.cpp
        StatisticalDiffusion/CollisionModel_CollisionStatistics.hpp
        StatisticalDiffusion/CollisionModel_StatisticalDiffusion.cpp
        StatisticalDiffusion/CollisionModel_StatisticalDiffusion.hpp
        StatisticalDiffusion/CollisionStatistic_default.hpp
        MDInteractions/CollisionModel_MDInteractions.cpp
        MDInteractions/CollisionModel_MDInteractionsExperimental.cpp
        MDInteractions/CollisionModel_MDInteractionsExperimental.hpp
        MDInteractions/CollisionModel_MDInteractionsTrajectorySampler.cpp
        MDInteractions/CollisionModel_MDInteractionsTrajectorySampler.hpp
        MDInteractions/CollisionModel_MDInteractions.hpp
        MDInteractions/CollisionModel_AbstractMDForceField.cpp
        MDInteractions/CollisionModel_AbstractMDForceField.hpp
        MDInteractions/CollisionModel_Atom.cpp
        MDInteractions/CollisionModel_Atom.hpp
        MDInteractions/CollisionModel_Molecule.cpp
        MDInteractions/CollisionModel_Molecule.hpp
        MDInteractions/CollisionModel_MolecularStructure.cpp
        MDInteractions/CollisionModel_MolecularStructure.hpp
        MDInteractions/CollisionModel_MDForceField_LJ12_6.cpp
        MDInteractions/CollisionModel_MDForceField_LJ12_6.hpp
        MDInteractions/CollisionModel_MDForceField_Buckingham.cpp
        MDInteractions/CollisionModel_MDForceField_Buckingham.hpp
        MDInteractions/CollisionModel_AtomInteractionMap.hpp
        MDInteractions/CollisionModel_MDTrajectoryWriter.cpp
        MDInteractions/CollisionModel_MDTrajectoryWriter.hpp)

add_library(collisionmodels STATIC ${SOURCE_FILES})
target_include_directories(collisionmodels PUBLIC . HardSoftSphere StatisticalDiffusion MDInteractions)
target_link_libraries(collisionmodels core spacecharge rs apputils)