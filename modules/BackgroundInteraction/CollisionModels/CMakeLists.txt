project(IDSimF_collisionModels)

set(SOURCE_FILES
        CollisionModel_MathFunctions.cpp
        CollisionModel_SpatialFieldFunctions.cpp
        CollisionModel_util.cpp
        CollisionModel_MultiCollisionModel.cpp
        HardSoftSphere/CollisionModel_HardSphere.cpp
        HardSoftSphere/CollisionModel_SoftSphere.cpp
        StatisticalDiffusion/CollisionModel_CollisionStatistics.cpp
        StatisticalDiffusion/CollisionModel_StatisticalDiffusion.cpp
        MDInteractions/CollisionModel_MDInteractions.cpp
        MDInteractions/CollisionModel_MDInteractionsTrajectorySampler.cpp
        MDInteractions/CollisionModel_AbstractMDForceField.cpp
        MDInteractions/CollisionModel_Atom.cpp
        MDInteractions/CollisionModel_Molecule.cpp
        MDInteractions/CollisionModel_MolecularStructure.cpp
        MDInteractions/CollisionModel_MDForceField_LJ12_6.cpp
        MDInteractions/CollisionModel_MDForceField_Buckingham.cpp
        MDInteractions/CollisionModel_MDTrajectoryWriter.cpp
        MDInteractions/CollisionModel_HDF5MDTrajectoryWriter.cpp)

add_library(collisionmodels STATIC ${SOURCE_FILES})
target_include_directories(collisionmodels PUBLIC . HardSoftSphere StatisticalDiffusion MDInteractions)
target_link_libraries(collisionmodels PRIVATE core apputils particlesimulation file_io spacecharge rs spdlog::spdlog)