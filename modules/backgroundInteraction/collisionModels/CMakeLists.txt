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
        CollisionModel_EmptyCollisionModel.cpp
        CollisionModel_EmptyCollisionModel.hpp
        CollisionModel_util.cpp
        CollisionModel_util.hpp
        CollisionModel_MultiCollisionModel.cpp
        CollisionModel_MultiCollisionModel.hpp
        CollisionModel_CollisionStatistics.cpp
        CollisionModel_CollisionStatistics.hpp
        CollisionStatistic_default.hpp)

add_library(collisionmodels STATIC ${SOURCE_FILES})
target_include_directories(collisionmodels PUBLIC .)
target_link_libraries(collisionmodels core btree rs)