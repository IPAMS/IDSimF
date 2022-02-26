project(IDSimF_spaceCharge)

set(SOURCE_FILES
        BTree/BTree_abstractNode.cpp
        BTree/BTree_abstractNode.hpp
        BTree/BTree_node.cpp
        BTree/BTree_node.hpp
        BTree/BTree_parallelNode.cpp
        BTree/BTree_parallelNode.hpp
        BTree/BTree_treeParticle.cpp
        BTree/BTree_treeParticle.hpp
        BTree/BTree_tree.cpp
        BTree/BTree_tree.hpp
        BTree/BTree_parallelTree.cpp
        BTree/BTree_parallelTree.hpp
        BTree/BTree_genericBaseNode.cpp
        BTree/BTree_genericBaseNode.hpp
        SC_generic.cpp
        SC_generic.hpp)

if(USE_FMM_3D)
    set (SOURCE_FILES ${SOURCE_FILES}
            FMM3D/FMM3D_fmmSolver.cpp
            FMM3D/FMM3D_fmmSolver.hpp)
endif()


add_library(spacecharge STATIC ${SOURCE_FILES})
target_include_directories(spacecharge PUBLIC . BTree FMM3D)
target_link_libraries(spacecharge core)