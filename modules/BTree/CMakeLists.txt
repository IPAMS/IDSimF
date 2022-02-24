project(IDSimF_bTree)

set(SOURCE_FILES
        BTree_abstractNode.cpp
        BTree_abstractNode.hpp
        BTree_node.cpp
        BTree_node.hpp
        BTree_parallelNode.cpp
        BTree_parallelNode.hpp
        BTree_treeParticle.cpp
        BTree_treeParticle.hpp
        BTree_tree.cpp
        BTree_tree.hpp
        BTree_parallelTree.cpp
        BTree_parallelTree.hpp
        BTree_genericBaseNode.cpp
        BTree_genericBaseNode.hpp)

add_library(btree STATIC ${SOURCE_FILES})
target_include_directories(btree PUBLIC .)
target_link_libraries(btree core)