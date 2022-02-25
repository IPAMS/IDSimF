project(test_BTree)

set(SOURCE_FILES
        test_main.cpp
        test_node.cpp
        test_parallelNode.cpp
        test_parallelTree.cpp
        test_tree.cpp
        test_treeParticle.cpp)

add_executable(test_BTree ${SOURCE_FILES})
target_include_directories(test_BTree PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_BTree core spacecharge particlesimulation)

add_test(NAME test_BTree COMMAND test_BTree)