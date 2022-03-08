project(test_SpaceCharge)

set(SOURCE_FILES
        test_main.cpp
        BTree/test_node.cpp
        BTree/test_parallelNode.cpp
        BTree/test_parallelTree.cpp
        BTree/test_tree.cpp
        BTree/test_treeParticle.cpp)

if(USE_FMM_3D)
    set(SOURCE_FILES ${SOURCE_FILES}
            FMM3D/test_fmm3DSolver.cpp)
endif()

if(USE_EXAFMMT)
    set(SOURCE_FILES ${SOURCE_FILES}
        ExaFMMt/test_exafmmtSolver.cpp)
endif()

add_executable(test_SpaceCharge ${SOURCE_FILES})
target_include_directories(test_SpaceCharge PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/catch
        ${CMAKE_SOURCE_DIR}/tests/util)

target_link_libraries(test_SpaceCharge core spacecharge particlesimulation)

add_test(NAME test_SpaceCharge COMMAND test_SpaceCharge)