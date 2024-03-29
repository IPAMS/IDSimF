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
        SC_generic.hpp
        SC_fullSumSolver.cpp
        SC_fullSumSolver.hpp)

if(USE_FMM_3D)
    set (SOURCE_FILES ${SOURCE_FILES}
            FMM3D/FMM3D_fmmSolver.cpp
            FMM3D/FMM3D_fmmSolver.hpp
            FMM3D/FMM3D_fmmSolver_C_interface.c
            FMM3D/FMM3D_fmmSolver_C_interface.h)
endif()

if(USE_EXAFMMT)
    set (SOURCE_FILES ${SOURCE_FILES}
            ExaFMMt/ExaFMMt_fmmSolver.cpp
            ExaFMMt/ExaFMMt_fmmSolver.hpp)
endif()


add_library(spacecharge STATIC ${SOURCE_FILES})
target_include_directories(spacecharge PUBLIC . BTree FMM3D ExaFMMt)
target_link_libraries(spacecharge core)
if(OpenMP_CXX_FOUND)
    target_link_libraries(spacecharge OpenMP::OpenMP_CXX)
endif()

if(USE_FMM_3D)
    target_include_directories(spacecharge PUBLIC ${FMM_3D_PATH}/c)
    target_link_libraries(spacecharge ${FMM_3D_PATH}/lib/libfmm3d.so)
endif()

if(USE_EXAFMMT)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/ExaFMMt/cmake/findFFTW")

    message(STATUS ${CMAKE_MODULE_PATH})
    target_include_directories(spacecharge PUBLIC ${EXAFMMT_PATH}/include)

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    find_package(FFTW REQUIRED)

    target_link_libraries (spacecharge
            ${BLAS_LIBRARIES}
            ${LAPACK_LIBRARIES}
            ${FFTW_LIBRARIES})
endif()
