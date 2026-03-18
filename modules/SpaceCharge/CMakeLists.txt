project(IDSimF_spaceCharge)

set(SOURCE_FILES
        BTree/BTree_abstractNode.cpp
        BTree/BTree_node.cpp
        BTree/BTree_parallelNode.cpp
        BTree/BTree_treeParticle.cpp
        BTree/BTree_tree.cpp
        BTree/BTree_parallelTree.cpp
        BTree/BTree_genericBaseNode.cpp
        SC_generic.cpp
        SC_fullSumSolver.cpp)

if(USE_FMM_3D)
    set (SOURCE_FILES ${SOURCE_FILES}
            FMM3D/FMM3D_fmmSolver.cpp
            FMM3D/FMM3D_fmmSolver_C_interface.c
            FMM3D/FMM3D_fmmSolver_C_interface.h)
endif()

if(USE_EXAFMMT)
    set (SOURCE_FILES ${SOURCE_FILES}
            ExaFMMt/ExaFMMt_fmmSolver.cpp)
endif()


add_library(spacecharge STATIC ${SOURCE_FILES})
target_include_directories(spacecharge PUBLIC . BTree FMM3D ExaFMMt)
target_link_libraries(spacecharge PRIVATE core)
if(OpenMP_CXX_FOUND)
    target_link_libraries(spacecharge PRIVATE OpenMP::OpenMP_CXX)
endif()

if(USE_FMM_3D)
    target_include_directories(spacecharge PUBLIC ${FMM_3D_PATH}/c)
    target_link_libraries(spacecharge PRIVATE ${FMM_3D_PATH}/lib/libfmm3d.so)
endif()

if(USE_EXAFMMT)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/ExaFMMt/cmake/findFFTW")

    message(STATUS ${CMAKE_MODULE_PATH})
    target_include_directories(spacecharge PUBLIC ${EXAFMMT_PATH}/include ${FFTW_INCLUDE_DIRS})

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    find_package(FFTW REQUIRED)

    target_link_libraries (spacecharge PRIVATE
            ${BLAS_LIBRARIES}
            ${LAPACK_LIBRARIES}
            ${FFTW_LIBRARIES})
endif()
