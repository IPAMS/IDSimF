cmake_minimum_required(VERSION 3.6)
project(IDSimF)

# use CTest extension
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wdelete-non-virtual-dtor -Wall -Wextra -Wsign-conversion -Wfloat-equal -Wfloat-conversion -Wimplicit-fallthrough)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wshorten-64-to-32) #-Wzero-as-null-pointer-constant)
endif()



        #if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#    add_compile_options(-Wreserved-id-macro)
#endif()

# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}") #-fsanitize=address -fno-omit-frame-pointer -O1")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")

option (USE_OpenMP "Use OpenMP" ON)
if(USE_OpenMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        message("OpenMP found")
    endif()
endif()

find_package(HDF5 COMPONENTS CXX HL)
include_directories(${HDF5_INCLUDE_DIRS})

#optional libraries / optional features-------------
option (USE_CPP_FSLIB "Use stdc++fs" OFF)
option (USE_EIGEN "Use Eigen based tests / experimental implementations" OFF)
#-------------------------------

# third party libraries ------------------------------------
add_subdirectory(libs/spdlog)
# add_subdirectory(libs/cxxopts)

add_subdirectory (modules/Core)
add_subdirectory (modules/BTree)
add_subdirectory (modules/RS)
add_subdirectory (modules/backgroundInteraction/collisionModels)
add_subdirectory (modules/particleSimulation)

add_subdirectory (tests/modules/core)
add_subdirectory (tests/modules/BTree)
add_subdirectory (tests/modules/RS)
add_subdirectory (tests/modules/collisionModels)
add_subdirectory (tests/modules/particleSimulation)
add_subdirectory (tests/simulations)
add_subdirectory (tests/applications)

add_subdirectory (tests/benchmarks/errorfunction)
add_subdirectory (tests/benchmarks/filewriter/HDF5filewriter)
add_subdirectory (tests/benchmarks/collision_hardsphere/runtime_benchmark)
add_subdirectory (tests/benchmarks/collision_hardsphere/maxwell_benchmark)
add_subdirectory (tests/benchmarks/collision_sds/runtime_benchmark)
add_subdirectory (tests/benchmarks/RS/runtime_benchmark)
add_subdirectory (tests/benchmarks/simionPA/runtime_benchmark)
add_subdirectory (tests/benchmarks/simionPA/samplingErrors_benchmark)
add_subdirectory (tests/benchmarks/interpolatedField/runtime_benchmark)
add_subdirectory (tests/benchmarks/simulations/simpleSpaceCharge_benchmark)
if (USE_EIGEN)
    add_subdirectory (tests/benchmarks/vector/vectorImplementationsCompare)
endif()

add_subdirectory (applications/util)
add_subdirectory (applications/basic/BT-spaceChargeMinimalSim)
add_subdirectory (applications/basic/BT-spaceChargeMinimalParallelSim)
add_subdirectory (applications/basic/BT-spaceChargeSimpleSim)
add_subdirectory (applications/basic/BT-staticSimionPASim)
add_subdirectory (applications/ionCollision/BT-quadrupoleCollisionCellSim)
add_subdirectory (applications/ionMobility/BT-RS-DMSSim)
add_subdirectory (applications/ionMobility/RS-DMSSimplifiedSim)
add_subdirectory (applications/ionMobility/BT-RS-IMSSim)
#add_subdirectory (applications/ionTransfer/BT-idealizedQuadSim)
add_subdirectory (applications/ionTransfer/BT-generalQuadSim)
add_subdirectory (applications/ionTraps/BT-QITSim)
add_subdirectory (applications/ionTraps/BT-generalTrapSim)
add_subdirectory (applications/ionTraps/BT-LITSim)
add_subdirectory (applications/ionTraps/BT-RS-reactiveQITSim)
add_subdirectory (applications/chemistry/RS-idealIsothermReactorSim)