cmake_minimum_required(VERSION 3.9)
project(IDSimF)

# use CTest extension
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wdelete-non-virtual-dtor -Wall -Wextra -Wsign-conversion -Wfloat-equal -Wfloat-conversion -Wimplicit-fallthrough -Wpedantic)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wshorten-64-to-32 -stdlib=libc++) #-Wzero-as-null-pointer-constant)
    add_link_options(-stdlib=libc++)
endif()

## get number of cores for parallel test runs:
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
    #set(CTEST_BUILD_FLAGS -j${N})
    set(N_THREADS ${N})
endif()

        #if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#    add_compile_options(-Wreserved-id-macro)
#endif()

# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}") #-fsanitize=address -fno-omit-frame-pointer -O1")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

option (USE_OpenMP "Use OpenMP" ON)
if(USE_OpenMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        message("OpenMP CXX found")
    endif()
endif()

find_package(HDF5 COMPONENTS CXX HL)
include_directories(${HDF5_INCLUDE_DIRS})

#optional libraries / optional features-------------
option (USE_CPP_FSLIB "Use stdc++fs" OFF)
option (USE_EIGEN "Use Eigen based tests / experimental implementations" OFF)

## FMM3d library ----
if(NOT DEFINED FMM_3D_PATH)
    set(USE_FMM_3D false)
    message(STATUS "FMM 3d is FALSE")
else()
    set(USE_FMM_3D true)
    message(STATUS "FMM 3d is TRUE, path is:")
    message(STATUS ${FMM_3D_PATH})
    add_compile_definitions(WITH_FMM_3d)
endif()

## exafmm-t library ----
if(NOT DEFINED EXAFMMT_PATH)
    set(USE_EXAFMMT false)
    message(STATUS "ExaFMM-t is FALSE")
else()
    set(USE_EXAFMMT true)
    message(STATUS "ExaFMM-t is TRUE, path is:")
    message(STATUS ${EXAFMMT_PATH})
    add_compile_definitions(WITH_EXAFMMT)
endif()


#-------------------------------

# third party libraries ------------------------------------
add_subdirectory(libs/spdlog)
# add_subdirectory(libs/cxxopts)

add_subdirectory (modules/Core)
add_subdirectory (modules/SpaceCharge)
add_subdirectory (modules/RS)
add_subdirectory (modules/BackgroundInteraction/CollisionModels)
add_subdirectory (modules/FileIO)
add_subdirectory (modules/Integration)
add_subdirectory (modules/ParticleSimulation)

add_subdirectory (tests/modules/Core)
add_subdirectory (tests/modules/SpaceCharge)
add_subdirectory (tests/modules/RS)
add_subdirectory (tests/modules/CollisionModels)
add_subdirectory (tests/modules/ParticleSimulation)
add_subdirectory (tests/modules/Integration)
add_subdirectory (tests/modules/FileIO)
add_subdirectory (tests/simulations)
add_subdirectory (tests/applications)

add_subdirectory (tests/benchmarks/errorfunction)
add_subdirectory (tests/benchmarks/filewriter/HDF5filewriter)
add_subdirectory (tests/benchmarks/collision_hardsphere/runtime_benchmark)
add_subdirectory (tests/benchmarks/collision_hardsphere/maxwell_benchmark)
add_subdirectory (tests/benchmarks/collision_sds/runtime_benchmark)
add_subdirectory (tests/benchmarks/collision_md/runtime_benchmark)
add_subdirectory (tests/benchmarks/RS/runtime_benchmark)
add_subdirectory (tests/benchmarks/simionPA/runtime_benchmark)
add_subdirectory (tests/benchmarks/simionPA/samplingErrors_benchmark)
add_subdirectory (tests/benchmarks/interpolatedField/runtime_benchmark)
add_subdirectory (tests/benchmarks/simulations/simpleSpaceCharge_benchmark)
add_subdirectory (tests/benchmarks/parallelization/rng_parallelization)
if (USE_EIGEN)
    add_subdirectory (tests/benchmarks/vector/vectorImplementationsCompare)
endif()

add_subdirectory (applications/util)
add_subdirectory (applications/basic/spaceChargeMinimalSim)
add_subdirectory (applications/basic/spaceChargeMinimalParallelSim)
add_subdirectory (applications/basic/spaceChargeSimpleSim)
add_subdirectory (applications/basic/staticSimionPASim)
add_subdirectory (applications/ionCollision/quadrupoleCollisionCellSim)
add_subdirectory (applications/ionMobility/util)
add_subdirectory (applications/ionMobility/DMSSim)
add_subdirectory (applications/ionMobility/DMSSimplifiedSim)
add_subdirectory (applications/ionMobility/IMSSim)
#add_subdirectory (applications/ionTransfer/BT-idealizedQuadSim)
add_subdirectory (applications/ionTransfer/generalQuadSim)
add_subdirectory (applications/ionTraps/QITSim)
add_subdirectory (applications/ionTraps/generalTrapSim)
add_subdirectory (applications/ionTraps/LITSim)
add_subdirectory (applications/ionTraps/reactiveQITSim)
add_subdirectory (applications/chemistry/idealIsothermReactorSim)