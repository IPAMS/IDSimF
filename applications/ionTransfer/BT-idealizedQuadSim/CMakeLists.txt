project(BT-idealizedQuadSimulation)
set(SOURCE_FILES
        BT-idealizedQuadSim.cpp)

add_executable(BT-idealizedQuadSimulation ${SOURCE_FILES})
target_link_libraries(BT-idealizedQuadSimulation core apputils btree collisionmodels particlesimulation)
