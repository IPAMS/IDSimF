project(idealizedQuadSimulation)
set(SOURCE_FILES
        idealizedQuadSim.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} core apputils spacecharge collisionmodels particlesimulation)
