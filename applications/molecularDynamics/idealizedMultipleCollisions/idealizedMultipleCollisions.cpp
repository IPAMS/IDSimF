/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2020 - Physical and Theoretical Chemistry /
 Institute of Pure and Applied Mass Spectrometry
 of the University of Wuppertal, Germany

 IDSimF is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 IDSimF is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with IDSimF.  If not, see <https://www.gnu.org/licenses/>.

 ------------
 idealizedMultipleCollisions.cpp

 Idealized, preconstructed molecular dynamics simulations of ion-neutral collisions for a defined set of start parameters.

 ****************************/

#include "CollisionModel_MDInteractionsPreconstructed.hpp"
#include "CollisionModel_Molecule.hpp"
#include "CollisionModel_Atom.hpp"
#include "Core_randomGenerators.hpp"
#include "Core_constants.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "appUtils_commandlineParser.hpp"
#include "FileIO_CSVReader.hpp"
#include "CollisionModel_MDForceField_Buckingham.hpp"
#include <iostream>

int main(int argc, const char * argv[]) {
    Core::globalRandomGeneratorPool = std::make_unique<Core::RandomGeneratorPool>();

    // open configuration, parse configuration file =========================================
    AppUtils::CommandlineParser cmdLineParser(argc, argv, "idealizedMultipleCollisions", "MD Simulation of a series of ion-neutral collisions", true);
    std::string projectName = cmdLineParser.resultName();
    AppUtils::logger_ptr logger = cmdLineParser.logger();

    std::string confFileName = cmdLineParser.confFileName();
    AppUtils::simConf_ptr simConf = cmdLineParser.simulationConfiguration();

    double backgroundTemperature_K = simConf->doubleParameter("background_temperature_K");
    double backgroundPartialPressures_Pa = simConf->doubleParameter("background_pressures_Pa");
    double collisionGasMasses_Amu = simConf->doubleParameter("collision_gas_masses_amu");
    double collisionGasDiameters_angstrom = simConf->doubleParameter("collision_gas_diameters_angstrom");
    double collisionGasPolarizability_m3 = simConf->doubleParameter("collision_gas_polarizability_m3");
    std::string collisionGasIdentifier = simConf->stringParameter("collision_gas_identifier");
    std::string particleIdentifier = simConf->stringParameter("particle_identifier");
    double subIntegratorIntegrationTime_s = simConf->doubleParameter("sub_integrator_integration_time_s");
    double subIntegratorStepSize_s = simConf->doubleParameter("sub_integrator_step_size_s");
    double collisionRadiusScaling = simConf->doubleParameter("collision_radius_scaling");
    double angleThetaScaling = simConf->doubleParameter("angle_theta_scaling");
    double spawnRadius_m = simConf->doubleParameter("spawn_radius_m");
    std::string trajectoryFile = simConf->stringParameter("trajectory_file");
    double trajectoryDistance_m = simConf->doubleParameter("trajectory_distance_m");
    unsigned int saveTrajectoryStartTimeStep = simConf->unsignedIntParameter("trajectory_start_time_step");
    
    // ======================================================================================

    //read position and velocity input 
    
    std::string startingConfiguration = simConf->pathRelativeToConfFile(simConf->stringParameter("starting_configuration"));
    FileIO::CSVReader startConfReader = FileIO::CSVReader();
    std::vector<std::vector<std::string>> startingConfigurationCollection = startConfReader.readCSVFile(startingConfiguration, ' ');

    //read molecular structure file
    
    std::string mdCollisionConfFile = simConf->pathRelativeToConfFile(simConf->stringParameter("md_configuration"));
    FileIO::MolecularStructureReader mdConfReader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = 
                                                                        mdConfReader.readMolecularStructure(mdCollisionConfFile);


    CollisionModel::MDForceField_Buckingham forceField(collisionGasPolarizability_m3);
    auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_Buckingham>(forceField);

    Core::Particle ion;
    ion.setMolecularStructure(molecularStructureCollection.at(particleIdentifier));

    size_t samples = 1;
    // double pi = 3.1415;
    // std::vector<Core::Vector> ionVelocities;
    // std::vector<Core::Vector> ionPositions;
    // std::vector<Core::Vector> ionRotations;
    // for(int i = 0; i < samples; i++){
    //     Core::Vector tmp1 = Core::Vector(0.0, 0.0, 0.0);
    //     ionVelocities.push_back(tmp1);
    //     ionPositions.push_back(Core::Vector(3.99e-10, 42.5e-10, 0.0));
    //     // Core::Vector tmp2 = Core::Vector(0.0, 0.0, rndSource->uniformRealRndValue()*1*pi/2-pi/4);
    //     Core::Vector tmp2 = Core::Vector(0, 0, 0);
    //     ionRotations.push_back(tmp2);
    // }
    Core::Vector ionVelocities = {0, 0, 0};
    for(size_t i = 0; i < samples; i++){
        ion.setVelocity(ionVelocities);
        CollisionModel::MDInteractionsModelPreconstructed mdSim = 
                                                            CollisionModel::MDInteractionsModelPreconstructed(backgroundPartialPressures_Pa, 
                                                                                                                backgroundTemperature_K, 
                                                                                                                collisionGasMasses_Amu, 
                                                                                                                collisionGasDiameters_angstrom,
                                                                                                                collisionGasIdentifier, 
                                                                                                                subIntegratorIntegrationTime_s, 
                                                                                                                subIntegratorStepSize_s, 
                                                                                                                collisionRadiusScaling, angleThetaScaling, 
                                                                                                                spawnRadius_m,
                                                                                                                std::move(forceFieldPtr),
                                                                                                                molecularStructureCollection);

        
        mdSim.setTrajectoryWriter(trajectoryFile, trajectoryDistance_m, 0);
        mdSim.updateModelTimestepParameters(saveTrajectoryStartTimeStep, 0);
        mdSim.modifyVelocity(ion);
    }
}

