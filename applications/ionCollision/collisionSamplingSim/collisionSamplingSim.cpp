/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2025 - Physical and Theoretical Chemistry /
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
 collisionSamplingSim.cpp

 Application to study the dynamics of individual ion / background gas particle collisions in detail

 ****************************/
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_inputFileUtilities.hpp"
#include "appUtils_ionDefinitionReading.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "appUtils_commandlineParser.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include "CollisionModel_MDInteractionsTrajectorySampler.hpp"
#include "CollisionModel_MDForceField_LJ12_6.hpp"
#include "CollisionModel_MDForceField_Buckingham.hpp"

int main(int argc, const char * argv[]) {
    try {
        // parse commandline / create conf and logger ===================================================
        AppUtils::CommandlineParser cmdLineParser(argc, argv, "BT-quadrupoleCollisionCellSim",
                "Simulation of a quadrupolar collision cell", true);
        std::string simResultBasename = cmdLineParser.resultName();
        AppUtils::logger_ptr logger = cmdLineParser.logger();
        AppUtils::simConf_ptr simConf = cmdLineParser.simulationConfiguration();

        double collisionGasPolarizability_m3 = simConf->doubleParameter("collision_gas_polarizability_m3");
        std::string collisionGasIdentifier = simConf->stringParameter("collision_gas_identifier");
        std::string particleIdentifier = simConf->stringParameter("particle_identifier");
        double gridSpacing_ang = simConf->doubleParameter("grid_spacing_angstrom");
        int gridSamples = simConf->intParameter("grid_samples");
        double subIntegratorIntegrationTime_s = simConf->doubleParameter("sub_integrator_integration_time_s");
        double subIntegratorStepSize_s = simConf->doubleParameter("sub_integrator_step_size_s");
        double collisionRadiusScaling = simConf->doubleParameter("collision_radius_scaling");
        //double angleThetaScaling = simConf->doubleParameter("angle_theta_scaling");
        //double spawnRadius_m = simConf->doubleParameter("spawn_radius_m");
        //bool saveTrajectory = simConf->boolParameter("save_trajectory");
        double trajectoryMinimalSampleInterval_s = simConf->doubleParameter("trajectory_minimal_sample_interval_s");
        double velocity_x = simConf->doubleParameter("velocity_x");
        //trajectoryDistance_m = simConf->doubleParameter("trajectory_distance_m");
        //saveTrajectoryStartTimeStep = simConf->unsignedIntParameter("trajectory_start_time_step");*/
        std::string potentialsFF = simConf->stringParameter("force_field");
        std::string potentialFunction = simConf->stringParameter("potential_function");

        //read molecular structure file
        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection;
        std::string mdCollisionConfFile = simConf->pathRelativeToConfFile(simConf->stringParameter("md_configuration"));
        FileIO::MolecularStructureReader mdConfReader = FileIO::MolecularStructureReader();
        molecularStructureCollection = mdConfReader.readMolecularStructure(mdCollisionConfFile);

        Core::Particle ion;
        ion.setMolecularStructure(molecularStructureCollection.at(particleIdentifier));
        ion.setVelocity({0, 0, 0});

        double gridSpacing_m = gridSpacing_ang*1e-10;
        for(int i = 0; i < gridSamples; i++) {
            CollisionModel::MDForceField_LJ12_6 forceField(collisionGasPolarizability_m3, potentialsFF);
            auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(forceField);

            Core::Vector particlePosition({-50e-10, i*gridSpacing_m, 0});
            Core::Vector particleVelocity({velocity_x,0,0});
            CollisionModel::MDInteractionsTrajectorySampler mdSim(
                CollisionModel::MDInteractionsTrajectorySampler::DIAMETER_N2,
                "N2",
                subIntegratorIntegrationTime_s,
                subIntegratorStepSize_s,
                collisionRadiusScaling,
                std::move(forceFieldPtr),
                molecularStructureCollection,
                particlePosition,
                particleVelocity);

            mdSim.setTrajectoryWriter(simResultBasename+"_MD_traj.txt", 1.0, 0, trajectoryMinimalSampleInterval_s);
            mdSim.updateModelTimestepParameters(1, 0);
            double dt = 2e-11;
            mdSim.modifyVelocity(ion);
            logger->info("i:{} ",i);
        }

        for(int i = 1; i < gridSamples; i++) {
            CollisionModel::MDForceField_LJ12_6 forceField(collisionGasPolarizability_m3, potentialsFF);
            auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(forceField);

            Core::Vector particlePosition({-50e-10, -i*gridSpacing_m, 0});
            Core::Vector particleVelocity({velocity_x,0,0});
            CollisionModel::MDInteractionsTrajectorySampler mdSim(
                CollisionModel::MDInteractionsTrajectorySampler::DIAMETER_N2,
                "N2",
                subIntegratorIntegrationTime_s,
                subIntegratorStepSize_s,
                collisionRadiusScaling,
                std::move(forceFieldPtr),
                molecularStructureCollection,
                particlePosition,
                particleVelocity);

            mdSim.setTrajectoryWriter(simResultBasename+"_MD_traj.txt", 1.0, 0, trajectoryMinimalSampleInterval_s);
            mdSim.updateModelTimestepParameters(1, 0);
            double dt = 2e-11;
            mdSim.modifyVelocity(ion);
            logger->info("i:{} ",i);
        }
    }
    catch(AppUtils::TerminatedWhileCommandlineParsing& terminatedMessage){
        return terminatedMessage.returnCode();
    }
    catch(const std::invalid_argument& ia){
        std::cout << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
}