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
#include "AppUtils_simulationConfiguration.hpp"
#include "AppUtils_inputFileUtilities.hpp"
#include "AppUtils_ionDefinitionReading.hpp"
#include "AppUtils_logging.hpp"
#include "AppUtils_stopwatch.hpp"
#include "AppUtils_signalHandler.hpp"
#include "AppUtils_commandlineParser.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include "CollisionModel_MDInteractionsTrajectorySampler.hpp"
#include "CollisionModel_MDForceField_LJ12_6.hpp"
#include "CollisionModel_MDForceField_Buckingham.hpp"
#include "Core_math.hpp"

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
        int maximumSteps = simConf->intParameter("maximum_step_number");
        double trajectoryMinimalSampleInterval_s = simConf->doubleParameter("trajectory_minimal_sample_interval_s");
        double velocity_x = simConf->doubleParameter("velocity_x");
        double angleGas_z_deg = simConf->doubleParameter("gas_particle_angle_z_deg");
        double angleGas_z_rad = Core::degToRad(angleGas_z_deg);
        Core::Vector anglesIon_deg = simConf->vector3dParameter("ion_angles_deg");
        Core::Vector anglesIon_rad = {
            Core::degToRad(anglesIon_deg.x()),
            Core::degToRad(anglesIon_deg.y()),
            Core::degToRad(anglesIon_deg.z()) };
        std::string potentialsFF = simConf->stringParameter("force_field");
        std::string potentialFunction = simConf->stringParameter("potential_function");
        bool ionIsFrozen = simConf->boolParameter("ion_is_frozen");

        //read molecular structure file
        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection;
        std::string mdCollisionConfFile = simConf->pathRelativeToConfFile(simConf->stringParameter("md_configuration"));
        FileIO::MolecularStructureReader mdConfReader = FileIO::MolecularStructureReader();
        molecularStructureCollection = mdConfReader.readMolecularStructure(mdCollisionConfFile);

        Core::Particle ion;
        ion.setMolecularStructure(molecularStructureCollection.at(particleIdentifier));

        std::unique_ptr<CollisionModel::AbstractMDForceField> forceFieldPtr;
        if(potentialFunction == "LJ") {
            forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(collisionGasPolarizability_m3, potentialsFF);
        }
        else if(potentialFunction == "Buckingham") {
            //we have to initialize buckingham force field:
            std::vector<Core::Particle*> particlePtrs = {&ion};
            auto buckinghamPtr = std::make_unique<CollisionModel::MDForceField_Buckingham>(collisionGasPolarizability_m3, potentialsFF);
            buckinghamPtr->populateInteractionTable(particlePtrs, molecularStructureCollection, collisionGasIdentifier);
            forceFieldPtr = std::move(buckinghamPtr);
        }
        CollisionModel::MDInteractionsTrajectorySampler mdSim(
            std::move(forceFieldPtr), molecularStructureCollection, logger);
        mdSim.setTrajectoryWriter(simResultBasename+"_MD_traj.txt", trajectoryMinimalSampleInterval_s);

        double gridSpacing_m = gridSpacing_ang*1e-10;
        for(int i = -gridSamples+1; i < gridSamples; i++) {
            //reset ion position:
            ion.setLocation({0, 0, 0});
            ion.setVelocity({0, 0, 0});

            Core::Vector gasParticlePosition({-50e-10, i*gridSpacing_m, 0});
            Core::Vector gasParticleVelocity({velocity_x,0,0});
            Core::Vector gasParticleRotationAngles({0,0,angleGas_z_rad});
            mdSim.calculateTrajectory(
                ion, anglesIon_rad,
                collisionGasIdentifier, gasParticlePosition, gasParticleVelocity, gasParticleRotationAngles,
                subIntegratorIntegrationTime_s, subIntegratorStepSize_s, maximumSteps, ionIsFrozen);

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
