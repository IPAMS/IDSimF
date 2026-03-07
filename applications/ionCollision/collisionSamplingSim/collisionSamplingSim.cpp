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
#include "FileIO_CSVReader.hpp"
#include "CollisionModel_MDInteractionsTrajectorySampler.hpp"
#include "CollisionModel_MDForceField_LJ12_6.hpp"
#include "CollisionModel_MDForceField_Buckingham.hpp"
#include "Core_math.hpp"

enum SamplerMode {
    yGrid, ///< single grid line on y axis
    initFile ///< file with initial conditions
};

struct InitialCondition {
    Core::Vector startPosition;
    Core::Vector startVelocity;
    Core::Vector startRotation;
};

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
        double subIntegratorIntegrationTime_s = simConf->doubleParameter("sub_integrator_integration_time_s");
        double subIntegratorStepSize_s = simConf->doubleParameter("sub_integrator_step_size_s");
        int maximumSteps = simConf->intParameter("maximum_step_number");
        double trajectoryMinimalSampleInterval_s = simConf->doubleParameter("trajectory_minimal_sample_interval_s");
        std::string potentialsFF = simConf->stringParameter("force_field");
        std::string potentialFunction = simConf->stringParameter("potential_function");
        bool rotationOn = simConf->boolParameter("rotation");
        bool ionIsFrozen = simConf->boolParameter("ion_is_frozen");
        bool useHDF5Writer = simConf->boolParameter("hdf5_writer");

        std::string samplerModeStr = simConf->stringParameter("sampler_mode");

        SamplerMode samplerMode;
        if (samplerModeStr == "y_grid") {
            samplerMode = yGrid;
        }
        else if (samplerModeStr == "init_file") {
            samplerMode = initFile;
        }
        else {
            throw(std::invalid_argument("Illegal sampler mode"));
        }


        //read molecular structure file
        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection;
        std::string mdCollisionConfFile = simConf->pathRelativeToConfFile(simConf->stringParameter("md_configuration"));
        FileIO::MolecularStructureReader mdConfReader = FileIO::MolecularStructureReader();
        molecularStructureCollection = mdConfReader.readMolecularStructure(mdCollisionConfFile);

        Core::Particle ion;
        ion.setMolecularStructure(molecularStructureCollection.at(particleIdentifier));

        std::unique_ptr<CollisionModel::AbstractMDForceField> forceFieldPtr;
        if(potentialFunction == "LJ") {
            forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(collisionGasPolarizability_m3, potentialsFF, rotationOn);
        }
        else if(potentialFunction == "Buckingham") {
            //we have to initialize buckingham force field:
            std::vector<Core::Particle*> particlePtrs = {&ion};
            auto buckinghamPtr = std::make_unique<CollisionModel::MDForceField_Buckingham>(collisionGasPolarizability_m3, potentialsFF, rotationOn);
            buckinghamPtr->populateInteractionTable(particlePtrs, molecularStructureCollection, collisionGasIdentifier);
            forceFieldPtr = std::move(buckinghamPtr);
        }
        CollisionModel::MDInteractionsTrajectorySampler mdSim(
            std::move(forceFieldPtr), rotationOn, molecularStructureCollection, logger);

        if (useHDF5Writer) {
            mdSim.setHDF5TrajectoryWriter(simResultBasename+"_MD_traj.h5",10, 0);
        }
        else {
            mdSim.setLegacyTrajectoryWriter(simResultBasename+"_MD_traj.txt", 10, trajectoryMinimalSampleInterval_s, 0);
        }


        Core::Vector anglesIon_deg = simConf->vector3dParameter("ion_angles_deg");
        Core::Vector anglesIon_rad = {
            Core::degToRad(anglesIon_deg.x()),
            Core::degToRad(anglesIon_deg.y()),
            Core::degToRad(anglesIon_deg.z()) };

        std::vector<InitialCondition> initialConditions;

        if (samplerMode == yGrid) {
            double gridSpacing_ang = simConf->doubleParameter("grid_spacing_angstrom");
            int gridSamples = simConf->intParameter("grid_samples");
            double velocity_x = simConf->doubleParameter("velocity_x");
            Core::Vector anglesGasParticle_deg = simConf->vector3dParameter("gas_particle_angles_deg");
            Core::Vector anglesGasParticle_rad = {
                Core::degToRad(anglesGasParticle_deg.x()),
                Core::degToRad(anglesGasParticle_deg.y()),
                Core::degToRad(anglesGasParticle_deg.z()) };

            double gridSpacing_m = gridSpacing_ang*1e-10;
            for(int i = -gridSamples+1; i < gridSamples; i++) {
                Core::Vector gasParticlePosition({-50e-10, i*gridSpacing_m, 0});
                Core::Vector gasParticleVelocity({velocity_x,0,0});
                Core::Vector gasParticleRotationAngles(anglesGasParticle_rad);

                initialConditions.emplace_back(
                    InitialCondition{gasParticlePosition, gasParticleVelocity, gasParticleRotationAngles}
                    );
            }
        }
        else if (samplerMode == initFile) {
            //read init file, with initial conditions
            FileIO::CSVReader initFileReader = FileIO::CSVReader();
            std::vector<std::vector<std::string>> stringVector = std::vector<std::vector<std::string>>();
            std::string initFileFn = simConf->pathRelativeToConfFile(simConf->stringParameter("init_file"));
            stringVector = initFileReader.readCSVFile(initFileFn, ';');

            std::vector<double> startPosX = initFileReader.extractDouble(stringVector, 0);
            std::vector<double> startPosY = initFileReader.extractDouble(stringVector, 1);
            std::vector<double> startPosZ = initFileReader.extractDouble(stringVector, 2);

            std::vector<double> startVelocityX = initFileReader.extractDouble(stringVector, 3);
            std::vector<double> startVelocityY = initFileReader.extractDouble(stringVector, 4);
            std::vector<double> startVelocityZ = initFileReader.extractDouble(stringVector, 5);

            for (size_t i=0; i<stringVector.size(); i++) {
                InitialCondition initCon;
                initCon.startPosition = {startPosX[i], startPosY[i], startPosZ[i]};
                initCon.startVelocity = {startVelocityX[i], startVelocityY[i], startVelocityZ[i]};
                initCon.startRotation = {0.0, 0.0, 0.0};
                initialConditions.emplace_back(initCon);
            }
        }

        for(size_t i = 0; i<initialConditions.size(); i++) {
            //reset ion position:
            ion.setLocation({0, 0, 0});
            ion.setVelocity({0, 0, 0});

            InitialCondition initCon = initialConditions.at(i);

            mdSim.calculateTrajectory(
                ion, anglesIon_rad,
                collisionGasIdentifier, initCon.startPosition, initCon.startVelocity, initCon.startRotation,
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
