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
 spaceChargeMinimalSim.cpp

 Minimal serial simulation of pure particle / particle interaction (space charge)

 ****************************/

#include "Core_particle.hpp"
#include "BTree_tree.hpp"
#include "FileIO_trajectoryExplorerJSONwriter.hpp"
#include "FileIO_trajectoryHDF5Writer.hpp"
#include "Integration_verletIntegrator.hpp"
#include "FileIO_ionCloudReader.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "appUtils_commandlineParser.hpp"
#include <iostream>
#include <vector>

int main(int argc, const char * argv[]) {
    try{
        // parse commandline / create conf and logger ===================================================
        AppUtils::CommandlineParser cmdLineParser(argc, argv, "BT-spaceChargeMinimalSim",
                "Basic space charge simulation (non parallel, mostly for testing purposes)", false);
        std::string simResultBasename = cmdLineParser.resultName();
        AppUtils::logger_ptr logger = cmdLineParser.logger();

        std::string confFileName = cmdLineParser.confFileName();
        AppUtils::simConf_ptr simConf = cmdLineParser.simulationConfiguration();


        // read basic simulation parameters =============================================================
        unsigned int timeSteps = simConf->unsignedIntParameter("sim_time_steps");
        unsigned int trajectoryWriteInterval = simConf->unsignedIntParameter("trajectory_write_interval");
        double dt = simConf->doubleParameter("dt");

        //read physical configuration ===================================================================
        double spaceChargeFactor = simConf->doubleParameter("space_charge_factor");


        //read ion configuration =======================================================================
        std::vector<std::unique_ptr<Core::Particle>> particles;
        std::vector<Core::Particle*> particlePtrs;

        if (simConf->isParameter("ion_cloud_init_file")) {
            std::string ionCloudFileName = simConf->pathRelativeToConfFile(
                    simConf->stringParameter("ion_cloud_init_file"));
            FileIO::IonCloudReader reader = FileIO::IonCloudReader();
            particles = reader.readIonCloud(ionCloudFileName);
            //prepare a vector of raw pointers
            for (const auto& part : particles) {
                particlePtrs.push_back(part.get());
            }
        }
        else {
            throw std::invalid_argument("missing configuration value: ion_cloud_init_file");
        }

        //prepare file writer ==============================================================================

        // function to add some additional exported parameters to the exported trajectory file:
        FileIO::partAttribTransformFctType additionalParameterTransformFct =
                [](Core::Particle* particle) -> std::vector<double> {
                    std::vector<double> result = {
                            particle->getVelocity().x(),
                            particle->getVelocity().y(),
                            particle->getVelocity().z()
                    };
                    return result;
                };

        std::vector<std::string> auxParamNames = {"velocity x", "velocity y", "velocity z"};

        auto hdf5Writer = std::make_unique<FileIO::TrajectoryHDF5Writer>(simResultBasename+"_trajectories.hd5");
        hdf5Writer->setParticleAttributes(auxParamNames, additionalParameterTransformFct);

        /*auto jsonWriter = std::make_unique<FileIO::TrajectoryExplorerJSONwriter>(
                projectName + "_trajectories.json");*/
        //jsonWriter->setScales(1000,1e6);


        // define functions for the trajectory integration ==================================================

        auto accelerationFunction =
                [spaceChargeFactor](
                        Core::Particle* particle, int /*particleIndex*/,
                        SpaceCharge::FieldCalculator& scFieldCalculator, double /*time*/, int /*timestep*/) -> Core::Vector {

                    double particleCharge = particle->getCharge();

                    Core::Vector spaceChargeForce(0, 0, 0);
                    if (spaceChargeFactor>0) {
                        spaceChargeForce =
                                scFieldCalculator.getEFieldFromSpaceCharge(*particle)*(particleCharge*spaceChargeFactor);
                    }
                    return (spaceChargeForce/particle->getMass());
                };

        auto postTimestepFunction =
                [trajectoryWriteInterval, &hdf5Writer, &logger](
                        Integration::AbstractTimeIntegrator* /*integrator*/,
                        std::vector<Core::Particle*>& particles, double time,
                        unsigned int timestep, bool lastTimestep)
                {

                    if (lastTimestep) {
                        hdf5Writer->writeTimestep(particles, time);

                        hdf5Writer->writeSplatTimes(particles);
                        hdf5Writer->finalizeTrajectory();
                        logger->info("finished ts:{} time:{:.2e}", timestep, time);
                    }

                    else if (timestep%trajectoryWriteInterval==0) {
                        logger->info("ts:{} time:{:.2e}", timestep, time);
                        hdf5Writer->writeTimestep(particles, time);
                    }
                };

        // simulate ===============================================================================================
        AppUtils::Stopwatch stopWatch;
        stopWatch.start();
        Integration::VerletIntegrator verletIntegrator(
                particlePtrs,
                accelerationFunction, postTimestepFunction);
        AppUtils::SignalHandler::setReceiver(verletIntegrator);
        verletIntegrator.run(timeSteps, dt);

        logger->info("elapsed secs (wall time) {}", stopWatch.elapsedSecondsWall());
        logger->info("elapsed secs (cpu time) {}", stopWatch.elapsedSecondsCPU());
        return EXIT_SUCCESS;
    }
    catch(AppUtils::TerminatedWhileCommandlineParsing& terminatedMessage){
        return terminatedMessage.returnCode();
    }
    catch(const FileIO::IonCloudFileException& ie)
    {
        std::cout << ie.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch(const std::invalid_argument& ia){
        std::cout << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
}