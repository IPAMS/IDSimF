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
 BT-spaceChargeSimpleSim.cpp

 Simple simulation of pure particle / particle interaction (space charge)

 ****************************/

#include "BTree_particle.hpp"
#include "BTree_tree.hpp"
#include "PSim_trajectoryExplorerJSONwriter.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_util.hpp"
#include "PSim_boxStartZone.hpp"
#include "PSim_verletIntegrator.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "appUtils_commandlineParser.hpp"
#include <iostream>
#include <vector>


int main(int argc, const char * argv[]) {
    try {
        // parse commandline / create conf and logger ===================================================
        AppUtils::CommandlineParser cmdLineParser(argc, argv, "BT-spaceChargeSimpleSim",
                "Basic pure space charge simulation (non parallel)", false);
        std::string simResultBasename = cmdLineParser.resultName();
        AppUtils::logger_ptr logger = cmdLineParser.logger();

        std::string confFileName = cmdLineParser.confFileName();
        AppUtils::simConf_ptr simConf = cmdLineParser.simulationConfiguration();


        // read basic simulation parameters =============================================================
        unsigned int timeSteps = simConf->unsignedIntParameter("sim_time_steps");
        int trajectoryWriteInterval = simConf->intParameter("trajectory_write_interval");
        double dt = simConf->doubleParameter("dt");

        //read physical configuration ===================================================================
        double spaceChargeFactor = simConf->doubleParameter("space_charge_factor");


        //read ion configuration =======================================================================
        std::vector<std::unique_ptr<BTree::Particle>> particles;
        std::vector<BTree::Particle*> particlePtrs;

        std::vector<unsigned int> nIons = simConf->unsignedIntVectorParameter("n_ions");
        std::vector<double> ionMasses = simConf->doubleVectorParameter("ion_masses");

        for (std::size_t i = 0; i<nIons.size(); i++) {
            unsigned int nParticles = nIons[i];
            double mass = ionMasses[i];
            ParticleSimulation::BoxStartZone startZone(Core::Vector(3.0, 3.0, 3.0)/1000.0);
            auto ions = startZone.getRandomParticlesInStartZone(nParticles, 1.0);
            /*auto ions = ParticleSimulation::util::getRandomIonsInBox(nParticles, 1.0,
                                                                     Core::Vector(-1.5, -1.5, -1.5) /
                                                                     1000.0,
            Core::Vector(3, 3, 3) / 1000.0);
             */
            for (unsigned int j = 0; j<nParticles; j++) {
                ions[j]->setMassAMU(mass);
                particlePtrs.push_back(ions[j].get());
                particles.push_back(std::move(ions[j]));
            }
        }



        //prepare file writer ==============================================================================

        // function to add some additional exported parameters to the exported trajectory file:
        ParticleSimulation::partAttribTransformFctType additionalParameterTransformFct =
                [](BTree::Particle* particle) -> std::vector<double> {
                    std::vector<double> result = {
                            particle->getVelocity().x(),
                            particle->getVelocity().y(),
                            particle->getVelocity().z()
                    };
                    return result;
                };

        std::vector<std::string> auxParamNames = {"velocity x", "velocity y", "velocity z"};

        auto hdf5Writer = std::make_unique<ParticleSimulation::TrajectoryHDF5Writer>(
                simResultBasename+"_trajectories.hd5");
        hdf5Writer->setParticleAttributes(auxParamNames, additionalParameterTransformFct);

        /*auto jsonWriter = std::make_unique<ParticleSimulation::TrajectoryExplorerJSONwriter>(
                projectName + "_trajectories.json");*/
        //jsonWriter->setScales(1000,1e6);


        // define functions for the trajectory integration ==================================================

        auto accelerationFunction =
                [spaceChargeFactor](
                        BTree::Particle* particle, int /*particleIndex*/,
                        BTree::Tree& tree, double /*time*/, int /*timestep*/) -> Core::Vector {

                    double particleCharge = particle->getCharge();

                    Core::Vector spaceChargeForce(0, 0, 0);
                    if (spaceChargeFactor>0) {
                        spaceChargeForce =
                                tree.computeEFieldFromTree(*particle)*(particleCharge*spaceChargeFactor);
                    }
                    return (spaceChargeForce/particle->getMass());
                };

        auto timestepWriteFunction =
                [trajectoryWriteInterval, &hdf5Writer, &logger](
                        std::vector<BTree::Particle*>& particles, BTree::Tree& /*tree*/, double time,
                        int timestep, bool lastTimestep) {

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
                accelerationFunction, timestepWriteFunction);
        AppUtils::SignalHandler::setReceiver(verletIntegrator);
        verletIntegrator.run(timeSteps, dt);
        stopWatch.stop();

        logger->info("elapsed secs (wall time) {}", stopWatch.elapsedSecondsWall());
        logger->info("elapsed secs (cpu time) {}", stopWatch.elapsedSecondsCPU());
        return EXIT_SUCCESS;
    }
    catch(AppUtils::TerminatedWhileCommandlineParsing& terminatedMessage){
        return terminatedMessage.returnCode();
    }
    catch(const std::invalid_argument& ia){
        std::cout << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
}