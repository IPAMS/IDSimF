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
 BT-RS-IMSSim.cpp

 Isothermic continuous field ion mobility spectrometry (IMS) transport and chemistry simulation

 ****************************/


#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_ConcentrationFileWriter.hpp"
#include "PSim_util.hpp"
#include "PSim_constants.hpp"
#include "PSim_verletIntegrator.hpp"
#include "PSim_velocityIntegrator.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "CollisionModel_StatisticalDiffusion.hpp"
#include "CollisionModel_MultiCollisionModel.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "json.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

enum IntegratorType{
    VERLET, SIMPLE, NO_INTEGRATOR
};
enum CollisionModelType{
    HS, SDS, NO_COLLISONS
};

std::string key_ChemicalIndex = "keyChemicalIndex";


int main(int argc, const char *argv[]){

    try {
        // open configuration, parse configuration file =========================================
        if (argc<=2) {
            std::cout << "Run abort: No run configuration or project name given." << std::endl;
            return EXIT_FAILURE;
        }
        std::string projectName = argv[2];
        std::cout << projectName << std::endl;
        auto logger = AppUtils::createLogger(projectName+".log");

        std::string confFileName = argv[1];
        AppUtils::SimulationConfiguration simConf(confFileName, logger);

        std::vector<int> nParticles = simConf.intVectorParameter("n_particles");
        int nSteps = simConf.intParameter("sim_time_steps");
        int concentrationWriteInterval = simConf.intParameter("concentrations_write_interval");
        int trajectoryWriteInterval = simConf.intParameter("trajectory_write_interval");
        bool writeVelocities = simConf.boolParameter("trajectory_write_velocities");
        double dt_s = simConf.doubleParameter("dt_s");
        double eFieldMagnitude = simConf.doubleParameter("electric_field_mag_Vm-1");
        double spaceChargeFactor = simConf.doubleParameter("space_charge_factor");

        double startWidthX_m = simConf.doubleParameter("start_width_x_mm")/1000.0;
        double startWidthYZ_m = simConf.doubleParameter("start_width_yz_mm")/1000.0;
        double stopPosX_m = simConf.doubleParameter("stop_position_x_mm")/1000.0;

        //read and check gas parameters:
        std::string transportModelType = simConf.stringParameter("transport_model_type");
        double backgroundTemperature_K = simConf.doubleParameter("background_temperature_K");

        std::vector<double> backgroundPartialPressures_Pa = simConf.doubleVectorParameter(
                "background_partial_pressures_Pa");
        std::vector<double> collisionGasMasses_Amu = simConf.doubleVectorParameter("collision_gas_masses_amu");
        std::vector<double> collisionGasDiameters_angstrom = simConf.doubleVectorParameter(
                "collision_gas_diameters_angstrom");

        std::size_t nBackgroundGases = backgroundPartialPressures_Pa.size();
        if (collisionGasMasses_Amu.size()!=nBackgroundGases || collisionGasDiameters_angstrom.size()!=nBackgroundGases) {
            throw std::invalid_argument("Inconsistent background gas configuration");
        }

        //compute additional gas parameters:
        double totalBackgroundPressure_Pa = std::accumulate(
                backgroundPartialPressures_Pa.begin(),
                backgroundPartialPressures_Pa.end(), 0.0);

        std::vector<double> collisionGasDiameters_m;
        std::transform(
                collisionGasDiameters_angstrom.begin(),
                collisionGasDiameters_angstrom.end(),
                std::back_inserter(collisionGasDiameters_m),
                [](double cgd) -> double { return cgd*1e-10; });

        // ======================================================================================

        //read and prepare chemical configuration ===============================================
        std::string rsConfFileName = simConf.pathRelativeToConfFile(simConf.stringParameter("reaction_configuration"));
        RS::ConfigFileParser parser = RS::ConfigFileParser();
        RS::Simulation rsSim = RS::Simulation(parser.parseFile(rsConfFileName));
        RS::SimulationConfiguration* rsSimConf = rsSim.simulationConfiguration();
        //prepare a map for retrieval of the substance index:
        std::map<RS::Substance*, int> substanceIndices;
        std::vector<RS::Substance*> discreteSubstances = rsSimConf->getAllDiscreteSubstances();
        std::vector<double> ionMobility; // = simConf.doubleVectorParameter("ion_mobility",confRoot);
        for (std::size_t i = 0; i<discreteSubstances.size(); i++) {
            substanceIndices.insert(std::pair<RS::Substance*, int>(discreteSubstances[i], i));
            ionMobility.push_back(discreteSubstances[i]->mobility());
        }

        // prepare file writer  =================================================================
        RS::ConcentrationFileWriter resultFilewriter(projectName+"_conc.csv");

        //prepare auxiliary parameters transform functions
        ParticleSimulation::partAttribTransformFctType additionalParamTFct;
        std::vector<std::string> auxParamNames;

        if (writeVelocities) {
            additionalParamTFct = [](BTree::Particle* particle) -> std::vector<double> {
                std::vector<double> result = {
                        particle->getFloatAttribute(key_ChemicalIndex),
                        particle->getVelocity().x(),
                        particle->getVelocity().y(),
                        particle->getVelocity().z()
                };
                return result;
            };
            auxParamNames = {"chemical id", "velocity x", "velocity y", "velocity z"};
        }
        else {
            additionalParamTFct = [](BTree::Particle* particle) -> std::vector<double> {
                std::vector<double> result = {
                        particle->getFloatAttribute(key_ChemicalIndex)
                };
                return result;
            };
            auxParamNames = {"chemical id"};
        }

        //init hdf5 filewriter
        std::string hdf5Filename = projectName+"_trajectories.hd5";
        ParticleSimulation::TrajectoryHDF5Writer hdf5Writer(hdf5Filename);
        hdf5Writer.setParticleAttributes(auxParamNames, additionalParamTFct);

        int ionsInactive = 0;

        //fixme: nAllParticles and nTotalParticles are the same parameter
        int nAllParticles = 0;
        for (const auto ni: nParticles) {
            nAllParticles += ni;
        }

        // init simulation  =====================================================================

        // create and add simulation particles:
        int nParticlesTotal = 0;
        std::vector<uniqueReactivePartPtr> particles;
        std::vector<BTree::Particle*> particlesPtrs;
        std::vector<std::vector<double>> trajectoryAdditionalParams;

        Core::Vector initCorner(0, 0, 0);
        Core::Vector initBoxSize(startWidthX_m, startWidthYZ_m, startWidthYZ_m);


        for (std::size_t i = 0; i<nParticles.size(); i++) {
            RS::Substance* subst = rsSimConf->substance(i);
            std::vector<Core::Vector> initialPositions =
                    ParticleSimulation::util::getRandomPositionsInBox(nParticles[i], initCorner, initBoxSize);
            for (int k = 0; k<nParticles[i]; k++) {
                uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(subst);

                particle->setLocation(initialPositions[k]);
                particlesPtrs.push_back(particle.get());
                rsSim.addParticle(particle.get(), nParticlesTotal);
                particles.push_back(std::move(particle));
                trajectoryAdditionalParams.emplace_back(std::vector<double>(1));
                nParticlesTotal++;
            }
        }

        RS::ReactionConditions reactionConditions = RS::ReactionConditions();
        reactionConditions.temperature = backgroundTemperature_K;
        reactionConditions.pressure = totalBackgroundPressure_Pa;
        reactionConditions.electricField = eFieldMagnitude;

        resultFilewriter.initFile(rsSimConf);
        // ======================================================================================

        //check which integrator type we have to setup:
        std::vector<std::string> verletTypes{"btree_SDS", "btree_HS"};
        auto vType = std::find(std::begin(verletTypes), std::end(verletTypes), transportModelType);

        IntegratorType integratorType;
        if (vType!=std::end(verletTypes)) {
            integratorType = VERLET;
            logger->info("Verlet type simulation");
        }
        else if (transportModelType=="simple") {
            integratorType = SIMPLE;
            logger->info("Simple transport simulation");
        }
        else if (transportModelType=="no_transport") {
            integratorType = NO_INTEGRATOR;
            logger->info("No transport simulation");
        }
        else {
            throw std::invalid_argument("illegal transport simulation type");
        }
        //===================================================================================================

        // define trajectory integration parameters / functions =================================
        double referencePressure_Pa = 100000;
        double referenceTemperature_K = 273.15;
        double backgroundPTRatio =
                referencePressure_Pa/totalBackgroundPressure_Pa*backgroundTemperature_K/referenceTemperature_K;

        auto accelerationFctVerlet =
                [eFieldMagnitude, spaceChargeFactor]
                        (BTree::Particle* particle, int /*particleIndex*/, BTree::Tree& tree, double /*time*/, int /*timestep*/) {
                    double particleCharge = particle->getCharge();

                    Core::Vector fieldForce(eFieldMagnitude*particleCharge, 0, 0);

                    if (spaceChargeFactor==0.0) {
                        return (fieldForce/particle->getMass());
                    }
                    else {
                        Core::Vector spaceChargeForce =
                                tree.computeEFieldFromTree(*particle)*(particleCharge*spaceChargeFactor);
                        return ((fieldForce+spaceChargeForce)/particle->getMass());
                    }
                };

        auto timestepWriteFctSimple =
                [&hdf5Writer, trajectoryWriteInterval, &ionsInactive, &logger]
                        (std::vector<BTree::Particle*>& particles, double time, int timestep, bool lastTimestep) {
                    if (lastTimestep) {
                        hdf5Writer.writeTimestep(particles, time);
                        hdf5Writer.writeSplatTimes(particles);
                        hdf5Writer.finalizeTrajectory();
                        logger->info("finished ts:{} time:{:.2e}", timestep, time);
                    }
                    else if (timestep%trajectoryWriteInterval==0) {
                        hdf5Writer.writeTimestep(particles, time);
                        logger->info("ts:{} time:{:.2e} splatted ions:{}", timestep, time, ionsInactive);
                    }
                };

        auto timestepWriteFctVerlet =
                [&timestepWriteFctSimple]
                        (std::vector<BTree::Particle*>& particles, BTree::Tree& /*tree*/, double time, int timestep,
                         bool lastTimestep) {
                    timestepWriteFctSimple(particles, time, timestep, lastTimestep);
                };

        auto otherActionsFunctionIMSSimple =
                [stopPosX_m, &ionsInactive]
                        (Core::Vector& newPartPos, BTree::Particle* particle, int /*particleIndex*/, double time,
                         int /*timestep*/) {
                    if (newPartPos.x()>=stopPosX_m) {
                        particle->setActive(false);
                        particle->setSplatTime(time);
                        ionsInactive++;
                    }
                };

        auto otherActionsFunctionIMSVerlet =
                [&otherActionsFunctionIMSSimple]
                        (Core::Vector& newPartPos, BTree::Particle* particle, int particleIndex,
                         BTree::Tree& /*tree*/, double time, int timestep) {
                    otherActionsFunctionIMSSimple(newPartPos, particle, particleIndex, time, timestep);
                };


        //define and init transport models and trajectory integrators:
        std::unique_ptr<CollisionModel::AbstractCollisionModel> collisionModelPtr;
        CollisionModelType collisionModelType = NO_COLLISONS;
        if (transportModelType=="btree_SDS") {
            //prepare SDS model
            if (nBackgroundGases!=1) {
                throw std::invalid_argument("SDS simulation requires a single collision gas");
            }

            //create sds collision model, if statistics file is given: create with custom statistics
            std::unique_ptr<CollisionModel::StatisticalDiffusionModel> collisionModel;
            if (simConf.isParameter("sds_collision_statistics")) {
                std::string sdsCollisionStatisticsFileName = simConf.pathRelativeToConfFile(
                        simConf.stringParameter("sds_collision_statistics"));

                logger->info("SDS with custom collision statistics file: {}", sdsCollisionStatisticsFileName);
                CollisionModel::CollisionStatistics cs(sdsCollisionStatisticsFileName);
                collisionModel =
                        std::make_unique<CollisionModel::StatisticalDiffusionModel>(
                                backgroundPartialPressures_Pa[0],
                                backgroundTemperature_K,
                                collisionGasMasses_Amu[0],
                                collisionGasDiameters_m[0],
                                cs
                        );
            }
            else {
                collisionModel =
                        std::make_unique<CollisionModel::StatisticalDiffusionModel>(
                                backgroundPartialPressures_Pa[0],
                                backgroundTemperature_K,
                                collisionGasMasses_Amu[0],
                                collisionGasDiameters_m[0]);
            }

            for (const auto& particle: particlesPtrs) {
                collisionModel->setSTPParameters(*particle);
            }
            collisionModelPtr = std::move(collisionModel);
            collisionModelType = SDS;
        }
        else if (transportModelType=="btree_HS") {
            //prepare multimodel with multiple Hard Sphere models (one per collision gas)
            std::vector<std::unique_ptr<CollisionModel::AbstractCollisionModel>> hsModels;
            for (std::size_t i = 0; i<nBackgroundGases; ++i) {
                auto hsModel = std::make_unique<CollisionModel::HardSphereModel>(
                        backgroundPartialPressures_Pa[i],
                        backgroundTemperature_K,
                        collisionGasMasses_Amu[i],
                        collisionGasDiameters_m[i]);
                hsModels.emplace_back(std::move(hsModel));
            }

            std::unique_ptr<CollisionModel::MultiCollisionModel> collisionModel =
                    std::make_unique<CollisionModel::MultiCollisionModel>(std::move(hsModels));

            collisionModelPtr = std::move(collisionModel);
            collisionModelType = HS;
        }

        //init trajectory simulation object:
        std::unique_ptr<ParticleSimulation::AbstractTimeIntegrator> trajectoryIntegrator = nullptr;
        if (integratorType==VERLET) {
            trajectoryIntegrator = std::make_unique<ParticleSimulation::VerletIntegrator>(
                    particlesPtrs,
                    accelerationFctVerlet, timestepWriteFctVerlet, otherActionsFunctionIMSVerlet,
                    ParticleSimulation::noFunction,
                    collisionModelPtr.get());
        }
        else if (integratorType==SIMPLE) {

            auto velocityFctSimple = [eFieldMagnitude, backgroundPTRatio](BTree::Particle* particle, int /*particleIndex*/,
                                                                          double /*time*/, int /*timestep*/) {
                double particleMobility = particle->getMobility();

                Core::Vector velocity(eFieldMagnitude*particleMobility*backgroundPTRatio, 0, 0);
                return velocity;
            };

            trajectoryIntegrator = std::make_unique<ParticleSimulation::VelocityIntegrator>(
                    particlesPtrs,
                    velocityFctSimple, timestepWriteFctSimple, otherActionsFunctionIMSSimple
            );
        }
        // ======================================================================================


        // simulate   ===========================================================================
        AppUtils::SignalHandler::registerSignalHandler(); //this method is used because trajectory integrator can be null
        AppUtils::Stopwatch stopWatch;
        stopWatch.start();

        for (int step = 0; step<nSteps; step++) {
            if (step%concentrationWriteInterval==0) {
                resultFilewriter.writeTimestep(rsSim);
            }
            if (step%trajectoryWriteInterval==0) {
                rsSim.logConcentrations(logger);
            }
            for (int i = 0; i<nParticlesTotal; i++) {
                bool reacted = rsSim.react(i, reactionConditions, dt_s);
                int substIndex = substanceIndices.at(particles[i]->getSpecies());
                particles[i]->setFloatAttribute(key_ChemicalIndex, substIndex);

                if (reacted && collisionModelType==SDS) {
                    //we had an reaction event: update the collision model parameters for the particle which are not
                    //based on location (mostly STP parameters in SDS)
                    collisionModelPtr->initializeModelParameters(*particles[i]);
                }
            }
            rsSim.advanceTimestep(dt_s);

            //terminate simualation loop if all particles are terminated or termination of the integrator was requested
            //from somewhere (e.g. signal from outside)
            if (trajectoryIntegrator !=nullptr) {
                trajectoryIntegrator->runSingleStep(dt_s);
                if (AppUtils::SignalHandler::isTerminationSignaled()){
                    break;
                }
            }
            if (ionsInactive>=nAllParticles){
                break;
            }
        }
        resultFilewriter.writeReactionStatistics(rsSim);
        if (trajectoryIntegrator) {
            trajectoryIntegrator->finalizeSimulation();
        }
        resultFilewriter.closeFile();

        stopWatch.stop();
        logger->info("----------------------");
        logger->info("Reaction Events:");
        rsSim.logReactionStatistics(logger);
        logger->info("----------------------");
        logger->info("total reaction events: {} ill events: {}", rsSim.totalReactionEvents(), rsSim.illEvents());
        logger->info("ill fraction: {}", rsSim.illEvents()/(double) rsSim.totalReactionEvents());

        logger->info("CPU time: {} s", stopWatch.elapsedSecondsCPU());
        logger->info("Finished in {} seconds (wall clock time)", stopWatch.elapsedSecondsWall());

        // ======================================================================================

        return 0;
    }
    catch(const std::invalid_argument& ia){
        std::cout << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
}
