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
 BT-RS-DMSSim.cpp

 Idealized plane electrode type differential ion mobility spectrometry (DMS) transport and chemistry simulation,
 including space chage and gas collision effects

 ****************************/

#include "Core_utils.hpp"
#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_ConcentrationFileWriter.hpp"
#include "PSim_util.hpp"
#include "PSim_constants.hpp"
#include "PSim_verletIntegrator.hpp"
#include "PSim_trajectoryExplorerJSONwriter.hpp"
#include "PSim_scalar_writer.hpp"
#include "CollisionModel_StatisticalDiffusion.hpp"
#include "CollisionModel_SpatialFieldFunctions.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "json.h"
#include <iostream>
#include <cmath>

const std::string key_ChemicalIndex = "keyChemicalIndex";
enum CVMode {STATIC_CV,AUTO_CV};
enum FlowMode {UNIFORM_FLOW,PARABOLIC_FLOW};
enum BackgroundTemperatureMode {ISOTHERM,LINEAR_GRADIENT};
enum CollisionType {SDS,NO_COLLISION};

int main(int argc, const char * argv[]) {

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

        std::vector<unsigned int> nParticles = simConf.unsignedIntVectorParameter("n_particles");
        unsigned int nSteps = simConf.unsignedIntParameter("sim_time_steps");
        unsigned int nStepsPerOscillation = simConf.unsignedIntParameter("sim_time_steps_per_sv_oscillation");
        int concentrationWriteInterval = simConf.intParameter("concentrations_write_interval");
        int trajectoryWriteInterval = simConf.intParameter("trajectory_write_interval");
        double spaceChargeFactor = simConf.doubleParameter("space_charge_factor");


        //geometric parameters:
        double startWidthX_m = simConf.doubleParameter("start_width_x_mm")/1000.0;
        double startWidthY_m = simConf.doubleParameter("start_width_y_mm")/1000.0;
        double startWidthZ_m = simConf.doubleParameter("start_width_z_mm")/1000.0;
        double electrodeDistance_m = simConf.doubleParameter("electrode_distance_mm")/1000.0;
        double electrodeLength_m = simConf.doubleParameter("electrode_length_mm")/1000.0;
        double electrodeHalfDistance_m = electrodeDistance_m/2.0;
        double electrodeHalfDistanceSquared_m = electrodeHalfDistance_m*electrodeHalfDistance_m;


        //background gas parameters:
        std::string collisionTypeStr = simConf.stringParameter("collision_model");
        CollisionType collisionType;
        if (collisionTypeStr=="SDS") {
            collisionType = SDS;
        }
        else if (collisionTypeStr=="none") {
            collisionType = NO_COLLISION;
        }
        else {
            throw std::invalid_argument("wrong configuration value: collision_model_type");
        }

        std::string flowModeStr = simConf.stringParameter("flow_mode");
        FlowMode flowMode;
        if (flowModeStr=="uniform") {
            flowMode = UNIFORM_FLOW;
        }
        else if (flowModeStr=="parabolic") {
            flowMode = PARABOLIC_FLOW;
        }
        else {
            throw std::invalid_argument("wrong configuration value: flow_mode");
        }

        //Define background temperature function for chemical reaction and collision model
        //BackgroundTemperatureMode backgroundTempMode;
        std::function<double(const Core::Vector&)> backgroundTemperatureFct;
        std::string backgroundTempStr = simConf.stringParameter("background_temperature_mode");
        if (backgroundTempStr=="isotherm") {
            //backgroundTempMode = ISOTHERM;
            double backgroundTemperature_K = simConf.doubleParameter("background_temperature_K");
            backgroundTemperatureFct = CollisionModel::getConstantDoubleFunction(backgroundTemperature_K);
        }
        else if (backgroundTempStr=="linear_gradient") {
            //backgroundTempMode = LINEAR_GRADIENT;
            double backgroundTemp_start = simConf.doubleParameter("background_temperature_start_K");
            double backgroundTemp_stop = simConf.doubleParameter("background_temperature_stop_K");
            double backgroundTemp_diff = backgroundTemp_stop-backgroundTemp_start;

            backgroundTemperatureFct =
                    [backgroundTemp_start,
                            backgroundTemp_stop,
                            backgroundTemp_diff,
                            electrodeLength_m](const Core::Vector& particleLocation) -> double {

                        if (particleLocation.x()>electrodeLength_m) {
                            return backgroundTemp_stop;
                        }
                        else {
                            return (backgroundTemp_diff/electrodeLength_m*particleLocation.x())+backgroundTemp_start;
                        }
                    };
        }
        else {
            throw std::invalid_argument("wrong configuration value: background_temperature_mode");
        }

        double backgroundPressure_Pa = simConf.doubleParameter("background_pressure_Pa");
        double gasVelocityX = simConf.doubleParameter("collision_gas_velocity_x_ms-1");
        double collisionGasMass_Amu = simConf.doubleParameter("collision_gas_mass_amu");
        double collisionGasDiameter_nm = simConf.doubleParameter("collision_gas_diameter_nm");


        //field parameters:
        std::string cvModeStr = simConf.stringParameter("cv_mode");
        CVMode cvMode;
        double meanZPos = 0.0; //variable used for automatic CV correction
        double cvRelaxationParameter = 0.0;
        if (cvModeStr=="static") {
            cvMode = STATIC_CV;
        }
        else if (cvModeStr=="auto") {
            cvMode = AUTO_CV;
            cvRelaxationParameter = simConf.doubleParameter("cv_relaxation_parameter");
        }
        else {
            throw std::invalid_argument("wrong configuration value: cv_mode");
        }

        double fieldSV_VPerM = simConf.doubleParameter("sv_Vmm-1")*1000.0;
        double fieldCV_VPerM = simConf.doubleParameter("cv_Vmm-1")*1000.0;
        double fieldFrequency = simConf.doubleParameter("sv_frequency_s-1");
        double fieldWavePeriod = 1.0/fieldFrequency;
        double field_h = 2.0;
        double field_F = 2.0;
        double field_W = (1/fieldWavePeriod)*2*M_PI;
        double fieldMagnitude = 0; //actual field magnitude

        double dt_s = fieldWavePeriod/nStepsPerOscillation;
        // ======================================================================================


        //read and prepare chemical configuration ===============================================
        RS::ConfigFileParser parser = RS::ConfigFileParser();
        std::string rsConfFileName = simConf.pathRelativeToConfFile(simConf.stringParameter("reaction_configuration"));
        RS::Simulation rsSim = RS::Simulation(parser.parseFile(rsConfFileName));
        RS::SimulationConfiguration* rsSimConf = rsSim.simulationConfiguration();
        //prepare a map for retrieval of the substance index:
        std::map<RS::Substance*, int> substanceIndices;
        std::vector<RS::Substance*> discreteSubstances = rsSimConf->getAllDiscreteSubstances();
        std::vector<double> ionMobility; // = simConf.doubleVectorParameter("ion_mobility");
        for (std::size_t i = 0; i<discreteSubstances.size(); ++i) {
            substanceIndices.insert(std::pair<RS::Substance*, int>(discreteSubstances[i], i));
            ionMobility.push_back(discreteSubstances[i]->mobility());
        }


        // prepare file writer  =================================================================
        RS::ConcentrationFileWriter resultFilewriter(projectName+"_conc.csv");

        auto jsonWriter = std::make_unique<ParticleSimulation::TrajectoryExplorerJSONwriter>(
                projectName+"_trajectories.json");
        jsonWriter->setScales(1000, 1);

        unsigned int ionsInactive = 0;
        unsigned int nAllParticles = 0;
        for (const auto ni: nParticles) {
            nAllParticles += ni;
        }

        std::unique_ptr<ParticleSimulation::Scalar_writer> cvFieldWriter;
        if (cvMode==AUTO_CV) {
            cvFieldWriter = std::make_unique<ParticleSimulation::Scalar_writer>(projectName+"_cv.csv");
        }

        std::unique_ptr<ParticleSimulation::Scalar_writer> voltageWriter;
        voltageWriter = std::make_unique<ParticleSimulation::Scalar_writer>(projectName+"_voltages.csv");


        // init simulation  =====================================================================

        // create and add simulation particles:
        unsigned int nParticlesTotal = 0;
        std::vector<uniqueReactivePartPtr> particles;
        std::vector<BTree::Particle*> particlesPtrs;
        std::vector<std::vector<double>> trajectoryAdditionalParams;

        Core::Vector initCorner(0, -startWidthY_m/2.0, -startWidthZ_m/2.0);
        Core::Vector initBoxSize(startWidthX_m, startWidthY_m, startWidthZ_m);

        for (std::size_t i = 0; i<nParticles.size(); i++) {
            RS::Substance* subst = rsSimConf->substance(i);
            std::vector<Core::Vector> initialPositions =
                    ParticleSimulation::util::getRandomPositionsInBox(nParticles[i], initCorner, initBoxSize);
            for (unsigned int k = 0; k<nParticles[i]; ++k) {
                uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(subst);

                particle->setLocation(initialPositions[k]);

                particlesPtrs.push_back(particle.get());
                rsSim.addParticle(particle.get(), nParticlesTotal);
                particles.push_back(std::move(particle));
                trajectoryAdditionalParams.push_back(std::vector<double>(1));
                nParticlesTotal++;
            }
        }

        RS::ReactionConditions reactionConditions = RS::ReactionConditions();
        reactionConditions.temperature = 0.0;//backgroundTemperature_K;
        reactionConditions.pressure = backgroundPressure_Pa;
        reactionConditions.electricField = 0.0;

        resultFilewriter.initFile(rsSimConf);
        // ======================================================================================


        // define trajectory integration parameters / functions =================================

        auto accelerationFct =
                [&fieldSV_VPerM, &fieldCV_VPerM, field_F, field_W, field_h, &fieldMagnitude, spaceChargeFactor]
                        (BTree::Particle* particle, int /*particleIndex*/, BTree::Tree& tree, double time, int /*timestep*/) {

                    double particleCharge = particle->getCharge();
                    double voltageSVgp = fieldSV_VPerM*0.6667; // V/m (1V/m peak to peak is 0.6667V/m ground to peak)
                    double voltageSVt = fieldCV_VPerM+(field_F*sin(field_W*time)
                            +sin(field_h*field_W*time-0.5*M_PI))
                            *voltageSVgp/(field_F+1);

                    fieldMagnitude = voltageSVt; //// / electrodeDistance_m;

                    Core::Vector fieldForce(0, 0, fieldMagnitude*particleCharge);

                    if (Core::isDoubleEqual(spaceChargeFactor, 0.0)) {
                        return (fieldForce/particle->getMass());
                    }
                    else {
                        Core::Vector spaceChargeForce =
                                tree.computeEFieldFromTree(*particle)*(particleCharge*spaceChargeFactor);
                        return ((fieldForce+spaceChargeForce)/particle->getMass());
                    }
                };

        ParticleSimulation::partAttribTransformFctType additionalParameterTransformFct =
                [=](BTree::Particle* particle) -> std::vector<double> {
                    std::vector<double> result = {particle->getFloatAttribute(key_ChemicalIndex)};
                    return result;
                };

        auto timestepWriteFct =
                [&jsonWriter, &voltageWriter, &additionalParameterTransformFct, trajectoryWriteInterval,
                        &rsSim, &resultFilewriter, concentrationWriteInterval, &fieldMagnitude, &logger]
                        (std::vector<BTree::Particle*>& particles, BTree::Tree& tree, double time, int timestep,
                         bool lastTimestep) {

                    if (timestep%concentrationWriteInterval==0) {
                        resultFilewriter.writeTimestep(rsSim);
                        voltageWriter->writeTimestep(fieldMagnitude, time);
                    }
                    if (lastTimestep) {
                        jsonWriter->writeTimestep(
                                particles, additionalParameterTransformFct, time, true);

                        jsonWriter->writeSplatTimes(particles);
                        jsonWriter->writeIonMasses(particles);
                        logger->info("finished ts:{} time:{:.2e}", timestep, time);
                    }

                    else if (timestep%trajectoryWriteInterval==0) {
                        Core::Vector cof = tree.getRoot()->getCenterOfCharge();
                        logger->info("ts:{}  time:{:.2e} average cloud position:({:.2e},{:.2e},{:.2e})", timestep, time,
                                cof.x(), cof.y(), cof.z());
                        rsSim.logConcentrations(logger);
                        jsonWriter->writeTimestep(
                                particles, additionalParameterTransformFct, time, false);
                    }
                };

        auto otherActionsFct = [electrodeHalfDistance_m, electrodeLength_m, &ionsInactive](
                Core::Vector& newPartPos, BTree::Particle* particle,
                int /*particleIndex*/, BTree::Tree& /*tree*/, double time, int /*timestep*/) {

            if (std::fabs(newPartPos.z())>=electrodeHalfDistance_m) {
                particle->setActive(false);
                particle->setSplatTime(time);
                ionsInactive++;
            }
            else if (newPartPos.x()>=electrodeLength_m) {
                particle->setActive(false);
                ionsInactive++;
            }
        };


        //define / gas interaction /  collision model:
        std::unique_ptr<CollisionModel::AbstractCollisionModel> collisionModelPtr;
        if (collisionType==SDS) {
            // prepare static pressure and temperature functions
            auto staticPressureFct = CollisionModel::getConstantDoubleFunction(backgroundPressure_Pa);

            std::function<Core::Vector(const Core::Vector&)> velocityFct;

            if (flowMode==UNIFORM_FLOW) {
                velocityFct =
                        [gasVelocityX](const Core::Vector& /*pos*/) {
                            return Core::Vector(gasVelocityX, 0.0, 0.0);
                        };
            }
            else if (flowMode==PARABOLIC_FLOW) {
                velocityFct =
                        [gasVelocityX, electrodeHalfDistanceSquared_m](const Core::Vector& pos) {
                            //parabolic profile is vX = 2 * Vavg * (1 - r^2 / R^2) with the radius / electrode distance R
                            double xVelo = gasVelocityX*2.0*(1-pos.z()*pos.z()/electrodeHalfDistanceSquared_m);
                            return Core::Vector(xVelo, 0.0, 0.0);
                        };
            }

            std::unique_ptr<CollisionModel::StatisticalDiffusionModel> collisionModel =
                    std::make_unique<CollisionModel::StatisticalDiffusionModel>(
                            staticPressureFct,
                            backgroundTemperatureFct,
                            velocityFct,
                            collisionGasMass_Amu,
                            collisionGasDiameter_nm*1e-9);

            for (const auto& particle: particlesPtrs) {
                particle->setDiameter(
                        CollisionModel::util::estimateCollisionDiameterFromMass(
                                particle->getMass()/Core::AMU_TO_KG
                        )*1e-9);
                collisionModel->setSTPParameters(*particle);
            }
            collisionModelPtr = std::move(collisionModel);
        }
        else if (collisionType==NO_COLLISION) {
            collisionModelPtr = nullptr;
        }

        //init trajectory simulation object:
        ParticleSimulation::VerletIntegrator verletIntegrator(
                particlesPtrs,
                accelerationFct, timestepWriteFct, otherActionsFct, ParticleSimulation::noFunction,
                collisionModelPtr.get());
        // ======================================================================================


        // simulate   ===========================================================================
        AppUtils::SignalHandler::setReceiver(verletIntegrator);
        AppUtils::Stopwatch stopWatch;
        stopWatch.start();

        for (unsigned int step = 0; step<nSteps; step++) {
            for (unsigned int i = 0; i<nParticlesTotal; i++) {
                reactionConditions.electricField = fieldMagnitude;
                reactionConditions.temperature = backgroundTemperatureFct(particles[i]->getLocation());

                bool reacted = rsSim.react(i, reactionConditions, dt_s);
                int substIndex = substanceIndices.at(particles[i]->getSpecies());
                particles[i]->setFloatAttribute(key_ChemicalIndex, substIndex);

                if (reacted && collisionModelPtr != nullptr) {
                    //we had an reaction event: update the collision model parameters for the particle which are not
                    //based on location (mostly STP parameters in SDS)
                    collisionModelPtr->initializeModelParameters(*particles[i]);
                }
            }
            rsSim.advanceTimestep(dt_s);
            verletIntegrator.runSingleStep(dt_s);

            //autocorrect compensation voltage, to minimize z drift (once for every single SV oscillation):
            if (cvMode==AUTO_CV && step%nStepsPerOscillation==0) {
                //calculate current mean z-position:
                double buf = 0.0;
                for (unsigned int i = 0; i<nParticlesTotal; i++) {
                    buf += particles[i]->getLocation().z();
                }
                double currentMeanZPos = buf/nParticlesTotal;

                //update cv value:
                double diffMeanZPos = meanZPos-currentMeanZPos;
                fieldCV_VPerM = fieldCV_VPerM+diffMeanZPos*cvRelaxationParameter;
                cvFieldWriter->writeTimestep(std::vector<double>{fieldCV_VPerM, currentMeanZPos},
                        rsSim.simulationTime());
                meanZPos = currentMeanZPos;
                logger->info("CV corrected ts:{} time:{:.2e} new CV:", step, rsSim.simulationTime(), fieldCV_VPerM);
            }

            //terminate simualation loops if all particles are terminated or termination of the integrator was requested
            //from somewhere (e.g. signal from outside)
            if (ionsInactive>=nAllParticles ||
                    verletIntegrator.runState()==ParticleSimulation::AbstractTimeIntegrator::IN_TERMINATION)
            {
                break;
            }
        }
        verletIntegrator.finalizeSimulation();
        resultFilewriter.closeFile();
        stopWatch.stop();

        logger->info("total reaction events: {} ill events: {}", rsSim.totalReactionEvents(), rsSim.illEvents());
        logger->info("ill fraction: {}", rsSim.illEvents()/(double) rsSim.totalReactionEvents());
        logger->info("CPU time: {} s", stopWatch.elapsedSecondsCPU());
        logger->info("Finished in {} seconds (wall clock time)", stopWatch.elapsedSecondsWall());
        // ======================================================================================

        return EXIT_SUCCESS;
    }
    catch(const std::invalid_argument& ia){
        std::cout << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
}
