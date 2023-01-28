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
 DMSSim.cpp

 Idealized plane electrode type differential ion mobility spectrometry (DMS) transport and chemistry simulation,
 including space chage and gas collision effects

 ****************************/

#include "Core_utils.hpp"
#include "Core_randomGenerators.hpp"
#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_ConcentrationFileWriter.hpp"
#include "PSim_util.hpp"
#include "PSim_constants.hpp"
#include "Integration_parallelVerletIntegrator.hpp"
#include "FileIO_trajectoryHDF5Writer.hpp"
#include "FileIO_scalar_writer.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include "CollisionModel_StatisticalDiffusion.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "CollisionModel_MDInteractions.hpp"
#include "CollisionModel_SpatialFieldFunctions.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "appUtils_commandlineParser.hpp"
#include "dmsSim_dmsFields.hpp"
#include "PSim_simionPotentialArray.hpp"
#include "PSim_particleStartSplatTracker.hpp"
#include <iostream>
#include <cmath>

const std::string key_ChemicalIndex = "keyChemicalIndex";
enum CollisionType {SDS, HS, MD, NO_COLLISION};

/*std::vector<double> &paVoltageFct (std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> potentialArrays,
                                   std::string WaveformFilename,
                                std::vector<double> phaseShift, double wavePeriod, double waveAmplitude, double time) {
    auto WaveForm = ParticleSimulation::SampledWaveform(WaveformFilename);
    std::vector<double> currentVoltages;
    for(size_t i=0; i<potentialArrays.size(); i++) {
        double period = std::fmod(time, wavePeriod) / wavePeriod;
        double shiftedPeriod = std::fmod(period + phaseShift[i], 1.0);
        currentVoltages.push_back(WaveForm.getInterpolatedValue(shiftedPeriod) * waveAmplitude);
    }
    return currentVoltages;
} */

int main(int argc, const char * argv[]) {

    try {
        // open configuration, parse configuration file =========================================
        AppUtils::CommandlineParser cmdLineParser(argc, argv, "BT-RS-DMSSim", "DMS Simulation with trajectories and chemistry", true);
        std::string projectName = cmdLineParser.resultName();
        AppUtils::logger_ptr logger = cmdLineParser.logger();

        std::string confFileName = cmdLineParser.confFileName();
        AppUtils::simConf_ptr simConf = cmdLineParser.simulationConfiguration();

        // optionally setting random generator seed manually (for debugging / reproduction purposes):
        if (simConf->isParameter("random_seed")) {
            unsigned int randomSeed = simConf->unsignedIntParameter("random_seed");
            Core::globalRandomGeneratorPool->setSeedForElements(randomSeed);
        }


        std::vector<unsigned int> nParticles = simConf->unsignedIntVectorParameter("n_particles");
        unsigned int nSteps = simConf->unsignedIntParameter("sim_time_steps");
        double dt_s = simConf->doubleParameter("dt_s");
        int concentrationWriteInterval = simConf->intParameter("concentrations_write_interval");
        int trajectoryWriteInterval = simConf->intParameter("trajectory_write_interval");
        double spaceChargeFactor = simConf->doubleParameter("space_charge_factor");
        double waveAmplitude = simConf->doubleParameter("wave_amplitude_V");
        double waveFrequency = simConf->doubleParameter("wave_frequency_hz");
        //double waveSpeed = simConf->doubleParameter("wave_Speed_m/s");

        //geometric parameters:
        double startWidthX_m = simConf->doubleParameter("start_width_x_mm")/1000.0;
        double startWidthY_m = simConf->doubleParameter("start_width_y_mm")/1000.0;
        double startWidthZ_m = simConf->doubleParameter("start_width_z_mm")/1000.0;


        //background gas parameters:
        std::string collisionTypeStr = simConf->stringParameter("collision_model");
        CollisionType collisionType;
        if (collisionTypeStr=="SDS") {
            collisionType = SDS;
        }
        else if (collisionTypeStr=="HS"){
            collisionType = HS;
        }
        else if (collisionTypeStr=="MD"){
            collisionType = MD;
        }
        else if (collisionTypeStr=="none") {
            collisionType = NO_COLLISION;
        }
        else {
            throw std::invalid_argument("wrong configuration value: collision_model_type");
        }

        double backgroundPressure_Pa = simConf->doubleParameter("background_pressure_Pa");
        double collisionGasMass_Amu = simConf->doubleParameter("collision_gas_mass_amu");
        double collisionGasDiameter_nm = simConf->doubleParameter("collision_gas_diameter_nm");
        double backgroundTemperature_K = simConf->doubleParameter("background_temperature_K");

        // ======================================================================================

        //read potential array configuration of the trap =================================================
        std::filesystem::path confBasePath = simConf->confBasePath();

        double paSpatialScale = simConf->doubleParameter("potential_array_scale");
        std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> potentialArrays;
        std::vector<std::string> potentialArraysNames = simConf->stringVectorParameter("potential_arrays");
        for (const auto& paName: potentialArraysNames) {
            std::filesystem::path paPath = confBasePath/paName;
            std::unique_ptr<ParticleSimulation::SimionPotentialArray> pa_pt =
                    std::make_unique<ParticleSimulation::SimionPotentialArray>(paPath, paSpatialScale);
            potentialArrays.push_back(std::move(pa_pt));
        }

        double potentialScale = 1.0/10000.0;

        std::array<double, 6> paBounds = potentialArrays[0]->getBounds();

        double wavePeriod = 1.0/waveFrequency;

        std::vector<double> phaseShift = simConf->doubleVectorParameter("phase_shift");
        std::string WaveformFilename = simConf->pathRelativeToConfFile(simConf->stringParameter("waveform"));
        auto WaveForm = ParticleSimulation::SampledWaveform(WaveformFilename);

        // ======================================================================================

        //read and prepare chemical configuration ===============================================
        RS::ConfigFileParser parser = RS::ConfigFileParser();
        std::string rsConfFileName = simConf->pathRelativeToConfFile(simConf->stringParameter("reaction_configuration"));
        RS::Simulation rsSim = RS::Simulation(parser.parseFile(rsConfFileName));
        RS::SimulationConfiguration* rsSimConf = rsSim.simulationConfiguration();

        //prepare a map for retrieval of the substance index:
        std::map<RS::Substance*, int> substanceIndices;
        std::vector<RS::Substance*> discreteSubstances = rsSimConf->getAllDiscreteSubstances();
        std::vector<double> ionMobility; // = simConf->doubleVectorParameter("ion_mobility");
        for (std::size_t i = 0; i<discreteSubstances.size(); ++i) {
            substanceIndices.insert(std::pair<RS::Substance*, int>(discreteSubstances[i], i));
            ionMobility.push_back(discreteSubstances[i]->lowFieldMobility());
        }

        //read molecular structure file
        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection;
        if(collisionType==MD){
            std::string mdCollisionConfFile = simConf->pathRelativeToConfFile(simConf->stringParameter("md_configuration"));
            FileIO::MolecularStructureReader mdConfReader = FileIO::MolecularStructureReader();
            molecularStructureCollection = mdConfReader.readMolecularStructure(mdCollisionConfFile);
        }


        // prepare file writer  =================================================================
        RS::ConcentrationFileWriter resultFilewriter(projectName+"_conc.csv");

        std::vector<std::string> auxParamNames = {"chemical id"};
        auto additionalParamTFct = [](Core::Particle* particle) -> std::vector<int> {
            std::vector<int> result = {
                    particle->getIntegerAttribute(key_ChemicalIndex)
            };
            return result;
        };

        std::string hdf5Filename = projectName+"_trajectories.hd5";
        FileIO::TrajectoryHDF5Writer trajectoryWriter(hdf5Filename);
        trajectoryWriter.setParticleAttributes(auxParamNames, additionalParamTFct);

        unsigned int ionsInactive = 0;
        unsigned int nAllParticles = 0;
        for (const auto ni: nParticles) {
            nAllParticles += ni;
        }

        std::unique_ptr<FileIO::Scalar_writer> voltageWriter;
        voltageWriter = std::make_unique<FileIO::Scalar_writer>(projectName+"_voltages.csv");


        // init simulation  =====================================================================

        // create and add simulation particles:
        unsigned int nParticlesTotal = 0;
        std::vector<uniqueReactivePartPtr> particles;
        std::vector<Core::Particle*> particlesPtrs;
        std::vector<std::vector<double>> trajectoryAdditionalParams;

        Core::Vector initCorner(2.5e-04, -startWidthY_m/2.0, -startWidthZ_m/2.0);
        Core::Vector initBoxSize(startWidthX_m, startWidthY_m, startWidthZ_m);

        for (std::size_t i = 0; i<nParticles.size(); i++) {
            RS::Substance* subst = rsSimConf->substance(i);
            int substIndex = substanceIndices.at(subst);
            std::vector<Core::Vector> initialPositions =
                    ParticleSimulation::util::getRandomPositionsInBox(nParticles[i], initCorner, initBoxSize);
            for (unsigned int k = 0; k<nParticles[i]; ++k) {
                uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(subst);

                particle->setLocation(initialPositions[k]);
                particle->setIntegerAttribute(key_ChemicalIndex, substIndex);
                particle->setFloatAttribute("effectiveField", 0.0);
                particle->setIndex(nParticlesTotal);

                particlesPtrs.push_back(particle.get());
                rsSim.addParticle(particle.get(), nParticlesTotal);
                particles.push_back(std::move(particle));
                trajectoryAdditionalParams.push_back(std::vector<double>(1));
                nParticlesTotal++;
            }
        }


        resultFilewriter.initFile(rsSimConf);
        // ======================================================================================


        // define trajectory integration parameters / functions =================================
        std::vector<double> totalFieldNow(potentialArrays.size(), 0.0);

        auto paVoltageFct = [&potentialArrays, &WaveForm, phaseShift,
                             wavePeriod, waveAmplitude, &totalFieldNow](double time)->std::vector<double> {
            //std::vector<double> currentVoltages;
            for(size_t i=0; i<potentialArrays.size(); i++) {
                double period = std::fmod(time, wavePeriod) / wavePeriod;
                double shiftedPeriod = std::fmod(period + phaseShift[i], 1.0);
                totalFieldNow[i]=(WaveForm.getInterpolatedValue(shiftedPeriod) * waveAmplitude);
            }
        };

        auto accelerationFct =
                [&potentialArrays, &totalFieldNow, potentialScale]
                        (Core::Particle* particle, int /*particleIndex*/, SpaceCharge::FieldCalculator &fieldCalculator,
                                double /*time*/, int /*timestep*/){
                    Core::Vector fEfield(0, 0, 0);
                    Core::Vector pos = particle->getLocation();
                    double particleCharge = particle->getCharge();

                    for(size_t i=0; i<potentialArrays.size(); i++) {
                        Core::Vector paField = potentialArrays[i]->getField(pos.x(), pos.y(), pos.z());
                        Core::Vector paEffectiveField = paField * totalFieldNow[i] * potentialScale;

                        fEfield = fEfield+paEffectiveField;
                    }
                    particle->setFloatAttribute("effectiveField", fEfield.magnitude());
                    return (fEfield*particleCharge/particle->getMass());
                };


        auto timestepWriteFct =
                [&trajectoryWriter, &voltageWriter, trajectoryWriteInterval, &rsSim, &resultFilewriter, concentrationWriteInterval,
                 &totalFieldNow, &logger]
                        (std::vector<Core::Particle*>& particles, double time, int timestep,
                         bool lastTimestep) {

                    if (timestep%concentrationWriteInterval==0) {
                        resultFilewriter.writeTimestep(rsSim);
                        voltageWriter->writeTimestep(totalFieldNow, time);
                    }
                    if (lastTimestep) {
                        trajectoryWriter.writeTimestep(particles, time);
                        trajectoryWriter.writeSplatTimes(particles);
                        trajectoryWriter.finalizeTrajectory();
                        logger->info("finished ts:{} time:{:.2e}", timestep, time);
                    }
                    else if (timestep%trajectoryWriteInterval==0) {
                        logger->info("ts:{}  time:{:.2e}",
                                timestep, time);
                        rsSim.logConcentrations(logger);
                        trajectoryWriter.writeTimestep(particles, time);
                    }
                };

        ParticleSimulation::ParticleStartSplatTracker startSplatTracker;
        auto particleStartMonitoringFct = [&startSplatTracker](Core::Particle* particle, double time) {
            startSplatTracker.particleStart(particle, time);
        };

        auto otherActionsFct = [&ionsInactive, &startSplatTracker](
                Core::Vector& /*newPartPos*/, Core::Particle* particle,
                int /*particleIndex*/,  double time, int /*timestep*/) {
            Core::Vector pos = particle->getLocation();
            if (pos.x() == 0) {
                particle->setActive(false);
                particle->setSplatTime(time);
                startSplatTracker.particleSplat(particle, time);
                ionsInactive++;
            }
        };

        //define / gas interaction /  collision model:
        std::unique_ptr<CollisionModel::AbstractCollisionModel> collisionModelPtr;
        if (collisionType==SDS || collisionType==HS || collisionType==MD) {
            // prepare static pressure and temperature functions

            if (collisionType==SDS){
                std::unique_ptr<CollisionModel::StatisticalDiffusionModel> collisionModel =
                        std::make_unique<CollisionModel::StatisticalDiffusionModel>(
                                backgroundPressure_Pa,
                                backgroundTemperature_K,
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
            else if (collisionType==HS){
                std::unique_ptr<CollisionModel::HardSphereModel> collisionModel =
                        std::make_unique<CollisionModel::HardSphereModel>(
                                backgroundPressure_Pa,
                                backgroundTemperature_K,
                                collisionGasMass_Amu,
                                collisionGasDiameter_nm*1e-9,
                                nullptr);
                collisionModelPtr = std::move(collisionModel);
            }
            else if (collisionType==MD){

                // collect additional config file parameters for MD model:
                double collisionGasPolarizability_m3 = simConf->doubleParameter("collision_gas_polarizability_m3");
                std::string collisionGasIdentifier = simConf->stringParameter("collision_gas_identifier");
                std::vector<std::string> particleIdentifiers = simConf->stringVectorParameter("particle_identifier");
                double subIntegratorIntegrationTime_s = simConf->doubleParameter("sub_integrator_integration_time_s");
                double subIntegratorStepSize_s = simConf->doubleParameter("sub_integrator_step_size_s");
                double collisionRadiusScaling = simConf->doubleParameter("collision_radius_scaling");
                double angleThetaScaling = simConf->doubleParameter("angle_theta_scaling");
                double spawnRadius_m = simConf->doubleParameter("spawn_radius_m");

                //construct MD model:
                std::unique_ptr<CollisionModel::MDInteractionsModel> collisionModel =
                        std::make_unique<CollisionModel::MDInteractionsModel>(
                            backgroundPressure_Pa,
                            backgroundTemperature_K,
                            collisionGasMass_Amu,
                            collisionGasDiameter_nm*1e-9,
                            collisionGasPolarizability_m3,
                            collisionGasIdentifier,
                            subIntegratorIntegrationTime_s,
                            subIntegratorStepSize_s,
                            collisionRadiusScaling,
                            angleThetaScaling,
                            spawnRadius_m,
                            molecularStructureCollection);

                // Set trajectory writing options:
                bool saveTrajectory = simConf->boolParameter("save_trajectory");

                if (saveTrajectory){
                    int saveTrajectoryStartTimeStep = simConf->intParameter("trajectory_start_time_step");
                    double trajectoryDistance_m = simConf->doubleParameter("trajectory_distance_m");
                    collisionModel->setTrajectoryWriter(projectName+"_md_trajectories.txt",
                            trajectoryDistance_m, saveTrajectoryStartTimeStep);
                }


                // Init particles with MD parameters:
                unsigned int particleIndex = 0;
                for (std::size_t i = 0; i<nParticles.size(); i++) {
                    for (unsigned int k = 0; k<nParticles[i]; k++) {
                        Core::Particle* particle = particlesPtrs[particleIndex];
                        particle->setMolecularStructure(molecularStructureCollection.at(particleIdentifiers[i]));
                        particle->setDiameter(particle->getMolecularStructure()->getDiameter());
                        particleIndex++;
                    }
                }

                collisionModelPtr = std::move(collisionModel);
            }
        }
        else if (collisionType==NO_COLLISION) {
            collisionModelPtr = nullptr;
        }

        //define reaction simulation functions:
        auto particlesHasReactedFct = [&collisionModelPtr, &substanceIndices](RS::ReactiveParticle* particle){
            //we had an reaction event: Count it (access to total counted value has to be synchronized)
            if (collisionModelPtr != nullptr) {
                collisionModelPtr->initializeModelParticleParameters(*particle);
            }
            int substIndex = substanceIndices.at(particle->getSpecies());
            particle->setIntegerAttribute(key_ChemicalIndex, substIndex);
        };

        auto reactionConditionsFct = [&totalFieldNow, backgroundTemperature_K, backgroundPressure_Pa]
                (RS::ReactiveParticle* particle, double /*time*/)->RS::ReactionConditions{
            RS::ReactionConditions reactionConditions = RS::ReactionConditions();

            reactionConditions.temperature = backgroundTemperature_K;
            //reactionConditions.electricField = totalFieldNow_VPerM;
            reactionConditions.pressure = backgroundPressure_Pa;
            return reactionConditions;
        };

        //init trajectory simulation object:
        Integration::ParallelVerletIntegrator verletIntegrator(
                particlesPtrs,
                accelerationFct, timestepWriteFct, otherActionsFct, ParticleSimulation::noFunction,
                collisionModelPtr.get());
        // ======================================================================================


        // simulate   ===========================================================================
        AppUtils::SignalHandler::setReceiver(verletIntegrator);
        AppUtils::Stopwatch stopWatch;
        stopWatch.start();

        for (unsigned int step = 0; step<nSteps; step++) {
            rsSim.performTimestep(reactionConditionsFct, dt_s, particlesHasReactedFct);
            rsSim.advanceTimestep(dt_s);
            paVoltageFct(rsSim.simulationTime());
            verletIntegrator.runSingleStep(dt_s);


            //terminate simulation loops if all particles are terminated or termination of the integrator was requested
            //from somewhere (e.g. signal from outside)
            if (ionsInactive>=nAllParticles ||
                    verletIntegrator.runState()==Integration::AbstractTimeIntegrator::IN_TERMINATION)
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
    catch(AppUtils::TerminatedWhileCommandlineParsing& terminatedMessage){
        return terminatedMessage.returnCode();
    }
    catch(const std::invalid_argument& ia){
        std::cout << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
}