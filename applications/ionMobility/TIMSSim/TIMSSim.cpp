/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2024 - Physical and Theoretical Chemistry /
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
 TIMSSim.cpp

 Trapped ion mobility spectrometry (TIMS) transport and chemistry simulation,
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
#include "FileIO_CSVReader.hpp"
#include "CollisionModel_StatisticalDiffusion.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "CollisionModel_MDInteractions.hpp"
#include "CollisionModel_MDForceField_LJ12_6.hpp"
#include "CollisionModel_SpatialFieldFunctions.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "appUtils_commandlineParser.hpp"
#include "dmsSim_dmsFields.hpp"
#include "PSim_simionPotentialArray.hpp"
#include "PSim_particleStartSplatTracker.hpp"
#include "CollisionModel_MultiCollisionModel.hpp"
#include <iostream>
#include <cmath>

const std::string key_ChemicalIndex = "keyChemicalIndex";
enum FlowMode {PARABOLIC_FLOW, UNIFORM_FLOW};
enum CollisionType {SDS, HS, MD, NO_COLLISION};

int main(int argc, const char * argv[]) {

    try {
        // open configuration, parse configuration file =========================================
        AppUtils::CommandlineParser cmdLineParser(argc, argv, "TIMSSim", "TIMS Simulation with trajectories and chemistry", true);
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

        double f_rf = simConf->doubleParameter("confining_RF_frequency_Hz");
        [[maybe_unused]] double omega = f_rf*2.0*M_PI;
        double V_rf = simConf->doubleParameter("confining_RF_amplitude_V");

        //geometric parameters:
        std::vector<double> startPosition_mm = simConf->doubleVectorParameter("start_box_position_mm");
        std::vector<double> startWidth_mm = simConf->doubleVectorParameter("start_box_dimensions_mm");

        double startWidthX_m = startWidth_mm[0]/1000.0;
        double startWidthY_m = startWidth_mm[1]/1000.0;
        double startWidthZ_m = startWidth_mm[2]/1000.0;


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

        FlowMode flowMode;
        if (simConf->isParameter("flow_mode")) {
            std::string flowModeStr = simConf->stringParameter("flow_mode");
            if (flowModeStr == "uniform") {
                flowMode = UNIFORM_FLOW;
            } else if (flowModeStr == "parabolic") {
                flowMode = PARABOLIC_FLOW;
            } else {
                throw std::invalid_argument("wrong configuration value: flow_mode");
            }
        }
        else {
            flowMode = UNIFORM_FLOW;
        }

        double backgroundPressure_Pa = simConf->doubleParameter("background_pressure_Pa");
        double gasVelocityX = simConf->doubleParameter("collision_gas_velocity_x_mscylinder-1");
        double collisionGasMass_Amu = simConf->doubleParameter("collision_gas_mass_amu");
        double collisionGasDiameter_nm = simConf->doubleParameter("collision_gas_diameter_nm");
        double backgroundTemperature_K = simConf->doubleParameter("background_temperature_K");

        double gateOpenTime = simConf->doubleParameter("gate_open_time_s");
        double gateOpenDuration = simConf->doubleParameter("gate_open_duration_s");
        double gateVoltage = simConf->doubleParameter("gate_voltage_V");
        double gradientStartTime = simConf->doubleParameter("gradient_start_time_s");
        double gradientVelocity = simConf->doubleParameter("gradient_ramp_velocity_V/ms")*1000;
        double gradientVoltage = simConf->doubleParameter("gradient_voltage_V");
        double gradientDuration = gradientVoltage / gradientVelocity;


        //Additional parameters needed for MD collision model
        std::string collisionGasIdentifier;
        std::vector<std::string> particleIdentifier;
        double collisionGasPolarizability_m3 = 0;
        double subIntegratorIntegrationTime_s = 0;
        double subIntegratorStepSize_s = 0;
        double collisionRadiusScaling = 0;
        double angleThetaScaling = 0;
        double spawnRadius_m = 0;
        bool saveTrajectory = false;
        if(collisionType==MD){
            collisionGasPolarizability_m3 = simConf->doubleParameter("collision_gas_polarizability_m3");
            collisionGasIdentifier = simConf->stringParameter("collision_gas_identifier");
            particleIdentifier = simConf->stringVectorParameter("particle_identifier");
            subIntegratorIntegrationTime_s = simConf->doubleParameter("sub_integrator_integration_time_s");
            subIntegratorStepSize_s = simConf->doubleParameter("sub_integrator_step_size_s");
            collisionRadiusScaling = simConf->doubleParameter("collision_radius_scaling");
            angleThetaScaling = simConf->doubleParameter("angle_theta_scaling");
            spawnRadius_m = simConf->doubleParameter("spawn_radius_m");
            saveTrajectory = simConf->boolParameter("save_trajectory");
        }

        // ======================================================================================

        //read potential array configuration of the trap =================================================
        std::filesystem::path confBasePath = simConf->confBasePath();

        FileIO::CSVReader csvReader = FileIO::CSVReader();
        std::vector<std::vector<std::string>> stringVector = std::vector<std::vector<std::string>>();
        std::string potentialConfFn = simConf->pathRelativeToConfFile(simConf->stringParameter("potential_configuration"));
        stringVector = csvReader.readCSVFile(potentialConfFn, ';');
        std::vector<std::string> PotentialArraysNames = csvReader.extractString(stringVector, 0);
        std::vector<double> DCVoltages = csvReader.extractDouble(stringVector, 1);
        std::vector<double> RFFactor = csvReader.extractDouble(stringVector, 2);
        std::vector<double> gradient = csvReader.extractDouble(stringVector, 3);
        std::vector<double> gate = csvReader.extractDouble(stringVector, 4);

        double paSpatialScale = simConf->doubleParameter("potential_array_scale");
        std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> PotentialArrays;
        //std::vector<std::string> PotentialArraysNames = simConf->stringVectorParameter("potential_arrays");
        for (const auto& paName: PotentialArraysNames) {
            std::filesystem::path paPath = confBasePath/paName;
            std::unique_ptr<ParticleSimulation::SimionPotentialArray> pa_pt =
                    std::make_unique<ParticleSimulation::SimionPotentialArray>(paPath, paSpatialScale);
            PotentialArrays.push_back(std::move(pa_pt));
        }

        double potentialScale = 1.0/10000.0;

        std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> flowField =
                simConf->readPotentialArrays("flow_field", paSpatialScale, true);

        std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> pressureField =
                simConf->readPotentialArrays("pressure_field", paSpatialScale, true);

        std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> temperatureField =
                simConf->readPotentialArrays("temperature_field", paSpatialScale, true);


        // defining simulation domain box (used for ion termination):
        std::array<std::array<double, 2>, 3> simulationDomainBoundaries{};
        if (simConf->isParameter("simulation_domain_boundaries_m")) {
            // get manual simulation domain boundaries from config file
            simulationDomainBoundaries = simConf->double3dBox("simulation_domain_boundaries_m");
        }
        else {
            // use minimum PA extent box as domain boundaries
            std::array<double, 6> minExtent = PotentialArrays[0]->getBounds();
            for (const auto& pa: PotentialArrays) {
                std::array<double, 6> paBounds = pa->getBounds();
                for (size_t i = 0; i<6; ++i) {
                    if (minExtent[i]<paBounds[i]) {
                        minExtent[i] = paBounds[i];
                    }
                }
            }
            simulationDomainBoundaries[0][0] = minExtent[0];
            simulationDomainBoundaries[0][1] = minExtent[1];
            simulationDomainBoundaries[1][0] = minExtent[2];
            simulationDomainBoundaries[1][1] = minExtent[3];
            simulationDomainBoundaries[2][0] = minExtent[4];
            simulationDomainBoundaries[2][1] = minExtent[5];
        }

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

        std::vector<std::string> integerParticleAttributesNames = {"chemical id"};
        auto additionalIntegerParamTFct = [](Core::Particle* particle) -> std::vector<int> {
            std::vector<int> result = {
                    particle->getIntegerAttribute(key_ChemicalIndex)
            };
            return result;
        };

        std::vector<std::string> auxParamNames = {"velocity x", "velocity y", "velocity z",
                                                  "ion velocity","kinetic energy (eV)", "effective Field (V/m)", "ion temperature (K)"};
        auto auxParamFct = [backgroundTemperature_K, collisionGasMass_Amu](Core::Particle* particle) -> std::vector<double> {
            double ionVelocity = particle->getVelocity().magnitude();
            double kineticEnergy_eV = 0.5*particle->getMass()*ionVelocity*ionVelocity*Core::JOULE_TO_EV;
            double ionTemperature = backgroundTemperature_K + (collisionGasMass_Amu*pow(ionVelocity,2))/3*1.381e-23;
            std::vector<double> result = {
                    particle->getVelocity().x(),
                    particle->getVelocity().y(),
                    particle->getVelocity().z(),
                    ionVelocity,
                    kineticEnergy_eV,
                    particle->getFloatAttribute("effectiveField"),
                    ionTemperature
            };
            if (!particle->isActive()) {
                result.push_back(particle->getSplatTime());
            }
            return result;
        };

        std::string hdf5Filename = projectName+"_trajectories.hd5";
        FileIO::TrajectoryHDF5Writer trajectoryWriter(hdf5Filename);
        trajectoryWriter.setParticleAttributes(integerParticleAttributesNames, additionalIntegerParamTFct);
        trajectoryWriter.setParticleAttributes(auxParamNames, auxParamFct);

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

        Core::Vector initCorner(startPosition_mm[0]/1000, startPosition_mm[1]/1000, startPosition_mm[2]/1000);
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

                if(collisionType==MD){
                    particle->setMolecularStructure(molecularStructureCollection.at(particleIdentifier[i]));
                    particle->setDiameter(particle->getMolecularStructure()->getDiameter());
                }

                particlesPtrs.push_back(particle.get());
                rsSim.addParticle(particle.get(), nParticlesTotal);
                particles.push_back(std::move(particle));
                trajectoryAdditionalParams.emplace_back(1);
                nParticlesTotal++;
            }
        }


        resultFilewriter.initFile(rsSimConf);
        // ======================================================================================


        // define trajectory integration parameters / functions =================================
        std::vector<double> totalFieldNow(PotentialArrays.size(), 0.0);

        auto paVoltageFct = [&PotentialArrays, &DCVoltages, &RFFactor, &gradient, &gate, &totalFieldNow,
                             omega, V_rf, gateOpenTime, gateOpenDuration, gateVoltage, gradientStartTime, gradientDuration, gradientVoltage, gradientVelocity]
                                     (double time){
            for(size_t i=0; i<PotentialArrays.size(); i++) {
                totalFieldNow[i] = DCVoltages[i] + sin(time*omega) * (V_rf * RFFactor[i]);

                if (time >= gateOpenTime && time <= (gateOpenTime + gateOpenDuration))
                    totalFieldNow[i] = totalFieldNow[i] + gate[i] * (gateVoltage);

                if (time >= gradientStartTime && time <= (gradientStartTime + gradientDuration)){
                    totalFieldNow[i] = totalFieldNow[i] + gradient[i] * ((time-gradientStartTime)*gradientVelocity);}
            }
        };


        auto accelerationFct = [&PotentialArrays, &totalFieldNow, potentialScale, spaceChargeFactor]
                (Core::Particle* particle, int /*particleIndex*/, SpaceCharge::FieldCalculator &scFieldCalculator,
                 double /*time*/, int /*timestep*/){
            Core::Vector fEfield(0, 0, 0);
            Core::Vector pos = particle->getLocation();
            double particleCharge = particle->getCharge();

            for(size_t i=0; i<PotentialArrays.size(); i++) {
                Core::Vector paField = PotentialArrays[i]->getField(pos.x(), pos.y(), pos.z());
                Core::Vector paEffectiveField = paField * totalFieldNow[i] * potentialScale;

                fEfield = fEfield + paEffectiveField;
            }

            particle->setFloatAttribute("effectiveField", fEfield.magnitude());

            if (Core::isDoubleEqual(spaceChargeFactor, 0.0)) {
                return (fEfield * particleCharge / particle->getMass());
            }
            else {
                Core::Vector spaceChargeForce =
                        scFieldCalculator.getEFieldFromSpaceCharge(*particle)*(spaceChargeFactor);
                return (((fEfield+spaceChargeForce)*particleCharge)/particle->getMass());
            }

        };


        auto timestepWriteFct =
                [&trajectoryWriter, &voltageWriter, trajectoryWriteInterval, &rsSim, &resultFilewriter, concentrationWriteInterval,
                        &totalFieldNow, &logger, &ionsInactive]
                        (
                                Integration::AbstractTimeIntegrator* /*integrator*/,
                                std::vector<Core::Particle*>& particles, double time, int timestep,
                                bool lastTimestep)
                {

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
                        logger->info("ts:{} time:{:.2e} splatted ions:{}",
                                     timestep, time, ionsInactive);
                        rsSim.logConcentrations(logger);
                        trajectoryWriter.writeTimestep(particles, time);
                    }
                };

        // Prepare ion start / stop tracker and ion start monitoring / ion termination functions
        ParticleSimulation::ParticleStartSplatTracker startSplatTracker;
        auto particleStartMonitoringFct = [&startSplatTracker](Core::Particle* particle, double time) {
            startSplatTracker.particleStart(particle, time);
        };

        auto otherActionsFct = [&simulationDomainBoundaries, &ionsInactive, &PotentialArrays, &V_rf, &startSplatTracker](
                Core::Vector& newPartPos, Core::Particle* particle,
                int /*particleIndex*/,  double time, int /*timestep*/) {
            //Core::Vector pos = particle->getLocation();
            if (newPartPos.x()<=simulationDomainBoundaries[0][0] ||
                newPartPos.x()>=simulationDomainBoundaries[0][1] ||
                newPartPos.y()<=simulationDomainBoundaries[1][0] ||
                newPartPos.y()>=simulationDomainBoundaries[1][1] ||
                newPartPos.z()<=simulationDomainBoundaries[2][0] ||
                newPartPos.z()>=simulationDomainBoundaries[2][1] ) {
                particle->setActive(false);
                particle->setSplatTime(time);
                startSplatTracker.particleSplat(particle, time);
                ionsInactive++;
            }
            if (PotentialArrays[0]->isElectrode(newPartPos.x(), newPartPos.y(), newPartPos.z())) {
                particle->setActive(false);
                particle->setSplatTime(time);
                startSplatTracker.particleSplat(particle, time);
                ionsInactive++;
            }
            Core::Vector PartVelocity(particle->getVelocity());
            if (Core::isDoubleEqual(V_rf, 0.0)) {
                double boundary = 0.001;
                if (newPartPos.y() > boundary) {
                    double new_y = -(boundary - (newPartPos.y() - boundary));
                    Core::Vector newPos(newPartPos.x(), new_y, newPartPos.z());
                    newPartPos = newPos;
                    Core::Vector newPartVelocity(PartVelocity.x(), -PartVelocity.y(), PartVelocity.z());
                    particle->setVelocity(newPartVelocity);
                }
                if (newPartPos.y() < -boundary) {
                    double new_y = +(boundary + (newPartPos.y() + boundary));
                    Core::Vector newPos(newPartPos.x(), new_y, newPartPos.z());
                    newPartPos = newPos;
                    Core::Vector newPartVelocity(PartVelocity.x(), -PartVelocity.y(), PartVelocity.z());
                    particle->setVelocity(newPartVelocity);
                }
                if (newPartPos.z() > boundary) {
                    double new_z = -(boundary - (newPartPos.z() - boundary));
                    Core::Vector newPos(newPartPos.x(), newPartPos.y(), new_z);
                    newPartPos = newPos;
                    Core::Vector newPartVelocity(PartVelocity.x(), PartVelocity.y(), -PartVelocity.z());
                    particle->setVelocity(newPartVelocity);
                }
                if (newPartPos.z() < -boundary) {
                    double new_z = +(boundary + (newPartPos.z() + boundary));
                    Core::Vector newPos(newPartPos.x(), newPartPos.y(), new_z);
                    newPartPos = newPos;
                    Core::Vector newPartVelocity(PartVelocity.x(), PartVelocity.y(), -PartVelocity.z());
                    particle->setVelocity(newPartVelocity);
                }
            }
        };

        //define / gas interaction /  collision model:
        std::unique_ptr<CollisionModel::AbstractCollisionModel> collisionModelPtr;
        if (collisionType==SDS || collisionType==HS || collisionType==MD) {
            // prepare static pressure and temperature functions

            auto pressureFct = CollisionModel::getVariableScalarFunction(*pressureField[0]);
            auto backgroundTemperatureFct = CollisionModel::getVariableScalarFunction(*temperatureField[0]);
            auto velocityFct = CollisionModel::getVariableVectorFunction(*flowField[0], *flowField[1], *flowField[2]);


            if (collisionType==SDS){
                std::unique_ptr<CollisionModel::StatisticalDiffusionModel> collisionModel =
                        std::make_unique<CollisionModel::StatisticalDiffusionModel>(
                                pressureFct,
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
            else if (collisionType==HS){
                std::unique_ptr<CollisionModel::HardSphereModel> collisionModel =
                        std::make_unique<CollisionModel::HardSphereModel>(
                                pressureFct,
                                velocityFct,
                                backgroundTemperatureFct,
                                collisionGasMass_Amu,
                                collisionGasDiameter_nm*1e-9,
                                nullptr);
                collisionModelPtr = std::move(collisionModel);
            }
            else if (collisionType==MD){
                CollisionModel::MDForceField_LJ12_6 forceField(collisionGasPolarizability_m3);
                auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(forceField);
                std::unique_ptr<CollisionModel::MDInteractionsModel> collisionModel =
                        std::make_unique<CollisionModel::MDInteractionsModel>(
                                pressureFct,
                                velocityFct,
                                backgroundTemperatureFct,
                                collisionGasMass_Amu,
                                collisionGasDiameter_nm*1e-9,
                                collisionGasIdentifier,
                                subIntegratorIntegrationTime_s,
                                subIntegratorStepSize_s,
                                collisionRadiusScaling,
                                angleThetaScaling,
                                spawnRadius_m,
                                std::move(forceFieldPtr),
                                molecularStructureCollection);

                if (saveTrajectory){
                    unsigned int saveTrajectoryStartTimeStep = simConf->unsignedIntParameter("trajectory_start_time_step");
                    double trajectoryDistance_m = simConf->doubleParameter("trajectory_distance_m");
                    collisionModel->setTrajectoryWriter(projectName+"_md_trajectories.txt",
                                                        trajectoryDistance_m, saveTrajectoryStartTimeStep);
                }
            }
        }
        else if (collisionType==NO_COLLISION) {
            collisionModelPtr = nullptr;
        }

        //define reaction simulation functions:
        auto particlesHasReactedFct = [&collisionModelPtr, &substanceIndices](RS::ReactiveParticle* particle){
            //we had a reaction event: Count it (access to total counted value has to be synchronized)
            if (collisionModelPtr != nullptr) {
                collisionModelPtr->initializeModelParticleParameters(*particle);
            }
            int substIndex = substanceIndices.at(particle->getSpecies());
            particle->setIntegerAttribute(key_ChemicalIndex, substIndex);
        };

        auto reactionConditionsFct = [backgroundTemperature_K, backgroundPressure_Pa]
                (RS::ReactiveParticle* particle, double /*time*/)->RS::ReactionConditions{
            RS::ReactionConditions reactionConditions = RS::ReactionConditions();

            reactionConditions.temperature = backgroundTemperature_K;
            reactionConditions.electricField = particle->getFloatAttribute("effectiveField");
            reactionConditions.pressure = backgroundPressure_Pa;
            return reactionConditions;
        };

        //init trajectory simulation object:
        Integration::ParallelVerletIntegrator verletIntegrator(
                particlesPtrs,
                accelerationFct, timestepWriteFct, otherActionsFct, particleStartMonitoringFct,
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
