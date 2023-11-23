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
#include "CollisionModel_MultiCollisionModel.hpp"
#include <iostream>
#include <cmath>

const std::string key_ChemicalIndex = "keyChemicalIndex";
enum CollisionType {SDS, HS, MD, MULTI_HS, HS_MD, NO_COLLISION};

int main(int argc, const char * argv[]) {

    try {
        // open configuration, parse configuration file =========================================
        AppUtils::CommandlineParser cmdLineParser(argc, argv, "BT-RS-TWIMSSim", "TWIMS Simulation with trajectories and chemistry", true);
        std::string projectName = cmdLineParser.resultName();
        AppUtils::logger_ptr logger = cmdLineParser.logger();

        std::string confFileName = cmdLineParser.confFileName();
        AppUtils::simConf_ptr simConf = cmdLineParser.simulationConfiguration();

        // optionally setting random generator seed manually (for debugging / reproduction purposes):
        if (simConf->isParameter("random_seed")) {
            unsigned int randomSeed = simConf->unsignedIntParameter("random_seed");
            Core::globalRandomGeneratorPool->setSeedForElements(randomSeed);
        }

        std::string simulationMode = simConf->stringParameter("Mode");

        std::vector<unsigned int> nParticles = simConf->unsignedIntVectorParameter("n_particles");
        unsigned int nSteps = simConf->unsignedIntParameter("sim_time_steps");
        double dt_s = simConf->doubleParameter("dt_s");
        int concentrationWriteInterval = simConf->intParameter("concentrations_write_interval");
        int trajectoryWriteInterval = simConf->intParameter("trajectory_write_interval");
        double spaceChargeFactor = simConf->doubleParameter("space_charge_factor");
        double waveAmplitude = simConf->doubleParameter("wave_amplitude_V");
        double waveFrequency = simConf->doubleParameter("wave_frequency_hz");

        double f_rf = simConf->doubleParameter("confining_RF_frequency_Hz");
        double omega = f_rf*2.0*M_PI;
        double V_rf = simConf->doubleParameter("confining_RF_amplitude_V");

        //geometric parameters:
        double startWidthX_m = simConf->doubleParameter("start_width_x_mm")/1000.0;
        double startWidthY_m = simConf->doubleParameter("start_width_y_mm")/1000.0;
        double startWidthZ_m = simConf->doubleParameter("start_width_z_mm")/1000.0;

        double guardVoltage = simConf->doubleParameter("guard_bias_V");

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
        else if (collisionTypeStr=="multi_HS") {
            collisionType = MULTI_HS;
        }
        else if (collisionTypeStr=="HS_MD") {
            collisionType = HS_MD;
        }
        else if (collisionTypeStr=="none") {
            collisionType = NO_COLLISION;
        }
        else {
            throw std::invalid_argument("wrong configuration value: collision_model_type");
        }

        //sets a flag on whether a collision type includes the MD model
        bool usesMDmodel = (collisionType==MD || collisionType==HS_MD);

        std::vector<double> backgroundPartialPressures_Pa = simConf->doubleVectorParameter(
                "background_partial_pressures_Pa");
        std::vector<double> collisionGasMasses_Amu = simConf->doubleVectorParameter("collision_gas_masses_amu");
        std::vector<double> collisionGasDiameters_angstrom = simConf->doubleVectorParameter(
                "collision_gas_diameters_angstrom");
        double backgroundTemperature_K = simConf->doubleParameter("background_temperature_K");

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

        //Additional parameters needed for MD collision model
        std::vector<std::string> collisionGasIdentifier;
        std::vector<std::string> particleIdentifier;
        std::vector<double> collisionGasPolarizability_m3;
        double subIntegratorIntegrationTime_s = 0;
        double subIntegratorStepSize_s = 0;
        double collisionRadiusScaling = 0;
        double angleThetaScaling = 0;
        double spawnRadius_m = 0;
        double trajectoryDistance_m = 0;
        bool saveTrajectory = false;
        int saveTrajectoryStartTimeStep = 0;
        if(usesMDmodel== true){
            collisionGasPolarizability_m3 = simConf->doubleVectorParameter("collision_gas_polarizability_m3");
            collisionGasIdentifier = simConf->stringVectorParameter("collision_gas_identifier");
            particleIdentifier = simConf->stringVectorParameter("particle_identifier");
            subIntegratorIntegrationTime_s = simConf->doubleParameter("sub_integrator_integration_time_s");
            subIntegratorStepSize_s = simConf->doubleParameter("sub_integrator_step_size_s");
            collisionRadiusScaling = simConf->doubleParameter("collision_radius_scaling");
            angleThetaScaling = simConf->doubleParameter("angle_theta_scaling");
            spawnRadius_m = simConf->doubleParameter("spawn_radius_m");
            saveTrajectory = simConf->boolParameter("save_trajectory");
            trajectoryDistance_m = simConf->doubleParameter("trajectory_distance_m");
            saveTrajectoryStartTimeStep = simConf->intParameter("trajectory_start_time_step");
        }

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
        if(usesMDmodel == true){
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
                                                  "ion velocity","kinetic energy (eV)", "effective Field (V/m)", "splat time"};
        auto auxParamFct = [](Core::Particle* particle) -> std::vector<double> {
            double ionVelocity = particle->getVelocity().magnitude();
            double kineticEnergy_eV = 0.5*particle->getMass()*ionVelocity*ionVelocity*Core::JOULE_TO_EV;
            std::vector<double> result = {
                    particle->getVelocity().x(),
                    particle->getVelocity().y(),
                    particle->getVelocity().z(),
                    ionVelocity,
                    kineticEnergy_eV,
                    particle->getFloatAttribute("effectiveField"),
                    //particle->getSplatTime()
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

        Core::Vector initCorner(2.0e-03, -startWidthY_m/2.0, -startWidthZ_m/2.0);
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

                if(usesMDmodel == true){
                    particle->setMolecularStructure(molecularStructureCollection.at(particleIdentifier[i]));
                    particle->setDiameter(particle->getMolecularStructure()->getDiameter());
                }

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
                wavePeriod, waveAmplitude, omega, V_rf, &totalFieldNow, simulationMode, guardVoltage](double time){
            if (simulationMode == "TWIMS") {
                for(size_t i=0; i<potentialArrays.size()-2; i++) {
                    double phase = std::fmod(time, wavePeriod) / wavePeriod;
                    double shiftedphase = std::fmod(phase + phaseShift[i], 1.0);
                    totalFieldNow[i]=(WaveForm.getInterpolatedValue(shiftedphase) * waveAmplitude);
                }

                totalFieldNow[potentialArrays.size()-2] = sin(time*omega) * V_rf;
                totalFieldNow[potentialArrays.size()-1] = -totalFieldNow[potentialArrays.size()-2];
            }

            if (simulationMode == "SLIM") {
                for(size_t i=0; i<potentialArrays.size()-3; i++) {
                    double phase = std::fmod(time, wavePeriod) / wavePeriod;
                    double shiftedphase = std::fmod(phase + phaseShift[i], 1.0);
                    totalFieldNow[i]=(WaveForm.getInterpolatedValue(shiftedphase) * waveAmplitude);
                }

                totalFieldNow[potentialArrays.size()-3] = sin(time*omega) * V_rf;
                totalFieldNow[potentialArrays.size()-2] = -totalFieldNow[potentialArrays.size()-2];

                totalFieldNow[potentialArrays.size()-1] = guardVoltage;
            }
        };

        auto accelerationFct =
                [&potentialArrays, &totalFieldNow, potentialScale, spaceChargeFactor]
                        (Core::Particle* particle, int /*particleIndex*/, SpaceCharge::FieldCalculator &scFieldCalculator,
                         double /*time*/, int timestep){
                    Core::Vector fEfield(0, 0, 0);
                    Core::Vector pos = particle->getLocation();
                    double particleCharge = particle->getCharge();

                    for(size_t i=0; i<potentialArrays.size(); i++) {
                        Core::Vector paField = potentialArrays[i]->getField(pos.x(), pos.y(), pos.z());
                        Core::Vector paEffectiveField = paField * totalFieldNow[i] * potentialScale;

                        fEfield = fEfield+paEffectiveField;
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
                        logger->info("ts:{} time:{:.2e} splatted ions:{}",
                                     timestep, time, ionsInactive);
                        rsSim.logConcentrations(logger);
                        trajectoryWriter.writeTimestep(particles, time);
                    }
                };


        auto otherActionsFct = [&ionsInactive, &potentialArrays, &V_rf](
                Core::Vector& newPartPos, Core::Particle* particle,
                int /*particleIndex*/,  double time, int /*timestep*/) {
            //Core::Vector pos = particle->getLocation();
            if (newPartPos.x() < 0.00001) {
                particle->setActive(false);
                ionsInactive++;
            }
            if (newPartPos.x() > 0.11200) {
                particle->setActive(false);
                particle->setSplatTime(time);
                ionsInactive++;
            }
            if (potentialArrays[0]->isElectrode(newPartPos.x(), newPartPos.y(), newPartPos.z())) {
                particle->setActive(false);
                particle->setSplatTime(time);
                ionsInactive++;
            }
            if (V_rf == 0.0) {
                double boundary = 0.001;
                if (newPartPos.y() > boundary) {
                    double new_y = -(boundary - (newPartPos.y() - boundary));
                    Core::Vector newPos(newPartPos.x(), new_y, newPartPos.z());
                    newPartPos = newPos;
                }
                if (newPartPos.y() < -boundary) {
                    double new_y = +(boundary + (newPartPos.y() + boundary));
                    Core::Vector newPos(newPartPos.x(), new_y, newPartPos.z());
                    newPartPos = newPos;
                }
                if (newPartPos.z() > boundary) {
                    double new_z = -(boundary - (newPartPos.z() - boundary));
                    Core::Vector newPos(newPartPos.x(), newPartPos.y(), new_z);
                    newPartPos = newPos;
                }
                if (newPartPos.z() < -boundary) {
                    double new_z = +(boundary + (newPartPos.z() + boundary));
                    Core::Vector newPos(newPartPos.x(), newPartPos.y(), new_z);
                    newPartPos = newPos;
                }
            }


            //method 1: box method
            /*double radius = 0.002;
            double box_length = sqrt(2) * radius;

            double new_y = newPartPos.y() + (box_length/2);
            double new_z = newPartPos.z() + (box_length/2);
            new_y = fmod(new_y, box_length);
            new_z = fmod(new_z, box_length);

            Core::Vector newPos(newPartPos.x(), (new_y - (box_length/2)), (new_z - (box_length/2)));
            newPartPos = newPos;
            */
            //method 2: circle radius method
            /*Core::Vector currentPoint = newPartPos;
            Core::Vector v(particle->getVelocity().x()*dt_s, particle->getVelocity().y()*dt_s, particle->getVelocity().z()*dt_s);
            Core::Vector nextPoint(currentPoint.x()+v.x(), currentPoint.y()*v.y(), currentPoint.z()+v.z());
            double radius = 0.001;

            //if (pow(nextPoint.y(),2)+ pow(nextPoint.z(),2) > pow(0.001, 2)) {
            while (pow(nextPoint.y(),2) + pow(nextPoint.z(),2) > pow(radius, 2)) {
                double a = (pow(v.y(),2)+ pow(v.z(),2));
                double b = 2*v.y()*currentPoint.y()+2*v.z()*currentPoint.z();
                double c = pow(currentPoint.y(),2) + pow(currentPoint.z(),2) - pow(radius,2);

                double r1 = (-(b) + sqrt(pow(b,2)-4*a*c)) / 2*a;
                double r2 = (-(b) - sqrt(pow(b,2)-4*a*c)) / 2*a;

                double r = 0;
                if (r1 > 0 && r1 < 1) {
                    r = r1;
                }
                else {
                    r = r2;
                }
                Core::Vector mirrorPoint(currentPoint.x()+r*v.x(), -(currentPoint.y()+r*v.y()), -(currentPoint.z()+r*v.z()));
                currentPoint = mirrorPoint;
                Core::Vector nextPoint(mirrorPoint.x()+(1-r)*v.x(), mirrorPoint.y()+(1-r)*v.y(), mirrorPoint.z()+(1-r)*v.z());
                v = (1-r) * v;
            }
            newPartPos = nextPoint;*/
        };

        //define / gas interaction /  collision model:
        std::unique_ptr<CollisionModel::AbstractCollisionModel> collisionModelPtr;
        if (collisionType==SDS || collisionType==HS || collisionType==MD || collisionType==MULTI_HS || collisionType==HS_MD) {
            // prepare static pressure and temperature functions

            if (collisionType==SDS){
                std::unique_ptr<CollisionModel::StatisticalDiffusionModel> collisionModel =
                        std::make_unique<CollisionModel::StatisticalDiffusionModel>(
                                backgroundPartialPressures_Pa[0],
                                backgroundTemperature_K,
                                collisionGasMasses_Amu[0],
                                collisionGasDiameters_m[0]);

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
                                backgroundPartialPressures_Pa[0],
                                backgroundTemperature_K,
                                collisionGasMasses_Amu[0],
                                collisionGasDiameters_m[0],
                                nullptr);
                collisionModelPtr = std::move(collisionModel);
            }
            else if (collisionType==MD){
                //prepare multimodel with multiple MD models (one per collision gas)
                std::vector<std::unique_ptr<CollisionModel::AbstractCollisionModel>> mdModels;
                for (std::size_t i = 0; i<nBackgroundGases; ++i) {
                    auto mdModel = std::make_unique<CollisionModel::MDInteractionsModel>(
                            backgroundPartialPressures_Pa[i],
                            backgroundTemperature_K,
                            collisionGasMasses_Amu[i],
                            collisionGasDiameters_m[i],
                            collisionGasPolarizability_m3[i],
                            collisionGasIdentifier[i],
                            subIntegratorIntegrationTime_s,
                            subIntegratorStepSize_s,
                            collisionRadiusScaling,
                            angleThetaScaling,
                            spawnRadius_m,
                            molecularStructureCollection);
                    if (saveTrajectory){
                        mdModel->setTrajectoryWriter(projectName+"_md_trajectories.txt",
                                                     trajectoryDistance_m, saveTrajectoryStartTimeStep);
                    }
                    mdModels.emplace_back(std::move(mdModel));
                }

                std::unique_ptr<CollisionModel::MultiCollisionModel> collisionModel =
                        std::make_unique<CollisionModel::MultiCollisionModel>(std::move(mdModels));

                collisionModelPtr = std::move(collisionModel);
            }
            else if (collisionType==MULTI_HS) {
                //prepare multiple Hard Sphere models (one per collision gas)
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
            }
            else if (collisionType==HS_MD) {
                std::vector<std::unique_ptr<CollisionModel::AbstractCollisionModel>> hsmdModels;
                for (std::size_t i=0; i<nBackgroundGases; ++i) {
                    if ((i%2)==0) {
                        auto cllModel = std::make_unique<CollisionModel::HardSphereModel>(
                                backgroundPartialPressures_Pa[i],
                                backgroundTemperature_K,
                                collisionGasMasses_Amu[i],
                                collisionGasDiameters_m[i]);
                        hsmdModels.emplace_back(std::move(cllModel));
                    }
                    else {
                        auto cllModel = std::make_unique<CollisionModel::MDInteractionsModel>(
                                backgroundPartialPressures_Pa[i],
                                backgroundTemperature_K,
                                collisionGasMasses_Amu[i],
                                collisionGasDiameters_m[i],
                                collisionGasPolarizability_m3[i],
                                collisionGasIdentifier[i],
                                subIntegratorIntegrationTime_s,
                                subIntegratorStepSize_s,
                                collisionRadiusScaling,
                                angleThetaScaling,
                                spawnRadius_m,
                                molecularStructureCollection);
                        if (saveTrajectory){
                            cllModel->setTrajectoryWriter(projectName+"_md_trajectories.txt",
                                                          trajectoryDistance_m, saveTrajectoryStartTimeStep);
                        }
                        hsmdModels.emplace_back(std::move(cllModel));
                    }
                }
                std::cout << "Collision model loop finished." << std::endl;
                std::unique_ptr<CollisionModel::MultiCollisionModel> collisionModel =
                        std::make_unique<CollisionModel::MultiCollisionModel>(std::move(hsmdModels));
                std::cout << "Collision model made." << std::endl;

                collisionModelPtr = std::move(collisionModel);
                std::cout << "Collision model done." << std::endl;
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

        auto reactionConditionsFct = [backgroundTemperature_K, backgroundPartialPressures_Pa]
                (RS::ReactiveParticle* particle, double /*time*/)->RS::ReactionConditions{
            RS::ReactionConditions reactionConditions = RS::ReactionConditions();

            reactionConditions.temperature = backgroundTemperature_K;
            reactionConditions.electricField = particle->getFloatAttribute("effectiveField");
            reactionConditions.pressure = backgroundPartialPressures_Pa[0];
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
