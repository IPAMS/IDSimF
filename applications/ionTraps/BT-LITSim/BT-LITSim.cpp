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
 BT-LitSim.cpp

 Ion trajectory simulation (including space charge and hard sphere collisions) in an RF Linear ion trap with
 arbitrary electrode geometry given by potential arrays and analytical potential along z (long) axis

 ****************************/


#include "BTree_particle.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_scalar_writer.hpp"
#include "PSim_util.hpp"
#include "PSim_verletIntegrator.hpp"
#include "PSim_parallelVerletIntegrator.hpp"
#include "PSim_sampledWaveform.hpp"
#include "PSim_particleStartSplatTracker.hpp"
#include "PSim_math.hpp"
#include "PSim_averageChargePositionWriter.hpp"
#include "PSim_inductionCurrentWriter.hpp"
#include "PSim_simionPotentialArray.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_ionDefinitionReading.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include <iostream>
#include <vector>
#include <ctime>
#include <filesystem>

enum IntegratorMode {VERLET,PARALLEL_VERLET};
enum RfAmplitudeMode {STATIC_RF,RAMPED_RF};
enum ExciteMode {RECTPULSE,SWIFT};
enum FftWriteMode {UNRESOLVED,MASS_RESOLVED};

const std::string key_spaceCharge_x = "keySpaceChargeX";
const std::string key_spaceCharge_y = "keySpaceChargeY";
const std::string key_spaceCharge_z = "keySpaceChargeZ";
const std::string key_trapForce_x = "keyTrapForceX";
const std::string key_trapForce_y = "keyTrapForceY";
const std::string key_trapForce_z = "keyTrapForceZ";

int main(int argc, const char * argv[]) {

    try {
        // read configuration file ======================================================================
        if (argc<=2) {
            std::cout << "Run abort: No run configuration or project name given." << std::endl;
            return EXIT_FAILURE;
        }
        std::string projectName = argv[2];
        std::cout << projectName << std::endl;
        auto logger = AppUtils::createLogger(projectName+".log");

        std::string confFileName = argv[1];
        AppUtils::SimulationConfiguration simConf(confFileName, logger);
        std::filesystem::path confBasePath = simConf.confBasePath();


        // read basic simulation parameters =============================================================
        std::string integratorMode_str = simConf.stringParameter("integrator_mode");
        IntegratorMode integratorMode;
        if (integratorMode_str=="parallel_verlet") {
            integratorMode = PARALLEL_VERLET;
        }
        else {
            throw std::invalid_argument("wrong configuration value: integrator mode");
        }

        unsigned int timeSteps = simConf.unsignedIntParameter("sim_time_steps");
        int trajectoryWriteInterval = simConf.intParameter("trajectory_write_interval");
        int fftWriteInterval = simConf.intParameter("fft_write_interval");
        double dt = simConf.doubleParameter("dt");

        std::string fftWriteMode_str = simConf.stringParameter("fft_write_mode");
        FftWriteMode fftWriteMode = UNRESOLVED;
        if (fftWriteMode_str=="unresolved") {
            fftWriteMode = UNRESOLVED;
        }
        else if (fftWriteMode_str=="mass_resolved") {
            fftWriteMode = MASS_RESOLVED;
        }

        //read potential array configuration of the trap =================================================
        double paSpatialScale = simConf.doubleParameter("potential_array_scaling");
        std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> potentialArrays;
        std::vector<std::string> potentialArraysNames = simConf.stringVectorParameter("potential_arrays");
        for (const auto& paName: potentialArraysNames) {
            std::filesystem::path paPath = confBasePath/paName;
            std::unique_ptr<ParticleSimulation::SimionPotentialArray> pa_pt =
                    std::make_unique<ParticleSimulation::SimionPotentialArray>(paPath, paSpatialScale);
            potentialArrays.push_back(std::move(pa_pt));
        }

        double axialPotentialCenter = simConf.doubleParameter("axial_potential_center");
        double axialPotentialFloowWidth = simConf.doubleParameter("axial_potential_floor_width");
        double axialPotentialGradient = simConf.doubleParameter("axial_potential_gradient_V/m");

        // SIMION fast adjust PAs use 10000 as normalized potential value, thus we have to scale everything with 1/10000
        double potentialScale = 1.0/10000.0;
        std::vector<double> potentialsFactorsDc = simConf.doubleVectorParameter("dc_potentials", potentialScale);
        std::vector<double> potentialFactorsRf = simConf.doubleVectorParameter("rf_potential_factors", potentialScale);
        std::vector<double> potentialFactorsExcite = simConf.doubleVectorParameter("excite_potential_factors",
                potentialScale);
        std::vector<double> detectionPAFactorsRaw = simConf.doubleVectorParameter("detection_potential_factors");
        std::vector<ParticleSimulation::SimionPotentialArray*> detectionPAs;

        std::vector<double> detectionPAFactors;
        for (size_t i = 0; i<detectionPAFactorsRaw.size(); ++i) {
            double paKey = detectionPAFactorsRaw[i];
            if (paKey!=0.0) {
                detectionPAs.emplace_back(potentialArrays[i].get());
                detectionPAFactors.emplace_back(paKey);
            }
        }

        // defining simulation domain box (used for ion termination):
        std::array<std::array<double, 2>, 3> simulationDomainBoundaries;
        if (simConf.isParameter("simulation_domain_boundaries")) {
            // get manual simulation domain boundaries from config file
            simulationDomainBoundaries = simConf.double3dBox("simulation_domain_boundaries");
        }
        else {
            // use minimum PA extent box as domain boundaries
            std::array<double, 6> minExtent = potentialArrays[0]->getBounds();
            for (const auto& pa: potentialArrays) {
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

        //read physical configuration ===================================================================
        double backgroundPressure = simConf.doubleParameter("background_pressure_Pa");
        double backgroundTemperature = simConf.doubleParameter("background_temperature_K");
        double spaceChargeFactor = simConf.doubleParameter("space_charge_factor");
        double collisionGasMassAmu = simConf.doubleParameter("collision_gas_mass_amu");
        double collisionGasDiameterM = simConf.doubleParameter("collision_gas_diameter_angstrom")*1e-10;


        //read rf configuration =========================================================================
        double f_rf = simConf.doubleParameter("f_rf"); //RF frequency 1e6;
        double omega = f_rf*2.0*M_PI; //RF angular frequencyf_rf* 2.0 * M_PI;

        RfAmplitudeMode rfMode;
        std::vector<double> V_0_ramp;
        double V_0 = 0.0;
        if (simConf.isParameter("V_rf_start")) {
            rfMode = RAMPED_RF;
            double V_rf_start = simConf.doubleParameter("V_rf_start");
            double V_rf_end = simConf.doubleParameter("V_rf_end");
            V_0_ramp = ParticleSimulation::linspace(V_rf_start, V_rf_end, timeSteps);
        }
        else {
            rfMode = STATIC_RF;
            V_0 = simConf.doubleParameter("V_rf");
        }
        std::vector<double> V_rf_export;


        //read excitation / swift configuration ========================================================
        ExciteMode exciteMode;
        std::unique_ptr<ParticleSimulation::SampledWaveform> swiftWaveForm;
        double excitePulseLength = 0.0;
        if (simConf.isParameter("excite_waveform_csv_file")) {
            exciteMode = SWIFT;
            std::string swiftFileName = simConf.stringParameter("excite_waveform_csv_file");
            swiftWaveForm = std::make_unique<ParticleSimulation::SampledWaveform>(swiftFileName);
            if (!swiftWaveForm->good()) {
                logger->error("swift transient file not accessible");
                return (0);
            }
        }
        else {
            exciteMode = RECTPULSE;
            excitePulseLength = simConf.doubleParameter("excite_pulse_length");
        }
        double excitePulsePotential = simConf.doubleParameter("excite_pulse_potential");

        //read ion configuration =======================================================================
        std::vector<std::unique_ptr<BTree::Particle>> particles;
        std::vector<BTree::Particle*> particlePtrs;
        AppUtils::readIonDefinition(particles, particlePtrs, simConf);

        // init additional ion parameters:
        for (const auto& particle: particles) {
            particle->setFloatAttribute(key_trapForce_x, 0.0);
            particle->setFloatAttribute(key_trapForce_y, 0.0);
            particle->setFloatAttribute(key_trapForce_z, 0.0);
            particle->setFloatAttribute(key_spaceCharge_x, 0.0);
            particle->setFloatAttribute(key_spaceCharge_y, 0.0);
            particle->setFloatAttribute(key_spaceCharge_z, 0.0);
        }

        // define functions for the trajectory integration ==================================================
        std::size_t ionsInactive = 0;
        auto trapFieldFunction =
                [exciteMode, rfMode, excitePulseLength, excitePulsePotential, omega, &swiftWaveForm, &V_0, &V_0_ramp,
                        &potentialArrays, &potentialsFactorsDc, &potentialFactorsRf, &potentialFactorsExcite]
                        (BTree::Particle* particle, int /*particleIndex*/, auto& /*tree*/, double time, unsigned int timestep)
                        -> Core::Vector {

                    Core::Vector pos = particle->getLocation();
                    double particleCharge = particle->getCharge();
                    double excitePotential = 0;
                    if (exciteMode==RECTPULSE) {
                        if (time<excitePulseLength) {
                            excitePotential = excitePulsePotential;
                        }
                    }
                    else {
                        excitePotential = swiftWaveForm->getValue(timestep)*excitePulsePotential;
                    }

                    if (rfMode==RAMPED_RF) {
                        V_0 = V_0_ramp[timestep];
                    }

                    Core::Vector fEfield(0, 0, 0);
                    double V_t = V_0*cos(omega*time);

                    for (size_t i = 0; i<potentialArrays.size(); ++i) {
                        Core::Vector paField = potentialArrays[i]->getField(pos.x(), pos.y(), pos.z());

                        Core::Vector paEffectiveField =
                                paField*potentialsFactorsDc.at(i)+
                                        paField*potentialFactorsRf.at(i)*V_t+
                                        paField*potentialFactorsExcite.at(i)*excitePotential;

                        fEfield = fEfield+paEffectiveField;
                    }

                    return fEfield*particleCharge;
                };

        double axialPotential_upperStart = axialPotentialCenter+axialPotentialFloowWidth;
        double axialPotential_lowerStart = axialPotentialCenter-axialPotentialFloowWidth;
        auto accelerationFunctionLIT_parallel =
                [spaceChargeFactor, &trapFieldFunction,
                        axialPotential_lowerStart, axialPotential_upperStart, axialPotentialGradient](
                        BTree::Particle* particle, int particleIndex,
                        auto& tree, double time, unsigned int timestep) -> Core::Vector {

                    Core::Vector pos = particle->getLocation();
                    double particleCharge = particle->getCharge();

                    Core::Vector rfForce = trapFieldFunction(particle, particleIndex, tree, time, timestep);
                    double z_force = 0;
                    if (pos.z()>axialPotential_upperStart) {
                        z_force = (axialPotential_upperStart-pos.z())*axialPotentialGradient*particleCharge;
                    }
                    else if (pos.z()<axialPotential_lowerStart) {
                        z_force = (axialPotential_lowerStart-pos.z())*axialPotentialGradient*particleCharge;
                    }

                    Core::Vector axialForce = {0.0, 0.0, z_force};

                    Core::Vector spaceChargeForce(0, 0, 0);
                    if (spaceChargeFactor>0) {
                        spaceChargeForce =
                                tree.computeEFieldFromTree(*particle)*(particleCharge*spaceChargeFactor);
                    }

                    //update the additional parameters for writing them later to the trajectory:
                    particle->setFloatAttribute(key_trapForce_x, rfForce.x());
                    particle->setFloatAttribute(key_trapForce_y, rfForce.y());
                    particle->setFloatAttribute(key_trapForce_z, rfForce.z());
                    particle->setFloatAttribute(key_spaceCharge_x, spaceChargeForce.x());
                    particle->setFloatAttribute(key_spaceCharge_y, spaceChargeForce.y());
                    particle->setFloatAttribute(key_spaceCharge_z, spaceChargeForce.z());

                    return ((rfForce+axialForce+spaceChargeForce)/particle->getMass());
                };

        // Prepare ion start / stop tracker and ion start monitoring / ion termination functions
        ParticleSimulation::ParticleStartSplatTracker startSplatTracker;
        auto particleStartMonitoringFct = [&startSplatTracker](BTree::Particle* particle, double time) {
            startSplatTracker.particleStart(particle, time);
        };

        auto otherActionsFunctionQIT = [&simulationDomainBoundaries, &ionsInactive, &potentialArrays, &startSplatTracker](
                Core::Vector& newPartPos, BTree::Particle* particle,
                int /*particleIndex*/, auto& /*tree*/, double time, int /*timestep*/) {
            // if the ion is out of the boundary box or ends up in an electrode:
            // Terminate the ion
            // (since all potential arrays of the simulation define the basis functions of a linear combination,
            // the electrode geometry has to be the same in all electrodes, thus check only the first one)
            if (newPartPos.x()<=simulationDomainBoundaries[0][0] ||
                    newPartPos.x()>=simulationDomainBoundaries[0][1] ||
                    newPartPos.y()<=simulationDomainBoundaries[1][0] ||
                    newPartPos.y()>=simulationDomainBoundaries[1][1] ||
                    newPartPos.z()<=simulationDomainBoundaries[2][0] ||
                    newPartPos.z()>=simulationDomainBoundaries[2][1] ||
                    potentialArrays.at(0)->isElectrode(newPartPos.x(), newPartPos.y(), newPartPos.z())) {
                particle->setActive(false);
                particle->setSplatTime(time);
                startSplatTracker.particleSplat(particle, time);
                ionsInactive++;
            }
        };

        //prepare file writers and data writing functions ==============================================================================
        auto avgPositionWriter = std::make_unique<ParticleSimulation::AverageChargePositionWriter>(
                projectName+"_averagePosition.txt");

        auto fftWriter = std::make_unique<ParticleSimulation::InductionCurrentWriter>(
                particlePtrs, projectName+"_fft.txt", detectionPAs, detectionPAFactors, paSpatialScale);
        auto ionsInactiveWriter = std::make_unique<ParticleSimulation::Scalar_writer>(projectName+"_ionsInactive.txt");

        ParticleSimulation::partAttribTransformFctType particleAttributesTransformFct =
                [](BTree::Particle* particle) -> std::vector<double> {
                    std::vector<double> result = {
                            particle->getVelocity().x(),
                            particle->getVelocity().y(),
                            particle->getVelocity().z(),
                            particle->getFloatAttribute(key_trapForce_x),
                            particle->getFloatAttribute(key_trapForce_y),
                            particle->getFloatAttribute(key_trapForce_z),
                            particle->getFloatAttribute(key_spaceCharge_x),
                            particle->getFloatAttribute(key_spaceCharge_y),
                            particle->getFloatAttribute(key_spaceCharge_z),
                            particle->getMass()/Core::AMU_TO_KG,
                            particle->getCharge()
                    };
                    return result;
                };

        std::vector<std::string> particleAttributesNames = {"velocity x", "velocity y", "velocity z",
                                                            "rf x", "rf y", "rf z",
                                                            "spacecharge x", "spacecharge y", "spacecharge z",
                                                            "mass", "charge"};

        ParticleSimulation::partAttribTransformFctTypeInteger integerParticleAttributesTransformFct =
                [](BTree::Particle* particle) -> std::vector<int> {
                    std::vector<int> result = {
                            particle->getIntegerAttribute("global index"),
                            particle->getIndex()
                    };
                    return result;
                };

        std::vector<std::string> integerParticleAttributesNames = {"global index", "index"};

        auto hdf5Writer = std::make_unique<ParticleSimulation::TrajectoryHDF5Writer>(projectName+"_trajectories.hd5");
        hdf5Writer->setParticleAttributes(particleAttributesNames, particleAttributesTransformFct);
        hdf5Writer->setParticleAttributes(integerParticleAttributesNames, integerParticleAttributesTransformFct);

        ParticleSimulation::AbstractTimeIntegrator* integratorPtr;
        auto timestepWriteFunction =
                [trajectoryWriteInterval, fftWriteInterval, fftWriteMode, &V_0, &V_rf_export, &ionsInactive,
                        &hdf5Writer, &startSplatTracker, &ionsInactiveWriter, &fftWriter, &integratorPtr, &logger](
                        std::vector<BTree::Particle*>& particles, auto& /*tree*/, double time, int timestep,
                        bool lastTimestep) {

                    // check if simulation should be terminated (if all particles are terminated)
                    if (ionsInactive>=particles.size() && particles.size()>0) {
                        integratorPtr->setTerminationState();
                    }

                    // process time step data and write / export results
                    if (timestep%fftWriteInterval==0) {
                        ionsInactiveWriter->writeTimestep(ionsInactive, time);
                        if (fftWriteMode==UNRESOLVED) {
                            fftWriter->writeTimestep(time);
                        }
                        else if (fftWriteMode==MASS_RESOLVED) {
                            //fftWriter->writeTimestepMassResolved(time);
                            std::stringstream ss;
                            ss << "Mass resolved induction current not implemented";
                            throw (std::runtime_error(ss.str()));
                        }
                    }

                    if (timestep%trajectoryWriteInterval==0 || lastTimestep) {
                        logger->info("ts:{} time:{:.2e} V_rf:{:.1f} ions existing:{} ions inactive:{}",
                                timestep, time, V_0, particles.size(), ionsInactive);
                        V_rf_export.emplace_back(V_0);
                        hdf5Writer->writeTimestep(particles, time);
                    }

                    if (lastTimestep) {
                        hdf5Writer->writeStartSplatData(startSplatTracker);
                        hdf5Writer->finalizeTrajectory();
                        logger->info("finished ts:{} time:{:.2e}", timestep, time);
                    }
                };

        CollisionModel::HardSphereModel hsModel(backgroundPressure,
                backgroundTemperature,
                collisionGasMassAmu,
                collisionGasDiameterM);

        // simulate ===============================================================================================
        AppUtils::Stopwatch stopWatch;
        stopWatch.start();

        if (integratorMode==PARALLEL_VERLET) {
            ParticleSimulation::ParallelVerletIntegrator verletIntegrator(
                    particlePtrs,
                    accelerationFunctionLIT_parallel, timestepWriteFunction, otherActionsFunctionQIT,
                    particleStartMonitoringFct,
                    &hsModel);

            integratorPtr = &verletIntegrator;

            AppUtils::SignalHandler::setReceiver(verletIntegrator);
            verletIntegrator.run(timeSteps, dt);
        }

        if (rfMode==RAMPED_RF) {
            hdf5Writer->writeNumericListDataset("V_rf", V_rf_export);
        }
        stopWatch.stop();

        logger->info("CPU time: {} s", stopWatch.elapsedSecondsCPU());
        logger->info("Finished in {} seconds (wall clock time)", stopWatch.elapsedSecondsWall());
        return 0;
    }
    catch(const ParticleSimulation::PotentialArrayException& pe)
    {
        std::cout << pe.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch(const std::invalid_argument& ia){
        std::cout << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
}
