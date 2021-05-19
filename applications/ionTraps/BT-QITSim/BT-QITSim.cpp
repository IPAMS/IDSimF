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
 BT-QITSim.cpp

 Ion trajectory simulation (including space charge and hard sphere collisions) in an idealized
 quadrupole ion trap (QIT)

 ****************************/

#include "BTree_particle.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_util.hpp"
#include "PSim_verletIntegrator.hpp"
#include "PSim_parallelVerletIntegrator.hpp"
#include "PSim_scalar_writer.hpp"
#include "PSim_sampledWaveform.hpp"
#include "PSim_math.hpp"
#include "PSim_averageChargePositionWriter.hpp"
#include "PSim_idealizedQitFFTWriter.hpp"
#include "PSim_particleStartSplatTracker.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_ionDefinitionReading.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include <iostream>
#include <vector>

enum IntegratorMode {VERLET,PARALLEL_VERLET};
enum GeometryMode {DEFAULT,SCALED,VARIABLE};
enum RfAmplitudeMode {STATIC_RF,RAMPED_RF};
enum RfWaveMode {SINE,SAMPLED};
enum FieldMode {BASIC, HIGHER_ORDERS};
enum ExciteMode {RECTPULSE,SWIFT};
enum FftWriteMode {UNRESOLVED,MASS_RESOLVED};

std::string key_spaceCharge_x = "keySpaceChargeX";
std::string key_spaceCharge_y = "keySpaceChargeY";
std::string key_spaceCharge_z = "keySpaceChargeZ";
std::string key_rf_x = "keyRFForceX";
std::string key_rf_y = "keyRFForceY";
std::string key_rf_z = "keyRFForceZ";
std::string key_velocity_x = "keyVeloX";
std::string key_velocity_y = "keyVeloY";
std::string key_velocity_z = "keyVeloZ";

//some constants:
constexpr double r_0_default = 10.0 / 1000.0;
constexpr double z_0_default = 7.0  / 1000.0;


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
        std::string confBasePath = simConf.confBasePath();

        // read basic simulation parameters =============================================================
        std::string integratorMode_str = simConf.stringParameter("integrator_mode");
        IntegratorMode integratorMode;
        if (integratorMode_str=="verlet") {
            integratorMode = VERLET;
        }
        else if (integratorMode_str=="parallel_verlet") {
            integratorMode = PARALLEL_VERLET;
        }
        else {
            throw std::invalid_argument("wrong configuration value: integrator mode");
        }

        int timeSteps = simConf.intParameter("sim_time_steps");
        int trajectoryWriteInterval = simConf.intParameter("trajectory_write_interval");
        int fftWriteInterval = simConf.intParameter("fft_write_interval");
        double dt = simConf.doubleParameter("dt");
        std::vector<int> nIons = std::vector<int>();
        std::vector<double> ionMasses = std::vector<double>();
        std::vector<double> ionCollisionDiameters_angstrom = std::vector<double>();
        std::string fftWriteMode_str = simConf.stringParameter("fft_write_mode");
        FftWriteMode fftWriteMode;
        if (fftWriteMode_str=="unresolved") {
            fftWriteMode = UNRESOLVED;
        }
        else if (fftWriteMode_str=="mass_resolved") {
            fftWriteMode = MASS_RESOLVED;
        }

        //read geometrical configuration of the trap======================================================
        double r_0 = 0.0;
        double z_0 = 0.0;

        std::string geometryMode_str = simConf.stringParameter("geometry_mode");
        GeometryMode geometryMode;
        if (geometryMode_str=="default") {
            geometryMode = DEFAULT;
            r_0 = r_0_default;
            z_0 = z_0_default;
        }
        else if (geometryMode_str=="variable") {
            geometryMode = VARIABLE;
            r_0 = simConf.doubleParameter("r_0");
            z_0 = simConf.doubleParameter("z_0");
        }
        else if (geometryMode_str=="scaled") {
            geometryMode = SCALED;
            double geomScale = simConf.doubleParameter("geometry_scale");
            r_0 = r_0_default*geomScale;
            z_0 = z_0_default*geomScale;
        }
        else {
            throw std::invalid_argument("wrong configuration value: geometry_mode");
        }

        double d_square = r_0*r_0+2*z_0*z_0;
        double d_square_2 = d_square/2.0;
        double U_0 = 0.0;


        //read physical configuration ===================================================================
        double maxIonRadius = simConf.doubleParameter("max_ion_radius");
        double backgroundPressure = simConf.doubleParameter("background_pressure_Pa");
        double backgroundTemperature = simConf.doubleParameter("background_temperature_K");
        double spaceChargeFactor = simConf.doubleParameter("space_charge_factor");
        double collisionGasMassAmu = simConf.doubleParameter("collision_gas_mass_amu");
        double collisionGasDiameterM = simConf.doubleParameter("collision_gas_diameter_angstrom")*1e-10;


        //read rf configuration =========================================================================
        std::string fieldMode_str = simConf.stringParameter("field_mode");
        FieldMode fieldMode;
        std::vector<double> higherFieldOrdersCoefficients;
        if (fieldMode_str=="basic") {
            fieldMode = BASIC;
        }
        else if (fieldMode_str=="higher_orders") {
            fieldMode = HIGHER_ORDERS;
            higherFieldOrdersCoefficients = simConf.doubleVectorParameter("field_higher_orders_coeffs");
        }
        else {
            throw std::invalid_argument("wrong configuration value: field_mode");
        }

        double f_rf = simConf.doubleParameter("f_rf"); //RF frequency 1e6;
        double omega = f_rf*2.0*M_PI; //RF angular frequencyf_rf* 2.0 * M_PI;

        // read sampled RF waveform
        RfWaveMode rfWaveMode;
        std::unique_ptr<ParticleSimulation::SampledWaveform> rfSampledWaveForm;
        if (simConf.isParameter("rf_waveform_csv_file")) {
            rfWaveMode = SAMPLED;
            std::string rfWaveformFileName = simConf.pathRelativeToConfFile(
                    simConf.stringParameter("rf_waveform_csv_file"));

            rfSampledWaveForm = std::make_unique<ParticleSimulation::SampledWaveform>(rfWaveformFileName);
            if (!rfSampledWaveForm->good()) {
                logger->error("rf waveform file not accessible");
                return EXIT_FAILURE;
            }
        }
        else {
            rfWaveMode = SINE;
        }

        // read RF amplitude configuration
        RfAmplitudeMode rfAmplitudeMode;
        std::vector<double> V_0_ramp;
        double V_0 = 0.0;
        if (simConf.isParameter("V_rf_start")) {
            rfAmplitudeMode = RAMPED_RF;
            double V_rf_start = simConf.doubleParameter("V_rf_start");
            double V_rf_end = simConf.doubleParameter("V_rf_end");
            V_0_ramp = ParticleSimulation::linspace(V_rf_start, V_rf_end, timeSteps);
        }
        else {
            rfAmplitudeMode = STATIC_RF;
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
                return EXIT_FAILURE;
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
            particle->setFloatAttribute(key_rf_x, 0.0);
            particle->setFloatAttribute(key_rf_y, 0.0);
            particle->setFloatAttribute(key_rf_z, 0.0);
            particle->setFloatAttribute(key_spaceCharge_x, 0.0);
            particle->setFloatAttribute(key_spaceCharge_y, 0.0);
            particle->setFloatAttribute(key_spaceCharge_z, 0.0);
        }

        // define functions for the trajectory integration ==================================================
        auto trapFieldFunction =
                [exciteMode, rfWaveMode, rfAmplitudeMode, fieldMode, excitePulseLength, excitePulsePotential,
                        spaceChargeFactor, omega, z_0, U_0, d_square_2,
                        &rfSampledWaveForm, &higherFieldOrdersCoefficients, &swiftWaveForm, &V_0, &V_0_ramp](
                        BTree::Particle* particle, int particleIndex,
                        auto& tree, double time, int timestep) -> Core::Vector {

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

                    if (rfAmplitudeMode==RAMPED_RF) {
                        V_0 = V_0_ramp[timestep];
                    }

                    double V_rf;
                    if (rfWaveMode==SINE) {
                        V_rf = V_0*cos(omega*time);
                    }
                    else if (rfWaveMode==SAMPLED) {
                        V_rf = V_0*rfSampledWaveForm->getValueLooped(timestep);
                    }

                    double a_ex = excitePotential/z_0*particleCharge;

                    Core::Vector rfForce(0, 0, 0);
                    if (fieldMode==BASIC) {
                        double a = (U_0+V_rf)/d_square_2*particleCharge;

                        rfForce = Core::Vector(
                                a*pos.x(),
                                a*pos.y(),
                                -2*a*pos.z()+a_ex);
                    }
                    else if (fieldMode==HIGHER_ORDERS) {
                        // transform to z-r space, calculate higher orders and transform force back to cartesian space
                        double r = std::sqrt(pos.x()*pos.x()+pos.y()*pos.y());
                        double z = pos.z();
                        double phi = std::atan2(pos.x(), pos.y());

                        //d_square_2 == r_0^2 for ideal electrode geometry
                        double phi_0 = (U_0+V_rf)/2.0*particleCharge;

                        // ideal quadropole field:
                        double E_2_r = 2.0*r/d_square_2;
                        double E_2_z = -4.0*z/d_square_2;

                        //hexapole field:
                        double r_0_3 = std::pow(d_square_2, 3.0/2.0);
                        double E_3_r = -6*r*z/r_0_3;
                        double E_3_z = -3*(r*r-2*z*z)/r_0_3;

                        //octapole field:
                        double r_0_4 = d_square_2*d_square_2;
                        double E_4_r = 12*(r*r*r-4*r*z*z)/r_0_4;
                        double E_4_z = (32*z*z*z-48*r*r*z)/r_0_4;

                        double E_r = phi_0*(E_2_r
                                +higherFieldOrdersCoefficients[0]*E_3_r
                                +higherFieldOrdersCoefficients[1]*E_4_r);

                        double E_z = phi_0*(E_2_z
                                +higherFieldOrdersCoefficients[0]*E_3_z
                                +higherFieldOrdersCoefficients[1]*E_4_z);

                        // transform back to cartesian space:
                        rfForce = Core::Vector(
                                E_r*std::sin(phi),
                                E_r*std::cos(phi),
                                E_z+a_ex);
                    }
                    return rfForce;
                };

        auto accelerationFunctionQIT =
                [exciteMode, rfAmplitudeMode, fieldMode, excitePulseLength, excitePulsePotential,
                        spaceChargeFactor, omega, z_0, U_0, d_square_2,
                        &higherFieldOrdersCoefficients, &swiftWaveForm, &V_0, &V_0_ramp, &trapFieldFunction](
                        BTree::Particle* particle, int particleIndex,
                        auto& tree, double time, int timestep) -> Core::Vector {

                    Core::Vector pos = particle->getLocation();
                    double particleCharge = particle->getCharge();

                    Core::Vector rfForce = trapFieldFunction(particle, particleIndex, tree, time, timestep);

                    Core::Vector spaceChargeForce(0, 0, 0);
                    if (spaceChargeFactor>0) {
                        spaceChargeForce =
                                tree.computeEFieldFromTree(*particle)*(particleCharge*spaceChargeFactor);
                    }

                    //update the additional parameters for writing them later to the trajectory:
                    particle->setFloatAttribute(key_rf_x, rfForce.x());
                    particle->setFloatAttribute(key_rf_y, rfForce.y());
                    particle->setFloatAttribute(key_rf_z, rfForce.z());
                    particle->setFloatAttribute(key_spaceCharge_x, spaceChargeForce.x());
                    particle->setFloatAttribute(key_spaceCharge_y, spaceChargeForce.y());
                    particle->setFloatAttribute(key_spaceCharge_z, spaceChargeForce.z());

                    return ((rfForce+spaceChargeForce)/particle->getMass());
                };

        auto accelerationFunctionQIT_parallel =
                [exciteMode, rfAmplitudeMode, fieldMode, excitePulseLength, excitePulsePotential,
                        spaceChargeFactor, omega, z_0, U_0, d_square_2,
                        &higherFieldOrdersCoefficients, &swiftWaveForm, &V_0, &V_0_ramp, &trapFieldFunction](
                        BTree::Particle* particle, int particleIndex,
                        auto& tree, double time, int timestep) -> Core::Vector {

                    Core::Vector pos = particle->getLocation();
                    double particleCharge = particle->getCharge();

                    Core::Vector rfForce = trapFieldFunction(particle, particleIndex, tree, time, timestep);

                    Core::Vector spaceChargeForce(0, 0, 0);
                    if (spaceChargeFactor>0) {
                        spaceChargeForce =
                                tree.computeEFieldFromTree(*particle)*(particleCharge*spaceChargeFactor);
                    }


                    //update the additional parameters for writing them later to the trajectory:
                    particle->setFloatAttribute(key_rf_x, rfForce.x());
                    particle->setFloatAttribute(key_rf_y, rfForce.y());
                    particle->setFloatAttribute(key_rf_z, rfForce.z());
                    particle->setFloatAttribute(key_spaceCharge_x, spaceChargeForce.x());
                    particle->setFloatAttribute(key_spaceCharge_y, spaceChargeForce.y());
                    particle->setFloatAttribute(key_spaceCharge_z, spaceChargeForce.z());

                    return ((rfForce+spaceChargeForce)/particle->getMass());
                };

        // Prepare ion start / stop tracker and ion start monitoring / ion termination functions
        ParticleSimulation::ParticleStartSplatTracker startSplatTracker;
        auto particleStartMonitoringFct = [&startSplatTracker](BTree::Particle* particle, double time) {
            startSplatTracker.particleStart(particle, time);
        };

        int ionsInactive = 0;
        auto otherActionsFunctionQIT = [maxIonRadius, &ionsInactive, &startSplatTracker]
                (Core::Vector& newPartPos, BTree::Particle* particle, int particleIndex,
                 auto& tree, double time, int timestep) {
            if (newPartPos.magnitude()>maxIonRadius) {
                particle->setActive(false);
                particle->setSplatTime(time);
                startSplatTracker.particleSplat(particle, time);
                ionsInactive++;
            }
        };

        //prepare file writers and data writing functions ==============================================================================
        auto avgPositionWriter = std::make_unique<ParticleSimulation::AverageChargePositionWriter>(
                projectName+"_averagePosition.txt");
        auto fftWriter = std::make_unique<ParticleSimulation::IdealizedQitFFTWriter>(particlePtrs,
                projectName+"_fft.txt");

        auto ionsInactiveWriter = std::make_unique<ParticleSimulation::Scalar_writer>(projectName+"_ionsInactive.txt");
        ParticleSimulation::partAttribTransformFctType additionalParameterTransformFct =
                [](BTree::Particle* particle) -> std::vector<double> {
                    std::vector<double> result = {
                            particle->getVelocity().x(),
                            particle->getVelocity().y(),
                            particle->getVelocity().z(),
                            particle->getFloatAttribute(key_rf_x),
                            particle->getFloatAttribute(key_rf_y),
                            particle->getFloatAttribute(key_rf_z),
                            particle->getFloatAttribute(key_spaceCharge_x),
                            particle->getFloatAttribute(key_spaceCharge_y),
                            particle->getFloatAttribute(key_spaceCharge_z),
                    };
                    return result;
                };

        std::vector<std::string> auxParamNames = {"velocity x", "velocity y", "velocity z",
                                                  "rf x", "rf y", "rf z",
                                                  "spacecharge x", "spacecharge y", "spacecharge z"};

        ParticleSimulation::partAttribTransformFctTypeInteger integerParticleAttributesTransformFct =
                [](BTree::Particle* particle) -> std::vector<int> {
                    std::vector<int> result = {
                            particle->getIntegerAttribute("global index")
                    };
                    return result;
                };

        std::vector<std::string> integerParticleAttributesNames = {"global index"};

        auto hdf5Writer = std::make_unique<ParticleSimulation::TrajectoryHDF5Writer>(projectName+"_trajectories.hd5");
        hdf5Writer->setParticleAttributes(auxParamNames, additionalParameterTransformFct);
        hdf5Writer->setParticleAttributes(integerParticleAttributesNames, integerParticleAttributesTransformFct);

        auto timestepWriteFunction =
                [trajectoryWriteInterval, fftWriteInterval, fftWriteMode, &V_0, &V_rf_export, &ionsInactive,
                        &hdf5Writer, &ionsInactiveWriter,
                        &fftWriter, &startSplatTracker, &logger](
                        std::vector<BTree::Particle*>& particles, auto& tree, double time, int timestep,
                        bool lastTimestep) {

                    if (timestep%fftWriteInterval==0) {
                        ionsInactiveWriter->writeTimestep(ionsInactive, time);
                        if (fftWriteMode==UNRESOLVED) {
                            fftWriter->writeTimestep(time);
                        }
                        else if (fftWriteMode==MASS_RESOLVED) {
                            fftWriter->writeTimestepMassResolved(time);
                        }
                    }

                    if (lastTimestep) {
                        V_rf_export.emplace_back(V_0);
                        hdf5Writer->writeStartSplatData(startSplatTracker);
                        hdf5Writer->writeTimestep(particles, time);
                        hdf5Writer->finalizeTrajectory();
                        logger->info("finished ts:{} time:{:.2e}", timestep, time);
                        std::cout << "finished ts:" << timestep << " time:" << time << std::endl;
                    }
                    else if (timestep%trajectoryWriteInterval==0) {
                        logger->info("ts:{} time:{:.2e} V_rf:{:.1f} ions existing:{} ions inactive:{}",
                                timestep, time, V_0, particles.size(), ionsInactive);
                        V_rf_export.emplace_back(V_0);
                        hdf5Writer->writeTimestep(particles, time);
                    }
                };

        CollisionModel::HardSphereModel hsModel(backgroundPressure,
                backgroundTemperature,
                collisionGasMassAmu,
                collisionGasDiameterM);


        // simulate ===============================================================================================
        AppUtils::Stopwatch stopWatch;
        stopWatch.start();

        if (integratorMode==VERLET) {
            ParticleSimulation::VerletIntegrator verletIntegrator(
                    particlePtrs,
                    accelerationFunctionQIT, timestepWriteFunction,
                    otherActionsFunctionQIT, particleStartMonitoringFct,
                    &hsModel);
            AppUtils::SignalHandler::setReceiver(verletIntegrator);
            verletIntegrator.run(timeSteps, dt);
        }
        else if (integratorMode==PARALLEL_VERLET) {
            ParticleSimulation::ParallelVerletIntegrator verletIntegrator(
                    particlePtrs,
                    accelerationFunctionQIT_parallel, timestepWriteFunction,
                    otherActionsFunctionQIT, particleStartMonitoringFct,
                    &hsModel);
            AppUtils::SignalHandler::setReceiver(verletIntegrator);
            verletIntegrator.run(timeSteps, dt);
        }

        if (rfAmplitudeMode==RAMPED_RF) {
            hdf5Writer->writeNumericListDataset("V_rf", V_rf_export);
        }

        stopWatch.stop();

        logger->info("CPU time: {} s", stopWatch.elapsedSecondsCPU());
        logger->info("Finished in {} seconds (wall clock time)", stopWatch.elapsedSecondsWall());
        return EXIT_SUCCESS;
    }
    catch(const std::invalid_argument& ia){
        std::cout << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
}
