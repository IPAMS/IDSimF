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
 QITSim.cpp

 Ion trajectory simulation (including space charge and hard sphere collisions) in an idealized
 quadrupole ion trap (QIT)

 ****************************/

#include "Core_particle.hpp"
#include "FileIO_trajectoryHDF5Writer.hpp"
#include "PSim_util.hpp"
#include "Integration_verletIntegrator.hpp"
#include "Integration_parallelVerletIntegrator.hpp"
#include "Integration_parallelRK4Integrator.hpp"
#ifdef WITH_FMM_3d
#include "Integration_fmmIntegrator.hpp"
#include "FMM3D_fmmSolver.hpp"
#endif

#ifdef WITH_EXAFMMT
#include "Integration_fmmIntegrator.hpp"
#include "ExaFMMt_fmmSolver.hpp"
#endif

#include "FileIO_scalar_writer.hpp"
#include "PSim_sampledWaveform.hpp"
#include "PSim_math.hpp"
#include "FileIO_averageChargePositionWriter.hpp"
#include "FileIO_idealizedQitFFTWriter.hpp"
#include "PSim_particleStartSplatTracker.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_ionDefinitionReading.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "appUtils_commandlineParser.hpp"
#include <iostream>
#include <vector>

enum IntegratorMode {VERLET, PARALLEL_VERLET, FMM3D_VERLET, EXAFMM_VERLET, PARALLEL_RUNGE_KUTTA4};
enum GeometryMode {DEFAULT, SCALED,VARIABLE};
enum RfAmplitudeMode {STATIC_RF, RAMPED_RF};
enum RfWaveMode {SINE, SAMPLED};
enum FieldMode {BASIC, HIGHER_ORDERS, AVERAGED};
enum ExciteMode {RECTPULSE, SWIFT};
enum FftWriteMode {UNRESOLVED, MASS_RESOLVED, ION_CLOUD_POSITION};

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
        // parse commandline / create conf and logger ===================================================
        AppUtils::CommandlineParser cmdLineParser(argc, argv, "BT-QITSim",
                "Idealized Quadrupole Ion Trap (QIT) simulation", true);
        std::string simResultBasename = cmdLineParser.resultName();
        AppUtils::logger_ptr logger = cmdLineParser.logger();
        AppUtils::simConf_ptr simConf = cmdLineParser.simulationConfiguration();
        std::filesystem::path confBasePath = simConf->confBasePath();

        // read basic simulation parameters =============================================================
        std::string integratorMode_str = simConf->stringParameter("integrator_mode");
        IntegratorMode integratorMode;
        if (integratorMode_str=="verlet") {
            integratorMode = VERLET;
        }
        else if (integratorMode_str=="parallel_verlet") {
            integratorMode = PARALLEL_VERLET;
        }
        else if (integratorMode_str=="RK4") {
            integratorMode = PARALLEL_RUNGE_KUTTA4;
        }
#ifdef WITH_FMM_3d
        else if (integratorMode_str=="FMM3D_verlet") {
            integratorMode = FMM3D_VERLET;
        }
#endif
#ifdef WITH_EXAFMMT
        else if (integratorMode_str=="ExaFMM_verlet") {
            integratorMode = EXAFMM_VERLET;
        }
#endif
        else {
            throw std::invalid_argument("wrong configuration value: integrator mode");
        }

        unsigned int timeSteps = simConf->unsignedIntParameter("sim_time_steps");
        unsigned int trajectoryWriteInterval = simConf->unsignedIntParameter("trajectory_write_interval");
        unsigned int fftWriteInterval = simConf->unsignedIntParameter("fft_write_interval");
        double dt = simConf->doubleParameter("dt");
        std::vector<int> nIons = std::vector<int>();
        std::vector<double> ionMasses = std::vector<double>();
        std::vector<double> ionCollisionDiameters_angstrom = std::vector<double>();
        std::string fftWriteMode_str = simConf->stringParameter("fft_write_mode");
        FftWriteMode fftWriteMode = UNRESOLVED;
        if (fftWriteMode_str=="unresolved") {
            fftWriteMode = UNRESOLVED;
        }
        else if (fftWriteMode_str=="mass_resolved") {
            fftWriteMode = MASS_RESOLVED;
        }
        else if (fftWriteMode_str=="ion_cloud_position") {
            fftWriteMode = ION_CLOUD_POSITION;
        }

        //read geometrical configuration of the trap======================================================
        double r_0 = 0.0;
        double z_0 = 0.0;

        std::string geometryMode_str = simConf->stringParameter("geometry_mode");
        //GeometryMode geometryMode;
        if (geometryMode_str=="default") {
            //geometryMode = DEFAULT;
            r_0 = r_0_default;
            z_0 = z_0_default;
        }
        else if (geometryMode_str=="variable") {
            //geometryMode = VARIABLE;
            r_0 = simConf->doubleParameter("r_0");
            z_0 = simConf->doubleParameter("z_0");
        }
        else if (geometryMode_str=="scaled") {
            //geometryMode = SCALED;
            double geomScale = simConf->doubleParameter("geometry_scale");
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
        double maxIonRadius = simConf->doubleParameter("max_ion_radius");
        double backgroundPressure = simConf->doubleParameter("background_gas_pressure_Pa");
        double backgroundTemperature = simConf->doubleParameter("background_gas_temperature_K");
        double spaceChargeFactor = simConf->doubleParameter("space_charge_factor");
        double collisionGasMassAmu = simConf->doubleParameter("collision_gas_mass_amu");
        double collisionGasDiameterM = simConf->doubleParameter("collision_gas_diameter_angstrom")*1e-10;


        //read rf configuration =========================================================================
        std::string fieldMode_str = simConf->stringParameter("field_mode");
        FieldMode fieldMode;
        std::vector<double> higherFieldOrdersCoefficients;
        if (fieldMode_str=="basic") {
            fieldMode = BASIC;
        }
        else if (fieldMode_str=="higher_orders") {
            fieldMode = HIGHER_ORDERS;
            higherFieldOrdersCoefficients = simConf->doubleVectorParameter("field_higher_orders_coeffs");
        }
        else if (fieldMode_str=="averaged") {
            fieldMode = AVERAGED;
        }
        else {
            throw std::invalid_argument("wrong configuration value: field_mode");
        }

        double f_rf = simConf->doubleParameter("f_rf"); //RF frequency 1e6;
        double omega = f_rf*2.0*M_PI; //RF angular frequencyf_rf* 2.0 * M_PI;

        // read sampled RF waveform
        RfWaveMode rfWaveMode;
        std::unique_ptr<ParticleSimulation::SampledWaveform> rfSampledWaveForm;
        if (simConf->isParameter("rf_waveform_csv_file")) {
            rfWaveMode = SAMPLED;
            std::string rfWaveformFileName = simConf->pathRelativeToConfFile(
                    simConf->stringParameter("rf_waveform_csv_file"));

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
        if (simConf->isParameter("V_rf_start")) {
            rfAmplitudeMode = RAMPED_RF;
            double V_rf_start = simConf->doubleParameter("V_rf_start");
            double V_rf_end = simConf->doubleParameter("V_rf_end");
            V_0_ramp = ParticleSimulation::linspace(V_rf_start, V_rf_end, static_cast<int>(timeSteps));
        }
        else {
            rfAmplitudeMode = STATIC_RF;
            V_0 = simConf->doubleParameter("V_rf");
        }
        std::vector<double> V_rf_export;


        //read excitation / swift configuration ========================================================

        ExciteMode exciteMode;
        std::unique_ptr<ParticleSimulation::SampledWaveform> swiftWaveForm;
        double excitePulseLength = 0.0;
        if (simConf->isParameter("excite_waveform_csv_file")) {
            exciteMode = SWIFT;
            std::string swiftFileName = simConf->stringParameter("excite_waveform_csv_file");
            swiftWaveForm = std::make_unique<ParticleSimulation::SampledWaveform>(swiftFileName);
            if (!swiftWaveForm->good()) {
                logger->error("swift transient file not accessible");
                return EXIT_FAILURE;
            }
        }
        else {
            exciteMode = RECTPULSE;
            excitePulseLength = simConf->doubleParameter("excite_pulse_length");
        }
        double excitePulsePotential = simConf->doubleParameter("excite_pulse_potential");



        //read ion configuration =======================================================================
        std::vector<std::unique_ptr<Core::Particle>> particles;
        std::vector<Core::Particle*> particlePtrs;
        AppUtils::readIonDefinition(particles, particlePtrs, *simConf);
        // init additional ion parameters:
        for (const auto& particle: particles) {
            particle->setFloatAttribute(key_rf_x, 0.0);
            particle->setFloatAttribute(key_rf_y, 0.0);
            particle->setFloatAttribute(key_rf_z, 0.0);
            particle->setFloatAttribute(key_spaceCharge_x, 0.0);
            particle->setFloatAttribute(key_spaceCharge_y, 0.0);
            particle->setFloatAttribute(key_spaceCharge_z, 0.0);
            ionMasses.emplace_back(particle->getMass()/Core::AMU_TO_KG);
        }

        // define functions for the trajectory integration ==================================================
        auto trapFieldFunction =
                [exciteMode, rfWaveMode, rfAmplitudeMode, fieldMode, excitePulseLength, excitePulsePotential,
                        omega, z_0, r_0, U_0, d_square_2,
                        &rfSampledWaveForm, &higherFieldOrdersCoefficients, &swiftWaveForm, &V_0, &V_0_ramp](
                        Core::Particle* particle, Core::Vector pPos, double time, unsigned int timestep) -> Core::Vector
                {
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
                    double a_ex = excitePotential/z_0*particleCharge;

                    if (fieldMode==AVERAGED){
                        // transform to z-r space, calculate higher orders and transform force back to cartesian space
                        double r = std::sqrt(pPos.x()*pPos.x()+pPos.y()*pPos.y());
                        double z = pPos.z();
                        double phi = std::atan2(pPos.x(), pPos.y());

                        double ponderomotiveForceFactor = -(particleCharge*particleCharge) / (2.0 * particle->getMass() * omega * omega);

                        double rz_0Factor = r_0*r_0 + 2*z_0*z_0;
                        double a_2_avg = V_0 * V_0 / (2*rz_0Factor*rz_0Factor);

                        double F_r = 8 * a_2_avg * r;
                        double F_z = 32 * a_2_avg * z;

                        // transform back to cartesian space:
                        Core::Vector averagedRFForce{
                                ponderomotiveForceFactor*F_r*std::sin(phi),
                                ponderomotiveForceFactor*F_r*std::cos(phi),
                                ponderomotiveForceFactor*F_z+a_ex};
                        //std::cout <<"V0 "<< V_0 << " forceFactor " << forceFactor<< "  ffRad "<< fieldFactorRadial <<"  xField: "<<xField<< "  zField: "<<zField<<std::endl;

                        return averagedRFForce;
                    }
                    else {
                        double V_rf = 0.0;
                        if (rfWaveMode==SINE) {
                            V_rf = V_0*cos(omega*time);
                        }
                        else if (rfWaveMode==SAMPLED) {
                            V_rf = V_0*rfSampledWaveForm->getValueLooped(timestep);
                        }

                        Core::Vector rfForce(0, 0, 0);
                        if (fieldMode==BASIC) {
                            double a = (U_0+V_rf)/d_square_2*particleCharge;

                            rfForce = Core::Vector(
                                    a*pPos.x(),
                                    a*pPos.y(),
                                    -2*a*pPos.z()+a_ex);
                        }
                        else if (fieldMode==HIGHER_ORDERS) {
                            // Source for this derivation of the field components:
                            // https://doi.org/10.1016/1044-0305(93)80017-S
                            // Nonlinear Resonance Effects During Ion Storage in a Quadrupole Ion Trap, Eades et.al.

                            // transform to z-r space, calculate higher orders and transform force back to cartesian space
                            double r = std::sqrt(pPos.x()*pPos.x()+pPos.y()*pPos.y());
                            double z = pPos.z();
                            double phi = std::atan2(pPos.x(), pPos.y());

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
                    }
                };

        auto spaceChargeAccelerationFct_RKIntegration =
                     [spaceChargeFactor](
                             Core::Particle* particle, int /*particleIndex*/,
                             SpaceCharge::FieldCalculator& scFieldCalculator, double /*time*/, unsigned int /*timestep*/) -> Core::Vector {

                     Core::Vector spaceChargeForce(0, 0, 0);
                     if (spaceChargeFactor>0) {
                         spaceChargeForce =
                                 scFieldCalculator.getEFieldFromSpaceCharge(*particle)*(particle->getCharge()*spaceChargeFactor);
                     }

                     return (spaceChargeForce/particle->getMass());
                };

        auto accelerationFctTrapField_RKIntegration =
                [&trapFieldFunction](Core::Particle* particle, Core::Vector position, Core::Vector /*velocity*/, double time, unsigned int timestep){
                    return trapFieldFunction(particle, position, time, timestep)/particle->getMass();
                };


        auto accelerationFctQIT_verletIntegration =
                [spaceChargeFactor, &trapFieldFunction](
                        Core::Particle* particle, int /*particleIndex*/,
                        SpaceCharge::FieldCalculator& scFieldCalculator, double time, unsigned int timestep) -> Core::Vector {

                    double particleCharge = particle->getCharge();

                    Core::Vector rfForce = trapFieldFunction(particle, particle->getLocation(), time, timestep);

                    Core::Vector spaceChargeForce(0, 0, 0);
                    if (spaceChargeFactor>0) {
                        spaceChargeForce =
                                scFieldCalculator.getEFieldFromSpaceCharge(*particle)*(particleCharge*spaceChargeFactor);
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
        auto particleStartMonitoringFct = [&startSplatTracker](Core::Particle* particle, double time) {
            startSplatTracker.particleStart(particle, time);
        };

        int ionsInactive = 0;
        auto otherActionsFunctionQIT = [maxIonRadius, &ionsInactive, &startSplatTracker]
                (Core::Vector& newPartPos, Core::Particle* particle, unsigned int /*particleIndex*/,
                 double time, unsigned int /*timestep*/) {
            if (newPartPos.magnitude()>maxIonRadius) {
                particle->setActive(false);
                particle->setSplatTime(time);
                startSplatTracker.particleSplat(particle, time);
                ionsInactive++;
            }
        };

        //prepare file writers and data writing functions ==============================================================================
        /*auto avgPositionWriter = std::make_unique<FileIO::AverageChargePositionWriter>(
                simResultBasename+"_averagePosition.txt");*/
        auto fftWriter = std::make_unique<FileIO::IdealizedQitFFTWriter>(particlePtrs,
                simResultBasename+"_fft.txt");

        auto ionsInactiveWriter = std::make_unique<FileIO::Scalar_writer>(simResultBasename+"_ionsInactive.txt");
        FileIO::partAttribTransformFctType additionalParameterTransformFct =
                [](Core::Particle* particle) -> std::vector<double> {
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

        FileIO::partAttribTransformFctTypeInteger integerParticleAttributesTransformFct =
                [](Core::Particle* particle) -> std::vector<int> {
                    std::vector<int> result = {
                            particle->getIntegerAttribute("global index")
                    };
                    return result;
                };

        std::vector<std::string> integerParticleAttributesNames = {"global index"};

        auto hdf5Writer = std::make_unique<FileIO::TrajectoryHDF5Writer>(simResultBasename+"_trajectories.hd5");
        hdf5Writer->setParticleAttributes(auxParamNames, additionalParameterTransformFct);
        hdf5Writer->setParticleAttributes(integerParticleAttributesNames, integerParticleAttributesTransformFct);

        auto timestepWriteFunction =
                [trajectoryWriteInterval, fftWriteInterval, fftWriteMode, &V_0, &V_rf_export, &ionsInactive,
                        &hdf5Writer, &ionsInactiveWriter,
                        &fftWriter, &startSplatTracker, &logger](
                        std::vector<Core::Particle*>& particles, double time, unsigned int timestep,
                        bool lastTimestep) {

                    if (timestep%fftWriteInterval==0) {
                        ionsInactiveWriter->writeTimestep(ionsInactive, time);
                        if (fftWriteMode==UNRESOLVED) {
                            fftWriter->writeTimestep(time);
                        }
                        else if (fftWriteMode==MASS_RESOLVED) {
                            fftWriter->writeTimestepMassResolved(time);
                        }
                        else if (fftWriteMode==ION_CLOUD_POSITION) {
                            fftWriter->writeTimestepAverageIonCloudPosition(time);
                        }
                    }

                    if (lastTimestep) {
                        V_rf_export.emplace_back(V_0);
                        hdf5Writer->writeStartSplatData(startSplatTracker);
                        hdf5Writer->writeTimestep(particles, time);
                        hdf5Writer->finalizeTrajectory();
                        logger->info("finished ts:{} time:{:.2e}", timestep, time);
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
            Integration::VerletIntegrator verletIntegrator(
                    particlePtrs,
                    accelerationFctQIT_verletIntegration, timestepWriteFunction,
                    otherActionsFunctionQIT, particleStartMonitoringFct,
                    &hsModel);
            AppUtils::SignalHandler::setReceiver(verletIntegrator);
            verletIntegrator.run(timeSteps, dt);
        }
        else if (integratorMode==PARALLEL_VERLET) {
            Integration::ParallelVerletIntegrator verletIntegrator(
                    particlePtrs,
                    accelerationFctQIT_verletIntegration, timestepWriteFunction,
                    otherActionsFunctionQIT, particleStartMonitoringFct,
                    &hsModel);
            AppUtils::SignalHandler::setReceiver(verletIntegrator);
            verletIntegrator.run(timeSteps, dt);
        }
        else if (integratorMode==PARALLEL_RUNGE_KUTTA4) {
            Integration::ParallelRK4Integrator rk4Integrator(
                    particlePtrs,
                    accelerationFctTrapField_RKIntegration, spaceChargeAccelerationFct_RKIntegration,
                    timestepWriteFunction, otherActionsFunctionQIT, particleStartMonitoringFct,
                    &hsModel);
            AppUtils::SignalHandler::setReceiver(rk4Integrator);
            rk4Integrator.run(timeSteps, dt);
        }
#ifdef WITH_FMM_3d
        else if (integratorMode==FMM3D_VERLET) {
            Integration::FMMVerletIntegrator<FMM3D::FMMSolver> integrator(
                    particlePtrs,
                    accelerationFctQIT_verletIntegration, timestepWriteFunction,
                    otherActionsFunctionQIT, particleStartMonitoringFct,
                    &hsModel);

            if (simConf->isParameter("FMM3D_precision")) {
                integrator.getFMMSolver()->setRequestedPrecision(simConf->doubleParameter("FMM3D_precision"));
            }

            AppUtils::SignalHandler::setReceiver(integrator);
            integrator.run(timeSteps, dt);
        }
#endif
#ifdef WITH_EXAFMMT
        else if (integratorMode==EXAFMM_VERLET) {
            Integration::FMMVerletIntegrator<ExaFMMt::FMMSolver> integrator(
                    particlePtrs,
                    accelerationFctQIT_verletIntegration, timestepWriteFunction,
                    otherActionsFunctionQIT, particleStartMonitoringFct,
                    &hsModel);

            if (simConf->isParameter("ExaFMM_order")) {
                integrator.getFMMSolver()->setExpansionOrder(simConf->intParameter("ExaFMM_order"));
            }

            AppUtils::SignalHandler::setReceiver(integrator);
            integrator.run(timeSteps, dt);
        }
#endif

        hdf5Writer->writeNumericListDataset("Particle Masses", ionMasses);

        if (rfAmplitudeMode==RAMPED_RF) {
            hdf5Writer->writeNumericListDataset("V_rf", V_rf_export);
        }
        stopWatch.stop();
        logger->info("CPU time: {} s", stopWatch.elapsedSecondsCPU());
        logger->info("Finished in {} seconds (wall clock time)", stopWatch.elapsedSecondsWall());
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
