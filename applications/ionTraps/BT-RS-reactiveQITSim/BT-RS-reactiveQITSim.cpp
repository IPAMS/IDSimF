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
 BT-RS-reactiveQITSim.cpp

 Ion trajectory simulation with reactive particles / ion chemistry
 (including space charge and hard sphere collisions) in an idealized
 quadrupole ion trap (QIT)

 ****************************/

#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_ConcentrationFileWriter.hpp"
#include "PSim_util.hpp"
#include "PSim_math.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_verletIntegrator.hpp"
#include "PSim_sampledWaveform.hpp"
#include "PSim_particleStartSplatTracker.hpp"
#include "PSim_scalar_writer.hpp"
#include "PSim_averageChargePositionWriter.hpp"
#include "PSim_idealizedQitFFTWriter.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "CollisionModel_MultiCollisionModel.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "appUtils_commandlineParser.hpp"
#include "json.h"
#include <iostream>
#include <vector>

enum GeometryMode {DEFAULT,SCALED,VARIABLE};
enum RfAmplitudeMode {STATIC_RF,RAMPED_RF};
enum ExciteMode {NOEXCITE,RECTPULSE,SWIFT,CONTINUOUSSINE};
enum FftWriteMode {OFF,UNRESOLVED,MASS_RESOLVED};

const std::string key_ChemicalIndex = "keyChemicalIndex";
const std::string key_Collisions_total = "keyBackgroundGasCollisionsTotal";
//some constants:
double r_0_default = 10.0 / 1000.0;
double z_0_default = 7.0  / 1000.0;

std::function<void(double,BTree::Particle&)> createCollisionCountFunction(std::string key){
    return [=](double /*collisionEnergy*/, BTree::Particle& ion)->void{
        int nCollisions = ion.getIntegerAttribute(key);
        ion.setIntegerAttribute(key, nCollisions+1);
    };
}

std::function<void(RS::CollisionConditions, BTree::Particle &)>
createCollisionReactionFunction(RS::Substance *collisionPartnerSubstance, RS::Simulation &rsSim, std::string key){

    return [collisionPartnerSubstance,key,&rsSim] (RS::CollisionConditions conditions, BTree::Particle& ion){
        rsSim.collisionReact(ion.getIndex(),collisionPartnerSubstance,conditions);
        int nCollisions = ion.getIntegerAttribute(key);
        ion.setIntegerAttribute(key, nCollisions+1);
    };
}


int main(int argc, const char * argv[]) {

    try {
        // parse commandline / create conf and logger ===================================================
        AppUtils::CommandlineParser cmdLineParser(argc, argv, "BT-RS-reactiveQITSim",
                "Quadrupole Ion Trap (QIT) simulation with chemical reactions", false);
        std::string simResultBasename = cmdLineParser.projectName();
        AppUtils::logger_ptr logger = cmdLineParser.logger();
        AppUtils::simConf_ptr simConf = cmdLineParser.simulationConfiguration();
        std::filesystem::path confBasePath = simConf->confBasePath();

        // read basic simulation parameters =============================================================
        int timeSteps = simConf->intParameter("sim_time_steps");
        int concentrationWriteInterval = simConf->intParameter("concentrations_write_interval");

        int trajectoryWriteInterval = simConf->intParameter("trajectory_write_interval");
        int fftWriteInterval = simConf->intParameter("fft_write_interval");
        double dt = simConf->doubleParameter("dt");
        std::string fftWriteMode_str = simConf->stringParameter("fft_mode");
        FftWriteMode fftWriteMode;
        if (fftWriteMode_str=="off") {
            fftWriteMode = OFF;
        }
        else if (fftWriteMode_str=="unresolved") {
            fftWriteMode = UNRESOLVED;
        }
        else if (fftWriteMode_str=="mass_resolved") {
            fftWriteMode = MASS_RESOLVED;
        }
        else {
            throw std::invalid_argument("wrong configuration value: fft_mode");
        }



        //read geometrical configuration of the trap======================================================
        double r_0 = 0.0;
        double z_0 = 0.0;

        std::string geometryMode_str = simConf->stringParameter("geometry_mode");
        if (geometryMode_str=="default") {
            r_0 = r_0_default;
            z_0 = z_0_default;
        }
        else if (geometryMode_str=="variable") {
            r_0 = simConf->doubleParameter("r_0");
            z_0 = simConf->doubleParameter("z_0");
        }
        else if (geometryMode_str=="scaled") {
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
        double maxIonRadius_m = simConf->doubleParameter("max_ion_radius_m");
        double spaceChargeFactor = simConf->doubleParameter("space_charge_factor");
        double backgroundTemperature_K = simConf->doubleParameter("background_temperature_K");

        //read collision gas / background gas configuration =============================================
        std::vector<std::string> collisionGasNames = simConf->stringVectorParameter("collision_gas_names");
        std::vector<double> partialPressures = simConf->doubleVectorParameter("partial_pressures_Pa");
        std::vector<double> collisionGasMasses_Amu = simConf->doubleVectorParameter("collision_gas_masses_amu");
        std::vector<double> collisionGasDiameters_m = simConf->doubleVectorParameter("collision_gas_diameters_angstrom");
        collisionGasDiameters_m[0] *= 1e-10;
        collisionGasDiameters_m[1] *= 1e-10;

        double totalBackgroundPressure_Pa = partialPressures[0]+partialPressures[1];

        double startWidth_m = simConf->doubleParameter("start_width_m");

        //read rf configuration =========================================================================
        double f_rf = simConf->doubleParameter("f_rf"); //RF frequency 1e6;
        double omega = f_rf*2.0*M_PI; //RF angular frequencyf_rf* 2.0 * M_PI;

        RfAmplitudeMode rfMode;
        std::vector<double> V_0_ramp;
        double V_0 = 0.0;
        if (simConf->isParameter("rf_ramp_start_V")) {
            rfMode = RAMPED_RF;
            int V_rf_waiting_ts = simConf->intParameter("rf_ramp_waiting_timesteps");
            double V_rf_start = simConf->doubleParameter("rf_ramp_start_V");
            double V_rf_stop = simConf->doubleParameter("rf_ramp_stop_V");

            auto vRf = ParticleSimulation::fillVector(V_rf_start, V_rf_waiting_ts);
            auto vRamp = ParticleSimulation::linspace(V_rf_start, V_rf_stop, timeSteps-V_rf_waiting_ts);
            vRf.insert(vRf.end(), vRamp.begin(), vRamp.end()); //concat vRf and vRamp
            V_0_ramp = vRf;
        }
        else {
            rfMode = STATIC_RF;
            V_0 = simConf->doubleParameter("rf_V");
        }


        //read excitation / swift configuration ========================================================
        ExciteMode exciteMode;
        std::unique_ptr<ParticleSimulation::SampledWaveform> swiftWaveForm;
        double excitePulseLength = 0.0;
        double exciteDivisor = 0.0;

        std::string exciteMode_str = simConf->stringParameter("excite_mode");
        if (exciteMode_str=="off") {
            exciteMode = NOEXCITE;
        }
        else if (exciteMode_str=="rect_pulse") {
            excitePulseLength = simConf->doubleParameter("excite_pulse_length");
            exciteMode = RECTPULSE;
        }
        else if (exciteMode_str=="waveform") {
            exciteMode = SWIFT;
            if (simConf->isParameter("excite_waveform_csv_file")) {
                exciteMode = SWIFT;
                std::string swiftFileName = simConf->stringParameter("excite_waveform_csv_file");
                swiftWaveForm = std::make_unique<ParticleSimulation::SampledWaveform>(swiftFileName);
                if (!swiftWaveForm->good()) {
                    logger->error("swift transient file not accessible");
                    return EXIT_FAILURE;
                }
            }
        }
        else if (exciteMode_str=="continuous_sine") {
            exciteMode = CONTINUOUSSINE;
            exciteDivisor = simConf->doubleParameter("excite_divisor");
        }
        else {
            throw std::invalid_argument("wrong configuration value: excite_mode");
        }
        double excitePotential = 0.0;
        if (exciteMode!=NOEXCITE) {
            excitePotential = simConf->doubleParameter("excite_potential");
        }

        //read and prepare chemical configuration ===============================================
        RS::ConfigFileParser parser = RS::ConfigFileParser();
        std::string rsFilePath = simConf->pathRelativeToConfFile(simConf->stringParameter("reaction_configuration"));
        RS::Simulation rsSim = RS::Simulation(parser.parseFile(rsFilePath));
        RS::SimulationConfiguration* rsSimConf = rsSim.simulationConfiguration();
        //prepare a map for retrieval of the substances from their index:
        std::map<RS::Substance*, int> substanceIndices;

        //prepare some parameter vectors required for the hdf5 trajectory writer:
        std::vector<std::string> discreteSubstanceNames;
        std::vector<double> discreteSubstanceMasses;
        std::vector<double> V_rf_export;

        std::vector<RS::Substance*> discreteSubstances = rsSimConf->getAllDiscreteSubstances();
        for (std::size_t i = 0; i<discreteSubstances.size(); i++) {
            substanceIndices.insert(std::pair<RS::Substance*, int>(discreteSubstances[i], i));
            discreteSubstanceNames.emplace_back(discreteSubstances[i]->name());
            discreteSubstanceMasses.emplace_back(discreteSubstances[i]->mass());
        }


        //read ion configuration =======================================================================
        unsigned int nParticlesTotal = 0;
        std::vector<uniqueReactivePartPtr> particles;
        std::vector<BTree::Particle*> particlePtrs;

        // read and init random ion box configuration
        std::vector<unsigned int> nIons = simConf->unsignedIntVectorParameter("n_ions");

        Core::Vector initCorner(-startWidth_m/2.0, -startWidth_m/2.0, -startWidth_m/2.0);
        Core::Vector initBoxSize(startWidth_m, startWidth_m, startWidth_m);

        for (std::size_t i = 0; i<nIons.size(); i++) {
            RS::Substance* subst = discreteSubstances[i];
            std::vector<Core::Vector> initialPositions =
                    ParticleSimulation::util::getRandomPositionsInBox(nIons[i], initCorner, initBoxSize);
            for (unsigned int k = 0; k<nIons[i]; k++) {
                uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(subst);

                particle->setIndex(nParticlesTotal);
                particle->setLocation(initialPositions[k]);
                particle->setIntegerAttribute(key_Collisions_total, 0);
                particlePtrs.push_back(particle.get());
                rsSim.addParticle(particle.get(), nParticlesTotal);
                particles.push_back(std::move(particle));
                nParticlesTotal++;
            }
        }

        RS::ReactionConditions reactionConditions = RS::ReactionConditions();
        reactionConditions.temperature = backgroundTemperature_K;
        reactionConditions.pressure = totalBackgroundPressure_Pa;

        //FIXME: set electric field for the individual reaction events?
        //reactionConditions.electricField = 0.0;
        //reactionConditions.kineticEnergy = 0.0;

        // define functions for the trajectory integration ==================================================

        // some variables for synchronization between calls of the acceleration function:
        unsigned int lastTimestep = 0;
        double parameter_a;
        double exciteCos;

        auto accelerationFunctionQIT =
                [exciteMode, rfMode, excitePulseLength, excitePotential,
                        spaceChargeFactor, omega, z_0, U_0, d_square_2,
                        &swiftWaveForm, exciteDivisor, &V_0, &V_0_ramp, &lastTimestep, &parameter_a, &exciteCos](
                        BTree::Particle* particle, int /*particleIndex*/,
                        BTree::Tree& tree, double time, unsigned int timestep) -> Core::Vector {

                    Core::Vector pos = particle->getLocation();
                    double particleCharge = particle->getCharge();
                    double tsExcitePotential = 0;
                    if (exciteMode==RECTPULSE) {
                        if (time<excitePulseLength) {
                            tsExcitePotential = excitePotential;
                        }
                    }
                    else if (exciteMode==SWIFT) {
                        tsExcitePotential = swiftWaveForm->getValue(timestep)*excitePotential;
                    }
                    else if (exciteMode==CONTINUOUSSINE) {
                        if ( (timestep+1) > lastTimestep) {
                            exciteCos = cos(omega/exciteDivisor*time)*excitePotential;
                        }
                        tsExcitePotential = exciteCos;
                    }

                    if (rfMode==RAMPED_RF) {
                        V_0 = V_0_ramp[timestep];
                    }

                    double a_ex = tsExcitePotential/z_0;

                    if ( (timestep+1) > lastTimestep) {
                        parameter_a = (U_0+(V_0*cos(omega*time)))/d_square_2;
                        lastTimestep = timestep;
                    }
                    Core::Vector rfForce = Core::Vector(
                            parameter_a*pos.x(),
                            parameter_a*pos.y(),
                            -2*parameter_a*pos.z()+a_ex)*particleCharge;

                    Core::Vector spaceChargeForce(0, 0, 0);
                    if (spaceChargeFactor>0) {
                        spaceChargeForce =
                                tree.computeEFieldFromTree(*particle)*(particleCharge*spaceChargeFactor);
                    }

                    //update the additional parameters for writing them later to the trajectory:
                    //particle->setAuxScalarParam(key_trapForce_x, rfForce.x());
                    return ((rfForce+spaceChargeForce)/particle->getMass());
                };

        // Prepare ion start / stop tracker and ion start monitoring / ion termination functions
        ParticleSimulation::ParticleStartSplatTracker startSplatTracker;
        auto particleStartMonitoringFct = [&startSplatTracker](BTree::Particle* particle, double time) {
            startSplatTracker.particleStart(particle, time);
        };

        unsigned int ionsInactive = 0;
        auto otherActionsFunctionQIT = [maxIonRadius_m, &ionsInactive, &startSplatTracker]
                (Core::Vector& newPartPos, BTree::Particle* particle,
                 int /*particleIndex*/, BTree::Tree& /*tree*/, double time, int /*timestep*/) {
            if (newPartPos.magnitude()>maxIonRadius_m) {

                particle->setActive(false);
                particle->setSplatTime(time);
                startSplatTracker.particleSplat(particle, time);
                ionsInactive++;
            }
        };

        //prepare file writers and data writing functions ==============================================================================
        auto avgPositionWriter = std::make_unique<ParticleSimulation::AverageChargePositionWriter>(
                simResultBasename+"_averagePosition.txt");
        std::unique_ptr<ParticleSimulation::IdealizedQitFFTWriter> fftWriter = nullptr;
        if (fftWriteMode!=OFF) {
            fftWriter = std::make_unique<ParticleSimulation::IdealizedQitFFTWriter>(particlePtrs,
                    simResultBasename+"_fft.txt");
        }
        auto ionsInactiveWriter = std::make_unique<ParticleSimulation::Scalar_writer>(simResultBasename+"_ionsInactive.txt");

        RS::ConcentrationFileWriter concentrationFilewriter(simResultBasename+"_concentrations.txt");
        concentrationFilewriter.initFile(rsSimConf);

        ParticleSimulation::partAttribTransformFctType additionalParameterTransformFct =
                [](BTree::Particle* particle) -> std::vector<double> {
                    double ionVelocity = particle->getVelocity().magnitude();
                    double kineticEnergy_eV = 0.5*particle->getMass()*ionVelocity*ionVelocity*Core::JOULE_TO_EV;
                    std::vector<double> result = {
                            particle->getVelocity().x(),
                            particle->getVelocity().y(),
                            particle->getVelocity().z(),
                            kineticEnergy_eV
                    };
                    return result;
                };

        std::vector<std::string> auxParamNames = {"velocity x", "velocity y", "velocity z", "kinetic energy (eV)"};

        ParticleSimulation::partAttribTransformFctTypeInteger integerParticleAttributesTransformFct =
                [](BTree::Particle* particle) -> std::vector<int> {
                    std::vector<int> result = {
                            particle->getIntegerAttribute("global index"),
                            particle->getIntegerAttribute(key_Collisions_total),
                            particle->getIntegerAttribute(key_ChemicalIndex)

                    };
                    return result;
                };

        std::vector<std::string> integerParticleAttributesNames = {"global index", "total collisions", "chemical id"};

        auto hdf5Writer = std::make_unique<ParticleSimulation::TrajectoryHDF5Writer>(
                simResultBasename+"_trajectories.hd5");
        hdf5Writer->setParticleAttributes(auxParamNames, additionalParameterTransformFct);
        hdf5Writer->setParticleAttributes(integerParticleAttributesNames, integerParticleAttributesTransformFct);

        auto timestepWriteFunction =
                [trajectoryWriteInterval, fftWriteInterval, fftWriteMode, &V_0, &V_rf_export, &ionsInactive,
                        &hdf5Writer, &avgPositionWriter, &ionsInactiveWriter, &fftWriter, &startSplatTracker, &logger](
                        std::vector<BTree::Particle*>& particles, BTree::Tree& tree, double time, int timestep,
                        bool lastTimestep) {

                    if (timestep%fftWriteInterval==0) {
                        avgPositionWriter->writeTimestep(tree, time);
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
                    }

                    else if (timestep%trajectoryWriteInterval==0) {
                        logger->info("ts:{} time:{:.2e} V_rf:{:.1f} ions existing:{} ions inactive:{}",
                                timestep, time, V_0, particles.size(), ionsInactive);
                        V_rf_export.emplace_back(V_0);
                        hdf5Writer->writeTimestep(particles, time);
                    }
                };

        //prepare background gas collision models and collision based chemical reactions
        std::vector<std::unique_ptr<CollisionModel::AbstractCollisionModel>> collisionModels;

        for (std::size_t i = 0; i<collisionGasNames.size(); ++i) {
            RS::Substance* collisionPartnerSubst = rsSimConf->substanceByName(collisionGasNames[i]);
            std::unique_ptr<CollisionModel::HardSphereModel> hsModel = std::make_unique<CollisionModel::HardSphereModel>(
                    partialPressures[i], backgroundTemperature_K,
                    collisionGasMasses_Amu[i], collisionGasDiameters_m[i],
                    createCollisionReactionFunction(collisionPartnerSubst, rsSim, key_Collisions_total));
            collisionModels.emplace_back(std::move(hsModel));
        }
        CollisionModel::MultiCollisionModel combinedCollisionModel(std::move(collisionModels));

        // simulate ===============================================================================================
        ParticleSimulation::VerletIntegrator verletIntegrator(
                particlePtrs,
                accelerationFunctionQIT, timestepWriteFunction, otherActionsFunctionQIT, particleStartMonitoringFct,
                &combinedCollisionModel);
        AppUtils::SignalHandler::setReceiver(verletIntegrator);

        AppUtils::Stopwatch stopWatch;
        stopWatch.start();

        for (int step = 0; step<timeSteps; ++step) {
            if (step%concentrationWriteInterval==0) {
                concentrationFilewriter.writeTimestep(rsSim);
            }
            if (step%trajectoryWriteInterval==0) {
                rsSim.logConcentrations(logger);
            }
            for (unsigned int i = 0; i<nParticlesTotal; i++) {
                if (particles[i]->isActive()) {
                    rsSim.react(i, reactionConditions, dt);
                    int substIndex = substanceIndices.at(particles[i]->getSpecies());
                    particles[i]->setIntegerAttribute(key_ChemicalIndex, substIndex);
                }
            }
            rsSim.advanceTimestep(dt);
            verletIntegrator.runSingleStep(dt);

            if (ionsInactive>=nParticlesTotal ||
                    verletIntegrator.runState()==ParticleSimulation::AbstractTimeIntegrator::IN_TERMINATION)
            {
                break;
            }
        }
        //write chemical and additional RF information to the trajectory:
        hdf5Writer->writeTrajectoryAttribute("Substance Names", discreteSubstanceNames);
        hdf5Writer->writeTrajectoryAttribute("Substance Masses", discreteSubstanceMasses);
        if (rfMode==RAMPED_RF) {
            hdf5Writer->writeNumericListDataset("V_rf", V_rf_export);
        }
        verletIntegrator.finalizeSimulation();
        concentrationFilewriter.closeFile();

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
