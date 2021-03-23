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

#include "json.h"
#include "appUtils_parameterParsing.hpp"
#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_ConcentrationFileWriter.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_scalar_writer.hpp"
#include "PSim_util.hpp"
#include "PSim_verletIntegrator.hpp"
#include "PSim_sampledWaveform.hpp"
#include "PSim_math.hpp"
#include "PSim_averageChargePositionWriter.hpp"
#include "PSim_idealizedQitFFTWriter.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "CollisionModel_MultiCollisionModel.hpp"
#include <iostream>
#include <vector>
#include <ctime>

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
    return [=](double collisionEnergy,BTree::Particle& ion)->void{
        int nCollisions = ion.getAuxScalarParam(key);
        ion.setAuxScalarParam(key, nCollisions+1);
    };
}

std::function<void(RS::CollisionConditions, BTree::Particle &)>
createCollisionReactionFunction(RS::Substance *collisionPartnerSubstance, RS::Simulation &rsSim, std::string key){

    return [collisionPartnerSubstance,key,&rsSim] (RS::CollisionConditions conditions, BTree::Particle& ion){
        rsSim.collisionReact(ion.getIndex(),collisionPartnerSubstance,conditions);
        int nCollisions = ion.getAuxScalarParam(key);
        ion.setAuxScalarParam(key, nCollisions+1);
    };
}


int main(int argc, const char * argv[]) {

    // read configuration file ======================================================================
    if (argc <2){
        std::cout << "no conf project name or conf file given"<<std::endl;
        return(1);
    }

    std::string confFileName = argv[1];
    std::cout << confFileName<<std::endl;

    std::string projectName = argv[2];
    std::cout << projectName<<std::endl;

    Json::Value confRoot = readConfigurationJson(confFileName);
    std::cout<<confRoot<<std::endl;


    // read basic simulation parameters =============================================================
    int timeSteps = intConfParameter("sim_time_steps", confRoot);
    int concentrationWriteInterval = intConfParameter("concentrations_write_interval",confRoot);

    int trajectoryWriteInterval = intConfParameter("trajectory_write_interval", confRoot);
    int fftWriteInterval = intConfParameter("fft_write_interval",confRoot);
    double dt = doubleConfParameter("dt", confRoot);
    std::vector<int> nIons = std::vector<int>();
    std::vector<double> ionMasses = std::vector<double>();
    std::string fftWriteMode_str = stringConfParameter("fft_mode",confRoot);
    FftWriteMode fftWriteMode;
    if(fftWriteMode_str== "off") {
        fftWriteMode = OFF;
    }else if(fftWriteMode_str == "unresolved"){
        fftWriteMode = UNRESOLVED;
    }else if(fftWriteMode_str== "mass_resolved") {
        fftWriteMode = MASS_RESOLVED;
    }
    else{
        throw std::invalid_argument("wrong configuration value: fft_mode");
    }



    //read geometrical configuration of the trap======================================================
    double r_0 = 0.0;
    double z_0 = 0.0;

    std::string geometryMode_str = stringConfParameter("geometry_mode",confRoot);
    GeometryMode geometryMode;
    if (geometryMode_str == "default"){
        geometryMode= DEFAULT;
        r_0 = r_0_default;
        z_0 = z_0_default;
    }else if(geometryMode_str== "variable") {
        geometryMode= VARIABLE;
        r_0 = doubleConfParameter("r_0",confRoot);
        z_0 = doubleConfParameter("z_0",confRoot);
    }
    else if(geometryMode_str== "scaled") {
        geometryMode= SCALED;
        double geomScale = doubleConfParameter("geometry_scale",confRoot);
        r_0 = r_0_default*geomScale;
        z_0 = z_0_default*geomScale;
    }
    else {
        throw std::invalid_argument("wrong configuration value: geometry_mode");
    }

    double d_square = r_0*r_0 + 2*z_0*z_0;
    double d_square_2 = d_square/2.0;
    double U_0 = 0.0;


    //read physical configuration ===================================================================
    double maxIonRadius_m = doubleConfParameter("max_ion_radius_m",confRoot);
    double spaceChargeFactor = doubleConfParameter("space_charge_factor", confRoot);
    double backgroundTemperature_K = doubleConfParameter("background_temperature_K", confRoot);

    //read collision gas / background gas configuration =============================================
    std::vector<std::string> collisionGasNames = stringVectorConfParameter("collision_gas_names", confRoot);
    std::vector<double> partialPressures = doubleVectorConfParameter("partial_pressures_Pa", confRoot);
    std::vector<double> collisionGasMasses_Amu = doubleVectorConfParameter("collision_gas_masses_amu", confRoot);
    std::vector<double> collisionGasDiameters_m = doubleVectorConfParameter("collision_gas_diameters_angstrom", confRoot);
    collisionGasDiameters_m[0]*=1e-10;
    collisionGasDiameters_m[1]*=1e-10;

    double totalBackgroundPressure_Pa = partialPressures[0] + partialPressures[1];

    double startWidth_m = doubleConfParameter("start_width_m",confRoot);

    //read rf configuration =========================================================================
    double f_rf= doubleConfParameter("f_rf", confRoot); //RF frequency 1e6;
    double omega = f_rf * 2.0 * M_PI; //RF angular frequencyf_rf* 2.0 * M_PI;

    RfAmplitudeMode rfMode;
    std::vector<double> V_0_ramp;
    double V_0 = 0.0;
    if (confRoot.isMember("rf_ramp_start_V")){
        rfMode = RAMPED_RF;
        int V_rf_waiting_ts = intConfParameter("rf_ramp_waiting_timesteps", confRoot);
        double V_rf_start = doubleConfParameter("rf_ramp_start_V", confRoot);
        double V_rf_stop = doubleConfParameter("rf_ramp_stop_V", confRoot);

        auto vRf = ParticleSimulation::fillVector(V_rf_start,V_rf_waiting_ts);
        auto vRamp = ParticleSimulation::linspace(V_rf_start,V_rf_stop,timeSteps-V_rf_waiting_ts);
        vRf.insert( vRf.end(), vRamp.begin(), vRamp.end() ); //concat vRf and vRamp
        V_0_ramp = vRf;
    }
    else {
        rfMode = STATIC_RF;
        V_0 = doubleConfParameter("rf_V", confRoot);
    }


    //read excitation / swift configuration ========================================================
    ExciteMode exciteMode;
    std::unique_ptr<ParticleSimulation::SampledWaveform> swiftWaveForm;
    double excitePulseLength = 0.0;
    double exciteDivisor = 0.0;


    std::string exciteMode_str = stringConfParameter("excite_mode",confRoot);
    if(exciteMode_str== "off") {
        exciteMode= NOEXCITE;
    }else if(exciteMode_str == "rect_pulse"){
        excitePulseLength = doubleConfParameter("excite_pulse_length", confRoot);
        exciteMode= RECTPULSE;
    }else if(exciteMode_str== "waveform") {
        exciteMode= SWIFT;
        if (confRoot.isMember("excite_waveform_csv_file") ){
            exciteMode = SWIFT;
            std::string swiftFileName = confRoot.get("excite_waveform_csv_file",0).asString();
            swiftWaveForm = std::make_unique<ParticleSimulation::SampledWaveform>(swiftFileName);
            if (! swiftWaveForm->good()){
                std::cout << "swift transient file not accessible"<<std::endl;
                return(0);
            }
        }
    }else if(exciteMode_str== "continuous_sine") {
        exciteMode= CONTINUOUSSINE;
        exciteDivisor = doubleConfParameter("excite_divisor", confRoot);
    }else{
        throw std::invalid_argument("wrong configuration value: excite_mode");
    }
    double excitePotential;
    if (exciteMode != NOEXCITE) {
        excitePotential = doubleConfParameter("excite_potential", confRoot);
    }

    //read and prepare chemical configuration ===============================================
    RS::ConfigFileParser parser = RS::ConfigFileParser();
    std::string rsFilePath = pathRelativeToConfFile(
            confFileName,stringConfParameter("reaction_configuration",confRoot) );
    RS::Simulation rsSim = RS::Simulation(parser.parseFile(rsFilePath));
    RS::SimulationConfiguration* simConf = rsSim.simulationConfiguration();
    //prepare a map for retrieval of the substances from their index:
    std::map<RS::Substance*,int> substanceIndices;

    //prepare some parameter vectors required for the hdf5 trajectory writer:
    std::vector<std::string> discreteSubstanceNames;
    std::vector<double> discreteSubstanceMasses;
    std::vector<double> V_rf_export;

    std::vector<RS::Substance*> discreteSubstances = simConf->getAllDiscreteSubstances();
    for (int i=0; i<discreteSubstances.size(); i++){
        substanceIndices.insert(std::pair<RS::Substance*,int>(discreteSubstances[i], i));
        discreteSubstanceNames.emplace_back(discreteSubstances[i]->name());
        discreteSubstanceMasses.emplace_back(discreteSubstances[i]->mass());
    }


    //read ion configuration =======================================================================
    int nParticlesTotal = 0;
    std::vector<uniqueReactivePartPtr>particles;
    std::vector<BTree::Particle*>particlePtrs;

    // read and init random ion box configuration
    if (confRoot.isMember("n_ions") == true) {
        Json::Value n_ions_json = confRoot.get("n_ions", 0);
        for (int i = 0; i < n_ions_json.size(); i++) {
            nIons.push_back(n_ions_json.get(i, 0.0).asInt());
        }
    } else {
        throw std::invalid_argument("missing configuration value: n_ions");
    }
    Core::Vector initCorner(-startWidth_m/2.0,-startWidth_m/2.0,-startWidth_m/2.0);
    Core::Vector initBoxSize(startWidth_m,startWidth_m,startWidth_m);

    for (int i=0; i<nIons.size();i++) {
        RS::Substance *subst = discreteSubstances[i];
        std::vector<Core::Vector> initialPositions =
                ParticleSimulation::util::getRandomPositionsInBox(nIons[i],initCorner,initBoxSize);
        for (int k = 0; k < nIons[i]; k++) {
            uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(subst);

            particle->setIndex(nParticlesTotal);
            particle->setLocation(initialPositions[k]);
            particle->setAuxScalarParam(key_Collisions_total,0);
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
    int lastTimestep = -1;
    double parameter_a;
    double exciteCos;

    auto accelerationFunctionQIT =
            [exciteMode, rfMode, excitePulseLength, excitePotential,
             spaceChargeFactor, omega, z_0, U_0, d_square_2,
             &swiftWaveForm, exciteDivisor, &V_0, &V_0_ramp, &lastTimestep, &parameter_a,&exciteCos](
                    BTree::Particle *particle, int particleIndex,
                    BTree::Tree &tree, double time, int timestep) -> Core::Vector{

                Core::Vector pos = particle->getLocation();
                double particleCharge = particle->getCharge();
                double tsExcitePotential = 0;
                if (exciteMode == RECTPULSE) {
                    if (time < excitePulseLength) {
                        tsExcitePotential = excitePotential;
                    }
                } else if (exciteMode == SWIFT) {
                    tsExcitePotential = swiftWaveForm->getValue(timestep) * excitePotential;
                } else if (exciteMode == CONTINUOUSSINE){
                    if (timestep > lastTimestep) {
                        exciteCos = cos(omega / exciteDivisor * time)* excitePotential;
                    }
                    tsExcitePotential = exciteCos;
                }

                if (rfMode == RAMPED_RF) {
                    V_0 = V_0_ramp[timestep];
                }

                double a_ex = tsExcitePotential / z_0;

                Core::Vector rfForce(0,0,0);
                if (timestep > lastTimestep) {
                    parameter_a = (U_0 + (V_0 * cos(omega * time))) / d_square_2;
                    lastTimestep = timestep;
                }
                rfForce = Core::Vector(
                        parameter_a * pos.x(),
                        parameter_a * pos.y(),
                        -2 * parameter_a * pos.z() + a_ex) *particleCharge;

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle) * (particleCharge * spaceChargeFactor);
                }

                //update the additional parameters for writing them later to the trajectory:
                //particle->setAuxScalarParam(key_trapForce_x, rfForce.x());
                return ((rfForce + spaceChargeForce) / particle->getMass());
            };

    int ionsInactive = 0;
    auto otherActionsFunctionQIT = [maxIonRadius_m, &ionsInactive](Core::Vector &newPartPos, BTree::Particle *particle,
                                                                 int particleIndex,
                                                                 BTree::Tree &tree, double time, int timestep){
        if (newPartPos.magnitude() > maxIonRadius_m) {
            particle->setActive(false);
            particle->setSplatTime(time);
            ionsInactive++;
        }
    };

    //prepare file writers and data writing functions ==============================================================================
    auto avgPositionWriter = std::make_unique<ParticleSimulation::AverageChargePositionWriter>(projectName+"_averagePosition.txt");
    std::unique_ptr<ParticleSimulation::IdealizedQitFFTWriter>fftWriter = nullptr;
    if (fftWriteMode != OFF){
        fftWriter = std::make_unique<ParticleSimulation::IdealizedQitFFTWriter>(particlePtrs, projectName + "_fft.txt");
    }
    auto ionsInactiveWriter = std::make_unique<ParticleSimulation::Scalar_writer>(projectName+ "_ionsInactive.txt");

    RS::ConcentrationFileWriter concentrationFilewriter(projectName+ "_concentrations.txt");
    concentrationFilewriter.initFile(simConf);

    ParticleSimulation::additionalPartParamFctType additionalParameterTransformFct =
            [](BTree::Particle *particle) -> std::vector<double>{
                double ionVelocity = particle->getVelocity().magnitude();
                double kineticEnergy_eV = 0.5* particle->getMass() * ionVelocity * ionVelocity * Core::JOULE_TO_EV;
                std::vector<double> result = {
                        particle->getVelocity().x(),
                        particle->getVelocity().y(),
                        particle->getVelocity().z(),
                        kineticEnergy_eV,
                        particle->getAuxScalarParam(key_Collisions_total),
                        particle->getAuxScalarParam(key_ChemicalIndex)
                };
                return result;
            };

    std::vector<std::string> auxParamNames = {"velocity x", "velocity y", "velocity z", "kinetic energy (eV)", "total collisions",
                                              "chemical id"};

    auto hdf5Writer = std::make_unique<ParticleSimulation::TrajectoryHDF5Writer>(
            projectName + "_trajectories.hd5", auxParamNames,additionalParameterTransformFct);

    auto timestepWriteFunction =
            [trajectoryWriteInterval, fftWriteInterval, fftWriteMode, &V_0, &V_rf_export, &ionsInactive, timeSteps,
             &hdf5Writer, &avgPositionWriter, &ionsInactiveWriter, &fftWriter](
                    std::vector<BTree::Particle*>& particles,BTree::Tree& tree, double time, int timestep, bool lastTimestep){

                if (timestep % fftWriteInterval == 0) {
                    avgPositionWriter->writeTimestep(tree, time);
                    ionsInactiveWriter->writeTimestep(ionsInactive, time);
                    if (fftWriteMode == UNRESOLVED){
                        fftWriter->writeTimestep(time);
                    }
                    else if(fftWriteMode == MASS_RESOLVED){
                        fftWriter->writeTimestepMassResolved(time);
                    }
                }

                if (lastTimestep) {
                    V_rf_export.emplace_back(V_0);
                    hdf5Writer->writeTimestep(particles,time);
                    hdf5Writer->finalizeTrajectory();
                    std::cout << "finished ts:" << timestep << " time:" << time << std::endl;
                }

                else if (timestep % trajectoryWriteInterval == 0) {
                    std::cout << "ts:" << timestep << " time:" << time << " V_rf:" << V_0 << " ions inactive:" << ionsInactive
                              << std::endl;
                    V_rf_export.emplace_back(V_0);
                    hdf5Writer->writeTimestep(particles,time);

                }
    };

    //prepare background gas collision models and collision based chemical reactions
    std::vector<std::unique_ptr<CollisionModel::AbstractCollisionModel>> collisionModels;

    for (int i=0; i<collisionGasNames.size(); ++i){
        RS::Substance* collisionPartnerSubst = simConf->substanceByName(collisionGasNames[i]);
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
            accelerationFunctionQIT, timestepWriteFunction, otherActionsFunctionQIT,
            combinedCollisionModel);

    clock_t begin = std::clock();

    for (int step=0; step<timeSteps; ++step) {
        if (step % concentrationWriteInterval == 0) {
            rsSim.printConcentrations();
            concentrationFilewriter.writeTimestep(rsSim);
        }
        for (int i = 0; i < nParticlesTotal; i++) {
            if (particles[i]->isActive()) {
                bool reacted = rsSim.react(i, reactionConditions, dt);
                int substIndex = substanceIndices.at(particles[i]->getSpecies());
                particles[i]->setAuxScalarParam(key_ChemicalIndex, substIndex);
            }
        }
        rsSim.advanceTimestep(dt);
        verletIntegrator.runSingleStep(dt);

        if (ionsInactive >= nParticlesTotal){
            break;
        }
    }
    //write chemical and additional RF information to the trajectory:
    hdf5Writer->writeTrajectoryAttribute("Substance Names",discreteSubstanceNames);
    hdf5Writer->writeTrajectoryAttribute("Substance Masses",discreteSubstanceMasses);
    if (rfMode == RAMPED_RF) {
        hdf5Writer->writeVectorDataset("V_rf",V_rf_export);
    }
    verletIntegrator.finalizeSimulation();

    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    std::cout << particles[0]->getLocation()<<std::endl;
    std::cout << "elapsed secs"<< elapsed_secs<<std::endl;
    return 0;
}
