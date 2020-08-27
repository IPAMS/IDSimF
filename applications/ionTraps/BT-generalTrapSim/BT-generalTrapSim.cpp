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
 BT-generalTrapSim.cpp

 Ion trajectory simulation (including space charge and hard sphere collisions) in an RF trap with
 arbitrary geometry given by potential arrays

 ****************************/

#include "json.h"
#include "parameterParsing.hpp"
#include "BTree_particle.hpp"
#include "BTree_tree.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_scalar_writer.hpp"
#include "PSim_util.hpp"
#include "PSim_verletIntegrator.hpp"
#include "PSim_parallelVerletIntegrator.hpp"
#include "PSim_sampledWaveform.hpp"
#include "PSim_math.hpp"
#include "PSim_averageChargePositionWriter.hpp"
#include "PSim_inductionCurrentWriter.hpp"
#include "PSim_ionCloudReader.hpp"
#include "PSim_simionPotentialArray.hpp"
#include "CollisionModel_HardSphere.hpp"
#include <iostream>
#include <vector>
#include <ctime>
#include <filesystem>

enum IntegratorMode {VERLET,PARALLEL_VERLET};
enum IonStartGeometry {BOX,CYLINDER};
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

    // read configuration file ======================================================================
    if (argc <2){
        std::cout << "no conf project name or conf file given"<<std::endl;
        return(0);
    }

    std::string confFileName = argv[1];
    std::filesystem::path confBasePath = std::filesystem::path(confFileName).parent_path();

    std::string projectName = argv[2];
    std::cout << projectName<<std::endl;

    Json::Value confRoot = readConfigurationJson(confFileName);

    // read basic simulation parameters =============================================================
    std::string integratorMode_str = stringConfParameter("integrator_mode", confRoot);
    IntegratorMode integratorMode;
    if (integratorMode_str == "verlet"){
        integratorMode = VERLET;
    }
    else if(integratorMode_str== "parallel_verlet"){
        integratorMode = PARALLEL_VERLET;
    }
    else {
        throw std::invalid_argument("wrong configuration value: integrator mode");
    }

    int timeSteps = intConfParameter("sim_time_steps", confRoot);
    int trajectoryWriteInterval = intConfParameter("trajectory_write_interval", confRoot);
    int fftWriteInterval = intConfParameter("fft_write_interval",confRoot);
    double dt = doubleConfParameter("dt", confRoot);
    std::vector<int> nIons = std::vector<int>();
    std::vector<double> ionMasses = std::vector<double>();
    std::vector<double> ionCollisionDiameters_angstrom = std::vector<double>();
    std::string fftWriteMode_str = stringConfParameter("fft_write_mode",confRoot);
    FftWriteMode fftWriteMode;
    if (fftWriteMode_str == "unresolved"){
        fftWriteMode = UNRESOLVED;
    }else if(fftWriteMode_str== "mass_resolved") {
        fftWriteMode = MASS_RESOLVED;
    }

    //read potential array configuration of the trap =================================================
    double paSpatialScale = doubleConfParameter("potential_array_scaling", confRoot);
    std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> potentialArrays;
    std::vector<std::string> potentialArraysNames = stringVectorConfParameter("potential_arrays", confRoot);
    for(const auto &paName: potentialArraysNames){
        std::filesystem::path paPath = confBasePath / paName;
        std::unique_ptr<ParticleSimulation::SimionPotentialArray> pa_pt =
                std::make_unique<ParticleSimulation::SimionPotentialArray>(paPath, paSpatialScale);
        potentialArrays.push_back(std::move(pa_pt));
    }


    // SIMION fast adjust PAs use 10000 as normalized potential value, thus we have to scale everything with 1/10000
    double potentialScale = 1.0 / 10000.0;
    std::vector<double> potentialsFactorsDc = doubleVectorConfParameter("dc_potentials", confRoot, potentialScale);
    std::vector<double> potentialFactorsRf = doubleVectorConfParameter("rf_potential_factors", confRoot, potentialScale);
    std::vector<double> potentialFactorsExcite = doubleVectorConfParameter("excite_potential_factors", confRoot, potentialScale);
    std::vector<double> detectionPAFactorsRaw = doubleVectorConfParameter("detection_potential_factors", confRoot);
    std::vector<ParticleSimulation::SimionPotentialArray*> detectionPAs;

    std::vector<double> detectionPAFactors;
    for(size_t i=0; i < detectionPAFactorsRaw.size(); ++i){
        double paKey = detectionPAFactorsRaw[i];
        if (paKey != 0.0){
            detectionPAs.emplace_back(potentialArrays[i].get());
            detectionPAFactors.emplace_back(paKey);
        }
    }

    // defining simulation domain box (used for ion termination):
    std::array<std::array<double,2>,3> simulationDomainBoundaries;
    if (confRoot.isMember("simulation_domain_boundaries")){
        // get manual simulation domain boundaries from config file
        simulationDomainBoundaries = double3dBox("simulation_domain_boundaries", confRoot);
    }
    else {
        // use minimum PA extent box as domain boundaries
        std::array<double, 6> minExtent = potentialArrays[0]->getBounds();
        for(const auto &pa: potentialArrays){
            std::array<double, 6> paBounds = pa->getBounds();
            for (size_t i=0; i<6; ++i){
                if (minExtent[i] < paBounds[i]){
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
    double backgroundPressure = doubleConfParameter("background_pressure_Pa", confRoot);
    double backgroundTemperature = doubleConfParameter("background_temperature_K", confRoot);
    double spaceChargeFactor = doubleConfParameter("space_charge_factor", confRoot);
    double collisionGasMassAmu = doubleConfParameter("collision_gas_mass_amu", confRoot);
    double collisionGasDiameterM = doubleConfParameter("collision_gas_diameter_angstrom", confRoot)*1e-10;


    //read rf configuration =========================================================================
    double f_rf= doubleConfParameter("f_rf", confRoot); //RF frequency 1e6;
    double omega = f_rf * 2.0 * M_PI; //RF angular frequencyf_rf* 2.0 * M_PI;

    RfAmplitudeMode rfMode;
    std::vector<double> V_0_ramp;
    double V_0 = 0.0;
    if (confRoot.isMember("V_rf_start")){
        rfMode = RAMPED_RF;
        double V_rf_start = doubleConfParameter("V_rf_start", confRoot);
        double V_rf_end = doubleConfParameter("V_rf_end", confRoot);
        V_0_ramp = ParticleSimulation::linspace(V_rf_start,V_rf_end,timeSteps);
    }
    else {
        rfMode = STATIC_RF;
        V_0 = doubleConfParameter("V_rf", confRoot);
    }
    std::vector<double> V_rf_export;


    //read excitation / swift configuration ========================================================
    ExciteMode exciteMode;
    std::unique_ptr<ParticleSimulation::SampledWaveform> swiftWaveForm;
    double excitePulseLength = 0.0;
    if (confRoot.isMember("excite_waveform_csv_file") ){
        exciteMode = SWIFT;
        std::string swiftFileName = confRoot.get("excite_waveform_csv_file",0).asString();
        swiftWaveForm = std::make_unique<ParticleSimulation::SampledWaveform>(swiftFileName);
        if (! swiftWaveForm->good()){
            std::cout << "swift transient file not accessible"<<std::endl;
            return(0);
        }
    }
    else {
        exciteMode = RECTPULSE;
        excitePulseLength = doubleConfParameter("excite_pulse_length", confRoot);
    }
    double excitePulsePotential = doubleConfParameter("excite_pulse_potential", confRoot);



    //read ion configuration =======================================================================
    std::vector<std::unique_ptr<BTree::Particle>>particles;
    std::vector<BTree::Particle*>particlePtrs;

    if (confRoot.isMember("ion_cloud_init_file") == true) {
        std::string ionCloudFileName = confRoot.get("ion_cloud_init_file", 0).asString();
        ParticleSimulation::IonCloudReader reader = ParticleSimulation::IonCloudReader();
        particles = reader.readIonCloud(ionCloudFileName);
        //prepare a vector of raw pointers
        for (const auto& part : particles){
            particlePtrs.push_back(part.get());
        }

    } else {
        // ions are not given in an init file, read and init random ion box configuration
        if (confRoot.isMember("n_ions") == true) {
            Json::Value n_ions_json = confRoot.get("n_ions", 0);
            for (int i = 0; i < n_ions_json.size(); i++) {
                nIons.push_back(n_ions_json.get(i, 0.0).asInt());
            }
        } else {
            throw std::invalid_argument("missing configuration value: n_ions");
        }

        if (confRoot.isMember("ion_masses") == true) {
            Json::Value ions_masses_json = confRoot.get("ion_masses", 0);
            for (int i = 0; i < ions_masses_json.size(); i++) {
                ionMasses.push_back(ions_masses_json.get(i, 0.0).asDouble());
            }
        } else {
            throw std::invalid_argument("missing configuration value: ion_masses");
        }

        if (confRoot.isMember("ion_collision_gas_diameters_angstrom") == true) {
            //Fixme: Use vector library functions
            Json::Value ions_collision_diams = confRoot.get("ion_collision_gas_diameters_angstrom", 0);
            for (int i = 0; i < ions_collision_diams.size(); i++) {
                ionCollisionDiameters_angstrom.push_back(ions_collision_diams.get(i, 0.0).asDouble());
            }
        } else {
            throw std::invalid_argument("missing configuration value: ion_collision_gas_diameters_angstrom");
        }

        double ions_tob_range = doubleConfParameter("ion_time_of_birth_range_s", confRoot);

        std::string ionStartGeom_str = stringConfParameter("ion_start_geometry",confRoot);
        IonStartGeometry ionStartGeom;

        double ionStartCylinder_radius;
        double ionStartCylinder_length;
        Core::Vector ionsBasePos;
        std::vector<double> basePos = doubleVectorConfParameter("ion_start_base_position_m", confRoot);
        ionsBasePos.x(basePos[0]);
        ionsBasePos.y(basePos[1]);
        ionsBasePos.z(basePos[2]);

        if (ionStartGeom_str == "box"){
            ionStartGeom = BOX;
        } else if (ionStartGeom_str == "cylinder"){
            ionStartGeom = CYLINDER;
            ionStartCylinder_radius = doubleConfParameter("ion_start_cylinder_radius_m", confRoot);
            ionStartCylinder_length = doubleConfParameter("ion_start_cylinder_length_m", confRoot);
        }

        for (int i = 0; i < nIons.size(); i++) {
            int nParticles = nIons[i];
            double mass = ionMasses[i];
            double collisionDiameter_m = ionCollisionDiameters_angstrom[i]*1e-10;
            std::vector<std::unique_ptr<BTree::Particle>> ions;

            if (ionStartGeom == BOX) {
                ions = ParticleSimulation::util::getRandomIonsInBox(
                        nParticles, 1.0,
                        Core::Vector(-1.5, -1.5, -1.5) / 1000.0,
                        Core::Vector(3, 3, 3) / 1000.0,
                        ions_tob_range);

            }
            else if (ionStartGeom == CYLINDER){
                ions = ParticleSimulation::util::getRandomIonsInCylinderXDirection(
                        nParticles, 1.0, ionStartCylinder_radius, ionStartCylinder_length, ions_tob_range);
            }

            for (int j = 0; j < nParticles; j++) {
                ions[j]->setMassAMU(mass);
                ions[j]->setDiameter(collisionDiameter_m);

                //shift all ions according to the base pos:
                Core::Vector buf = ions[j]->getLocation() + ionsBasePos;
                ions[j]->setLocation(buf);
                particlePtrs.push_back(ions[j].get());
                particles.push_back(std::move(ions[j]));
            }
        }
    }


    // define functions for the trajectory integration ==================================================
    int ionsInactive = 0;
    auto trapFieldFunction =
                 [exciteMode, rfMode, excitePulseLength, excitePulsePotential,
                         spaceChargeFactor, omega, &swiftWaveForm, &V_0, &V_0_ramp, &ionsInactive,
                         &potentialArrays, &potentialsFactorsDc, &potentialFactorsRf, &potentialFactorsExcite]
                         (BTree::Particle *particle, int particleIndex,auto& tree, double time, int timestep)
                         -> Core::Vector{

                     Core::Vector pos = particle->getLocation();
                     double particleCharge = particle->getCharge();
                     double excitePotential = 0;
                     if (exciteMode == RECTPULSE) {
                         if (time < excitePulseLength) {
                             excitePotential = excitePulsePotential;
                         }
                     } else {
                         excitePotential = swiftWaveForm->getValue(timestep) * excitePulsePotential;
                     }

                     if (rfMode == RAMPED_RF) {
                         V_0 = V_0_ramp[timestep];
                     }

                     Core::Vector fEfield(0, 0, 0);
                     double V_t = V_0 * cos(omega * time);


                     for (size_t i = 0; i < potentialArrays.size(); ++i) {
                         Core::Vector paField = potentialArrays[i]->getField(pos.x(), pos.y(), pos.z());

                         Core::Vector paEffectiveField =
                                 paField * potentialsFactorsDc.at(i) +
                                 paField * potentialFactorsRf.at(i) * V_t +
                                 paField * potentialFactorsExcite.at(i) * excitePotential;

                         fEfield = fEfield + paEffectiveField;
                     }

                     return fEfield * particleCharge;
                 };


    auto accelerationFunctionQIT =
            [spaceChargeFactor, &trapFieldFunction](
                    BTree::Particle *particle, int particleIndex,
                    auto& tree, double time, int timestep) -> Core::Vector{

                Core::Vector pos = particle->getLocation();
                double particleCharge = particle->getCharge();

                Core::Vector trapForce = trapFieldFunction(particle, particleIndex, tree, time, timestep);

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle) * (particleCharge * spaceChargeFactor);
                }

                //update the additional parameters for writing them later to the trajectory:
                particle->setAuxScalarParam(key_trapForce_x, trapForce.x());
                particle->setAuxScalarParam(key_trapForce_y, trapForce.y());
                particle->setAuxScalarParam(key_trapForce_z, trapForce.z());
                particle->setAuxScalarParam(key_spaceCharge_x, spaceChargeForce.x());
                particle->setAuxScalarParam(key_spaceCharge_y, spaceChargeForce.y());
                particle->setAuxScalarParam(key_spaceCharge_z, spaceChargeForce.z());

                return ((trapForce + spaceChargeForce) / particle->getMass());
            };

    auto accelerationFunctionQIT_parallel =
            [spaceChargeFactor, &trapFieldFunction](
                    BTree::Particle *particle, int particleIndex,
                    auto& tree, double time, int timestep) -> Core::Vector{

                Core::Vector pos = particle->getLocation();
                double particleCharge = particle->getCharge();

                Core::Vector rfForce = trapFieldFunction(particle, particleIndex, tree, time, timestep);

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle) * (particleCharge * spaceChargeFactor);
                }

                //update the additional parameters for writing them later to the trajectory:
                particle->setAuxScalarParam(key_trapForce_x, rfForce.x());
                particle->setAuxScalarParam(key_trapForce_y, rfForce.y());
                particle->setAuxScalarParam(key_trapForce_z, rfForce.z());
                particle->setAuxScalarParam(key_spaceCharge_x, spaceChargeForce.x());
                particle->setAuxScalarParam(key_spaceCharge_y, spaceChargeForce.y());
                particle->setAuxScalarParam(key_spaceCharge_z, spaceChargeForce.z());

                return ((rfForce + spaceChargeForce) / particle->getMass());
            };

    auto otherActionsFunctionQIT = [&simulationDomainBoundaries, &ionsInactive](
            Core::Vector& newPartPos, BTree::Particle* particle,
            int particleIndex,
            auto& tree, double time, int timestep) {
        // if the ion is out of the boundary box: Terminate the ion
        if ( newPartPos.x() <= simulationDomainBoundaries[0][0] ||
             newPartPos.x() >= simulationDomainBoundaries[0][1] ||
             newPartPos.y() <= simulationDomainBoundaries[1][0] ||
             newPartPos.y() >= simulationDomainBoundaries[1][1] ||
             newPartPos.z() <= simulationDomainBoundaries[2][0] ||
             newPartPos.z() >= simulationDomainBoundaries[2][1] )
        {
            particle->setActive(false);
            particle->setSplatTime(time);
            ionsInactive++;
        }
    };

    //prepare file writers and data writing functions ==============================================================================
    auto avgPositionWriter = std::make_unique<ParticleSimulation::AverageChargePositionWriter>(projectName+"_averagePosition.txt");

    auto fftWriter = std::make_unique<ParticleSimulation::InductionCurrentWriter>(
            particlePtrs, projectName + "_fft.txt",detectionPAs,detectionPAFactors,paSpatialScale);
    auto ionsInactiveWriter = std::make_unique<ParticleSimulation::Scalar_writer>(projectName+ "_ionsInactive.txt");

    ParticleSimulation::additionalPartParamFctType additionalParameterTransformFct =
            [](BTree::Particle *particle) -> std::vector<double>{
                std::vector<double> result = {
                        particle->getVelocity().x(),
                        particle->getVelocity().y(),
                        particle->getVelocity().z(),
                        particle->getAuxScalarParam(key_trapForce_x),
                        particle->getAuxScalarParam(key_trapForce_y),
                        particle->getAuxScalarParam(key_trapForce_z),
                        particle->getAuxScalarParam(key_spaceCharge_x),
                        particle->getAuxScalarParam(key_spaceCharge_y),
                        particle->getAuxScalarParam(key_spaceCharge_z),
                };
                return result;
            };

    std::vector<std::string> auxParamNames = {"velocity x","velocity y","velocity z",
                                              "rf x","rf y","rf z",
                                              "spacecharge x","spacecharge y","spacecharge z"};

    auto hdf5Writer = std::make_unique<ParticleSimulation::TrajectoryHDF5Writer>(
            projectName + "_trajectories.hd5", auxParamNames, additionalParameterTransformFct);

    auto timestepWriteFunction =
            [trajectoryWriteInterval, fftWriteInterval, fftWriteMode, &V_0, &V_rf_export, &ionsInactive, timeSteps,
             &hdf5Writer, &additionalParameterTransformFct, &avgPositionWriter, &ionsInactiveWriter, &fftWriter](
                    std::vector<BTree::Particle*>& particles, auto& tree, double time, int timestep, bool lastTimestep){

                if (timestep % fftWriteInterval == 0) {
                    //avgPositionWriter->writeTimestep(tree, time);
                    ionsInactiveWriter->writeTimestep(ionsInactive, time);
                    if (fftWriteMode == UNRESOLVED){
                        fftWriter->writeTimestep(time);
                    }
                    else if(fftWriteMode == MASS_RESOLVED){
                        //fftWriter->writeTimestepMassResolved(time);
                        std::stringstream ss;
                        ss << "Mass resolved induction current not implemented";
                        throw (std::runtime_error(ss.str()));
                    }
                }

                if (lastTimestep) {
                    V_rf_export.emplace_back(V_0);
                    hdf5Writer->writeTimestep(particles,time);
                    hdf5Writer->finalizeTrajectory();
                    std::cout << "finished ts:" << timestep << " time:" << time << std::endl;
                }
                else if (timestep % trajectoryWriteInterval == 0) {

                    std::cout << "ts:" << timestep << " time:" << time << " V_rf:" << V_0
                              <<" ions existing: "<<particles.size()<< " ions inactive: "
                              << ionsInactive << std::endl;

                    if (particles.size() > 0) {
                        V_rf_export.emplace_back(V_0);
                        hdf5Writer->writeTimestep(particles, time);
                    }
                }
    };

    CollisionModel::HardSphereModel hsModel(backgroundPressure,
                                            backgroundTemperature,
                                            collisionGasMassAmu,
                                            collisionGasDiameterM);



    // simulate ===============================================================================================
    clock_t begin = std::clock();

    if (integratorMode == VERLET) {
        ParticleSimulation::VerletIntegrator verletIntegrator(
                particlePtrs,
                accelerationFunctionQIT, timestepWriteFunction, otherActionsFunctionQIT,
                hsModel);
        verletIntegrator.run(timeSteps, dt);
    }
    else if (integratorMode == PARALLEL_VERLET) {
        ParticleSimulation::ParallelVerletIntegrator verletIntegrator(
                particlePtrs,
                accelerationFunctionQIT_parallel, timestepWriteFunction, otherActionsFunctionQIT,
                hsModel);
        verletIntegrator.run(timeSteps, dt);
    }


    if (rfMode == RAMPED_RF) {
        hdf5Writer->writeVectorDataset("V_rf",V_rf_export);
    }

    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    std::cout << particles[0]->getLocation()<<std::endl;
    std::cout << "elapsed secs"<< elapsed_secs<<std::endl;
    return 0;
}
