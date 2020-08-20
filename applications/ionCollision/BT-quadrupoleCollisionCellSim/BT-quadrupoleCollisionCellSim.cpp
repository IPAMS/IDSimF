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
 BT-quadrupoleCollisionCellSim.cpp

 Ion trajectory simulation (including space charge and hard sphere collisions) in an
 quadrupole collision cell with hard sphere collision dynamics

 ****************************/

#include "Core_randomGenerators.hpp"
#include "BTree_particle.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_util.hpp"
#include "PSim_parallelVerletIntegrator.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "json.h"
#include "parameterParsing.hpp"
#include "inputFileUtilities.hpp"
#include <iostream>
#include <vector>
#include <ctime>

enum IonStartGeometry {BOX,CYLINDER,VECTOR};

// constants:
const double rho_per_pa = 2.504e20; //(particles / m^3) / Pa

int main(int argc, const char * argv[]) {

    // read configuration file ======================================================================
    if (argc <2){
        std::cout << "no conf project name or conf file given"<<std::endl;
        return(0);
    }

    std::string confFileName = argv[1];
    std::string confBasePath = confFileBasePath(confFileName);
    std::cout << confFileName << std::endl;

    Json::Value confRoot = readConfigurationJson(confFileName);
    std::cout<<confRoot<<std::endl;

    std::string projectName = argv[2];
    std::cout << projectName << std::endl;

    // read basic simulation parameters =============================================================
    int timeSteps = intConfParameter("sim_time_steps",confRoot);
    int trajectoryWriteInterval = intConfParameter("trajectory_write_interval", confRoot);
    double dt = doubleConfParameter("dt", confRoot);

    // read physical and geometrical simulation parameters
    double spaceChargeFactor = doubleConfParameter("space_charge_factor", confRoot);
    double collisionGasMass_amu = doubleConfParameter("collision_gas_mass_amu", confRoot);
    double collisionGasDiameter_m = doubleConfParameter("collision_gas_diameter_angstrom", confRoot)*1e-10;
    double backgroundGasTemperature_K = doubleConfParameter("background_gas_temperature_K",confRoot);
    double backgroundGasPressure_pa = doubleConfParameter("background_gas_pressure_Pa",confRoot);

    double V_rf = doubleConfParameter("V_rf", confRoot); //volts, RF voltage
    double freq_rf = doubleConfParameter("frequency_rf", confRoot); //Hz, RF Frequency
    double omega_rf = freq_rf * M_PI *2.0;

    // read ion configuration ========================
    std::vector<int> nIons = intVectorConfParameter("n_ions",confRoot);
    std::vector<double> ionMasses = doubleVectorConfParameter("ion_masses",confRoot);
    std::vector<double> ionCharges = doubleVectorConfParameter("ion_charges",confRoot);
    std::vector<double> ionCollisionDiameters_angstrom =
            doubleVectorConfParameter("ion_collision_gas_diameters_angstrom",confRoot);

//    std::vector<double> centerPos = doubleVectorConfParameter("ion_start_box_center_position_m", confRoot);
//    Core::Vector ionStartBoxCenterPosition_m = {centerPos[0], centerPos[1], centerPos[2]};
//    std::vector<double> boxSize = doubleVectorConfParameter("ion_start_box_size_m", confRoot);
//    Core::Vector ionStartBoxSize_m = {boxSize[0], boxSize[1], boxSize[2]};
//    Core::Vector ionStartCornerPosition_m = ionStartBoxCenterPosition_m - (ionStartBoxSize_m*0.5);


    //read potential arrays and potential array configuration =================================================
    // Note that fast adjust PAs are expected here
    std::vector<std::string> potentialArraysNames = stringVectorConfParameter("potential_arrays", confRoot);
    double potentialArrayScale = doubleConfParameter("potential_array_scale", confRoot);
    std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> potentialArrays =
            AppUtils::readPotentialArrayFiles(potentialArraysNames, confBasePath, potentialArrayScale, true);

    //scaling factor of 0.1 because SIMION uses a value of 10000 in  Fast Adjust PAs and  mm to m is 1000 = 0.1
    std::vector<double> potentialsDc = doubleVectorConfParameter("dc_potentials", confRoot);
    std::vector<double> potentialFactorsRf = doubleVectorConfParameter("rf_potential_factors", confRoot);


    // defining simulation domain box (used for ion termination):
    std::array<std::array<double,2>,3> simulationDomainBoundaries;
    if (confRoot.isMember("simulation_domain_boundaries")){
        // get manual simulation domain boundaries from config file
        simulationDomainBoundaries = double3dBox("simulation_domain_boundaries", confRoot);
    } else {
        // TODO: use minimal Potential Array bounds as simulation domain
        throw std::invalid_argument("missing configuration value: simulation_domain_boundaries");
    }

    std::string ionStartGeom_str = stringConfParameter("ion_start_geometry",confRoot);
    IonStartGeometry ionStartGeom;

    double ionStartCylinder_radius;
    double ionStartCylinder_length;

    std::vector<double> centerPos;
    Core::Vector ionStartBoxCenterPosition_m;
    std::vector<double> boxSize;
    Core::Vector ionStartBoxSize_m;
    Core::Vector ionStartCornerPosition_m;

    std::vector<double> vectorStartingPosition;

    if (ionStartGeom_str == "box"){
        ionStartGeom = BOX;
        centerPos = doubleVectorConfParameter("ion_start_box_center_position_m", confRoot);
        ionStartBoxCenterPosition_m = {centerPos[0], centerPos[1], centerPos[2]};
        boxSize = doubleVectorConfParameter("ion_start_box_size_m", confRoot);
        ionStartBoxSize_m = {boxSize[0], boxSize[1], boxSize[2]};
        ionStartCornerPosition_m = ionStartBoxCenterPosition_m - (ionStartBoxSize_m*0.5);

    } else if (ionStartGeom_str == "cylinder"){
        ionStartGeom = CYLINDER;
        ionStartCylinder_radius = doubleConfParameter("ion_start_cylinder_radius_m", confRoot);
        ionStartCylinder_length = doubleConfParameter("ion_start_cylinder_length_m", confRoot);

    } else if (ionStartGeom_str == "vector"){
        ionStartGeom = VECTOR;
        vectorStartingPosition = doubleVectorConfParameter("ion_start_vector_position_m", confRoot);

    }

    //init ions:
    std::vector<std::unique_ptr<BTree::Particle>>particles;
    std::vector<BTree::Particle*>particlePtrs;

    for (int i=0; i<nIons.size(); ++i){
        int nParticles = nIons[i];
        double mass = ionMasses[i];
        double charge = ionCharges[i];
        double collisionDiameter_m = ionCollisionDiameters_angstrom[i]*1e-10;
        std::vector<std::unique_ptr<BTree::Particle>> ions;

        if (ionStartGeom == BOX) {
            ions =  ParticleSimulation::util::getRandomIonsInBox(nParticles, charge, ionStartCornerPosition_m, ionStartBoxSize_m);
        }
        else if (ionStartGeom == CYLINDER){
            ions = ParticleSimulation::util::getRandomIonsInCylinderXDirection(
                    nParticles, charge, ionStartCylinder_radius, ionStartCylinder_length);
        }else if(ionStartGeom == VECTOR){
            ions = ParticleSimulation::util::getIonOnLineVector(nParticles, charge, vectorStartingPosition[0], vectorStartingPosition[1], vectorStartingPosition[2]);
        }


        for (int j=0; j<nParticles; ++j){
            ions[j]->setMassAMU(mass);
            ions[j]->setChargeElementary(charge);
            ions[j]->setDiameter(collisionDiameter_m);

            particlePtrs.push_back(ions[j].get());
            particles.push_back(std::move(ions[j]));
        }
    }


    //init gas collision models:
    CollisionModel::HardSphereModel hsModel = CollisionModel::HardSphereModel(
            backgroundGasPressure_pa, backgroundGasTemperature_K,
            collisionGasMass_amu, collisionGasDiameter_m);

    // define functions for the trajectory integration ==================================================
    auto accelerationFunction =
            [spaceChargeFactor, omega_rf, V_rf, &potentialArrays, &potentialsDc, &potentialFactorsRf](
                    BTree::Particle* particle, int particleIndex, BTree::ParallelTree& tree, double time,
                    int timestep) -> Core::Vector {

                Core::Vector pos = particle->getLocation();
                double particleCharge = particle->getCharge();

                double V_t = cos(omega_rf * time) * V_rf;

                Core::Vector eField(0, 0, 0);
                for (size_t i = 0; i < potentialArrays.size(); i++) {

                    Core::Vector paFieldRaw = potentialArrays[i]->getField(pos.x(), pos.y(), pos.z());
                    Core::Vector paEffectiveField =
                                    paFieldRaw * potentialsDc.at(i) +
                                    paFieldRaw * potentialFactorsRf.at(i) * V_t;

                    eField = eField + paEffectiveField;
                }
                //std::cout << "eField: " << eField << std::endl;
                //std::cout << "Pos: " << pos << std::endl;
                Core::Vector spaceChargeField(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeField = tree.computeEFieldFromTree(*particle) * spaceChargeFactor;
                }
//              std::cout << particle->getVelocity() << std::endl;
                return ((eField + spaceChargeField) * particleCharge / particle->getMass());

            };


    ParticleSimulation::additionalPartParamFctType additionalParameterTransformFct =
            [](BTree::Particle *particle) -> std::vector<double>{
                std::vector<double> result = {
                        particle->getVelocity().x(),
                        particle->getVelocity().y(),
                        particle->getVelocity().z()
                };

                return result;
            };


    //prepare file writers ==============================================================================
    std::vector<std::string> auxParamNames = {"velocity x","velocity y","velocity z"};
    auto hdf5Writer = std::make_unique<ParticleSimulation::TrajectoryHDF5Writer>(
            projectName + "_trajectories.hd5", auxParamNames, additionalParameterTransformFct);

    int ionsInactive = 0;
    auto timestepWriteFunction =
            [trajectoryWriteInterval, &ionsInactive, &additionalParameterTransformFct, &hdf5Writer](
                    std::vector<BTree::Particle*>& particles, auto& tree,
                    double time, int timestep, bool lastTimestep){

                if (lastTimestep) {
                    hdf5Writer->writeTimestep(particles,time);
                    hdf5Writer->writeSplatTimes(particles);
                    hdf5Writer->finalizeTrajectory();
                    std::cout << "finished ts:" << timestep << " time:" << time << std::endl;
                }
                else if (timestep % trajectoryWriteInterval == 0) {

                    std::cout << "ts:" << timestep << " time:" << time << " V_rf:"
                              <<" ions existing: "<<particles.size()<< " ions inactive: "
                              << ionsInactive << std::endl;

                    if (particles.size() > 0) {
                        hdf5Writer->writeTimestep(particles, time);
                    }
                }
            };

    auto otherActionsFunction = [&simulationDomainBoundaries, &ionsInactive, &potentialArrays](
            Core::Vector& newPartPos, BTree::Particle* particle,
            int particleIndex,
            auto& tree, double time, int timestep) {
        // if the ion is out of the boundary box: Terminate the ion
        if ( newPartPos.x() <= simulationDomainBoundaries[0][0] ||
                newPartPos.x() >= simulationDomainBoundaries[0][1] ||
                newPartPos.y() <= simulationDomainBoundaries[1][0] ||
                newPartPos.y() >= simulationDomainBoundaries[1][1] ||
                newPartPos.z() <= simulationDomainBoundaries[2][0] ||
                newPartPos.z() >= simulationDomainBoundaries[2][1] ||
                potentialArrays.at(0)->isElectrode(newPartPos.x(), newPartPos.y(), newPartPos.z() )  )
        {
            std::cout << particle->getVelocity() << std::endl;
            particle->setActive(false);
            particle->setSplatTime(time);
            ionsInactive++;
        }
    };

    // simulate ===============================================================================================
    clock_t begin = std::clock();
    ParticleSimulation::ParallelVerletIntegrator verletIntegrator(
            particlePtrs,
            accelerationFunction, timestepWriteFunction, otherActionsFunction,
            hsModel);
    verletIntegrator.run(timeSteps,dt);

    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    std::cout << particles[0]->getLocation()<<std::endl;
    std::cout << "elapsed secs"<< elapsed_secs<<std::endl;
}