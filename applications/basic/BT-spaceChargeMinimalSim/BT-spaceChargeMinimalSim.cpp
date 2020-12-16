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
 BT-spaceChargeMinimalSim.cpp

 Minimal serial simulation of pure particle / particle interaction (space charge)

 ****************************/

#include "json.h"
#include "appUtils_parameterParsing.hpp"
#include "BTree_particle.hpp"
#include "BTree_tree.hpp"
#include "PSim_trajectoryExplorerJSONwriter.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_verletIntegrator.hpp"
#include "PSim_ionCloudReader.hpp"
#include "CollisionModel_EmptyCollisionModel.hpp"
#include <filesystem>
#include <iostream>
#include <vector>
#include <ctime>

int main(int argc, const char * argv[]) {

    // read configuration file ======================================================================
    if (argc <2){
        std::cout << "no conf project name or conf file given"<<std::endl;
        return(1);
    }

    std::string confFileName = argv[1];
    Json::Value confRoot = readConfigurationJson(confFileName);

    std::string projectName = argv[2];
    std::cout << projectName<<std::endl;


    // read basic simulation parameters =============================================================
    int timeSteps = intConfParameter("sim_time_steps", confRoot);
    int trajectoryWriteInterval = intConfParameter("trajectory_write_interval", confRoot);
    double dt = doubleConfParameter("dt", confRoot);

    //read physical configuration ===================================================================
    double spaceChargeFactor = doubleConfParameter("space_charge_factor", confRoot);


    //read ion configuration =======================================================================
    std::vector<std::unique_ptr<BTree::Particle>>particles;
    std::vector<BTree::Particle*>particlePtrs;

    if (confRoot.isMember("ion_cloud_init_file") == true) {
        std::string ionCloudFileName = pathRelativeToConfFile(
                                        confFileName,
                                        confRoot.get("ion_cloud_init_file", 0).asString());
        ParticleSimulation::IonCloudReader reader = ParticleSimulation::IonCloudReader();
        particles = reader.readIonCloud(ionCloudFileName);
        //prepare a vector of raw pointers
        for (const auto& part : particles){
            particlePtrs.push_back(part.get());
        }
    } else {
        throw std::invalid_argument("missing configuration value: ion_cloud_init_file");
    }

    //prepare file writer ==============================================================================

    // function to add some additional exported parameters to the exported trajectory file:
    ParticleSimulation::additionalPartParamFctType additionalParameterTransformFct =
            [](BTree::Particle *particle) -> std::vector<double>{
                std::vector<double> result = {
                        particle->getVelocity().x(),
                        particle->getVelocity().y(),
                        particle->getVelocity().z()
                };
                return result;
            };

    std::vector<std::string> auxParamNames = {"velocity x","velocity y","velocity z"};

    auto hdf5Writer = std::make_unique<ParticleSimulation::TrajectoryHDF5Writer>(
            projectName + "_trajectories.hd5", auxParamNames,additionalParameterTransformFct);

    /*auto jsonWriter = std::make_unique<ParticleSimulation::TrajectoryExplorerJSONwriter>(
            projectName + "_trajectories.json");*/
    //jsonWriter->setScales(1000,1e6);


    // define functions for the trajectory integration ==================================================

    auto accelerationFunction =
            [spaceChargeFactor](
                    BTree::Particle *particle, int particleIndex,
                    BTree::Tree &tree, double time, int timestep) -> Core::Vector{

                Core::Vector pos = particle->getLocation();
                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle) * (particleCharge * spaceChargeFactor);
                }
                return (spaceChargeForce / particle->getMass());
            };

    auto timestepWriteFunction =
            [trajectoryWriteInterval, &hdf5Writer, &additionalParameterTransformFct](
                    std::vector<BTree::Particle *> &particles, BTree::Tree &tree, double time,
                    int timestep, bool lastTimestep){

                if (lastTimestep) {
                    hdf5Writer->writeTimestep(particles,time);

                    hdf5Writer->writeSplatTimes(particles);
                    hdf5Writer->finalizeTrajectory();
                    std::cout << "finished ts:" << timestep << " time:" << time << std::endl;
                }

                else if (timestep % trajectoryWriteInterval == 0) {

                    std::cout << "ts:" << timestep << " time:" << time << std::endl;
                    hdf5Writer->writeTimestep(particles,time);
                }
    };


    //an empty other actions function (to do nothing additionally in a timestep)
    auto otherActionsFunctionQIT = [](Core::Vector &newPartPos, BTree::Particle *particle,
                                                                 int particleIndex,
                                                                 BTree::Tree &tree, double time, int timestep){
    };

    CollisionModel::EmptyCollisionModel emptyCollisionModel;


    // simulate ===============================================================================================
    clock_t begin = std::clock();
    ParticleSimulation::VerletIntegrator verletIntegrator(
            particlePtrs,
            accelerationFunction, timestepWriteFunction, otherActionsFunctionQIT,
            emptyCollisionModel);
    verletIntegrator.run(timeSteps, dt);

    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    std::cout << particles[0]->getLocation()<<std::endl;
    std::cout << "elapsed secs"<< elapsed_secs<<std::endl;
    return 0;
}