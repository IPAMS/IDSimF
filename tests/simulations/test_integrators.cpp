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
 test_integrators.cpp

 More complex and comparative testing of verlet integrators implementations

 ****************************/

#include "catch.hpp"
#include "BTree_parallelNodeOriginal.hpp"
#include "BTree_tree.hpp"
#include "BTree_parallelTree.hpp"
#include "BTree_parallelTreeOriginal.hpp"
#include "BTree_particle.hpp"
#include "PSim_verletIntegrator.hpp"
#include "PSim_parallelVerletIntegrator.hpp"
#include "PSim_util.hpp"
#include "CollisionModel_EmptyCollisionModel.hpp"
#include <iostream>
#include <cmath>
#include <numeric>

void prepareIons(std::vector<std::unique_ptr<BTree::Particle>> &particles,
                 std::vector<BTree::Particle*> &particlePtrs, int nIons){

    for (int i=0; i<nIons; i++){
        double posy = i*1.0/nIons;
        std::unique_ptr<BTree::Particle> newIon = std::make_unique<BTree::Particle>(Core::Vector(0,posy,0), 1.0);
        newIon -> setMassAMU(100);
        particlePtrs.push_back(newIon.get());
        particles.push_back(std::move(newIon));
    }
}

TEST_CASE("Compare results of serial and parallel varlet integrators with a line of charged particles", "[Simulation]"){

    int nIons = 200;
    int timeSteps = 1000;
    double dt = 3.0e-4;
    double spaceChargeFactor = 1.0;

    // define functions for the trajectory integration ==================================================
    auto accelerationFunctionSerial =
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


    auto accelerationFunctionParallel =
            [spaceChargeFactor](
                    BTree::Particle *particle, int particleIndex,
                    BTree::ParallelTreeOriginal &tree, double time, int timestep,std::vector<BTree::ParallelNodeOriginal*> &MyNodes) -> Core::Vector{

                Core::Vector pos = particle->getLocation();
                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle,MyNodes) * (particleCharge * spaceChargeFactor);
                }
                return (spaceChargeForce / particle->getMass());
            };

    auto accelerationFunctionParallelNew =
            [spaceChargeFactor](
                    BTree::Particle *particle, int particleIndex,
                    BTree::ParallelTree &tree, double time, int timestep) -> Core::Vector{

                Core::Vector pos = particle->getLocation();
                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle) * (particleCharge * spaceChargeFactor);
                }
                return (spaceChargeForce / particle->getMass());
            };


    std::vector<std::unique_ptr<BTree::Particle>> particlesSerial;
    std::vector<BTree::Particle*>particlePtrsSerial;
    std::vector<std::unique_ptr<BTree::Particle>> particlesParallelNew;
    std::vector<BTree::Particle*>particlePtrsParallelNew;

    prepareIons(particlesSerial, particlePtrsSerial, nIons);
    prepareIons(particlesParallelNew, particlePtrsParallelNew, nIons);


    // simulate ===============================================================================================
    ParticleSimulation::VerletIntegrator verletIntegratorSerial(
            particlePtrsSerial, accelerationFunctionSerial);

    ParticleSimulation::ParallelVerletIntegrator verletIntegratorParallelNew(
            particlePtrsParallelNew, accelerationFunctionParallelNew);

    verletIntegratorSerial.run(timeSteps, dt);
    verletIntegratorParallelNew.run(timeSteps, dt);

    std::vector<double> diffMags;
    for (int i=0; i<nIons; ++i){
        diffMags.push_back( (particlesSerial[i]->getLocation() - particlesParallelNew[i]->getLocation()).magnitude() );
    }
    double sum = std::accumulate(diffMags.begin(), diffMags.end(), 0.0);
    double maximumDiff = *std::max_element(diffMags.begin(), diffMags.end());

    for (int i=0; i<nIons; ++i){
        std::cout <<
            particlesSerial[i]->getLocation()<< " | " <<
            particlesParallelNew[i]->getLocation()<< " | " <<
            (particlesSerial[i]->getLocation() - particlesParallelNew[i]->getLocation()).magnitude()
        << std::endl;
    }

    REQUIRE(sum <= 1e-12);
    REQUIRE(maximumDiff <= 1e-14);
}