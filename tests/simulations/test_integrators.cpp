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

#include "BTree_tree.hpp"
#include "BTree_parallelTree.hpp"
#include "Core_particle.hpp"
#include "Integration_verletIntegrator.hpp"
#include "Integration_parallelVerletIntegrator.hpp"
#include "PSim_util.hpp"
#include "catch.hpp"
#include <iostream>
#include <cmath>
#include <numeric>

void prepareIons(std::vector<std::unique_ptr<Core::Particle>> &particles,
                 std::vector<Core::Particle*> &particlePtrs, unsigned int nIons){

    for (unsigned int i=0; i < nIons; ++i){
        double posy = i*1.0/nIons;
        std::unique_ptr<Core::Particle> newIon = std::make_unique<Core::Particle>(Core::Vector(0,posy,0), 1.0);
        newIon -> setMassAMU(100);
        particlePtrs.push_back(newIon.get());
        particles.push_back(std::move(newIon));
    }
}

TEST_CASE("Compare results of serial and parallel varlet integrators with a line of charged particles", "[Simulation]"){

    unsigned int nIons = 200;
    unsigned int timeSteps = 1000;
    double dt = 3.0e-4;
    double spaceChargeFactor = 1.0;

    // define functions for the trajectory integration ==================================================
    auto accelerationFunctionSerial =
            [spaceChargeFactor](
                    Core::Particle *particle, int /*particleIndex*/,
                    BTree::Tree &tree, double /*time*/, int /*timestep*/) -> Core::Vector{

                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle) * (particleCharge * spaceChargeFactor);
                }
                return (spaceChargeForce / particle->getMass());
            };

    auto accelerationFunctionParallelNew =
            [spaceChargeFactor](
                    Core::Particle *particle, int /*particleIndex*/,
                    BTree::ParallelTree &tree, double /*time*/, int /*timestep*/) -> Core::Vector{

                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle) * (particleCharge * spaceChargeFactor);
                }
                return (spaceChargeForce / particle->getMass());
            };


    std::vector<std::unique_ptr<Core::Particle>> particlesSerial;
    std::vector<Core::Particle*>particlePtrsSerial;
    std::vector<std::unique_ptr<Core::Particle>> particlesParallelNew;
    std::vector<Core::Particle*>particlePtrsParallelNew;

    prepareIons(particlesSerial, particlePtrsSerial, nIons);
    prepareIons(particlesParallelNew, particlePtrsParallelNew, nIons);


    // simulate ===============================================================================================
    Integration::VerletIntegrator verletIntegratorSerial(
            particlePtrsSerial, accelerationFunctionSerial);

    Integration::ParallelVerletIntegrator verletIntegratorParallelNew(
            particlePtrsParallelNew, accelerationFunctionParallelNew);

    verletIntegratorSerial.run(timeSteps, dt);
    verletIntegratorParallelNew.run(timeSteps, dt);

    std::vector<double> diffMags;
    for (unsigned int i=0; i<nIons; ++i){
        diffMags.push_back( (particlesSerial[i]->getLocation() - particlesParallelNew[i]->getLocation()).magnitude() );
    }
    double sum = std::accumulate(diffMags.begin(), diffMags.end(), 0.0);
    double maximumDiff = *std::max_element(diffMags.begin(), diffMags.end());

    for (unsigned int i=0; i<nIons; ++i){
        std::cout <<
            particlesSerial[i]->getLocation()<< " | " <<
            particlesParallelNew[i]->getLocation()<< " | " <<
            (particlesSerial[i]->getLocation() - particlesParallelNew[i]->getLocation()).magnitude()
        << std::endl;
    }

    CHECK(sum <= 1e-12);
    CHECK(maximumDiff <= 1e-14);
}