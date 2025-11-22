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
#include "Integration_fullSumVerletIntegrator.hpp"
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

TEST_CASE("Compare results of serial and parallel varlet integrators and full sum with a line of charged particles", "[Simulation]"){

    unsigned int nIons = 200;
    unsigned int timeSteps = 1000;
    double dt = 3.0e-4;
    double spaceChargeFactor = 1.0;

    // define functions for the trajectory integration ==================================================
    auto accelerationFunction =
            [spaceChargeFactor](
                    Core::Particle *particle, int /*particleIndex*/,
                    SpaceCharge::FieldCalculator &fieldCalculator, double /*time*/, int /*timestep*/) -> Core::Vector{

                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            fieldCalculator.getEFieldFromSpaceCharge(*particle) * (particleCharge * spaceChargeFactor);
                }
                return (spaceChargeForce / particle->getMass());
            };

    std::vector<std::unique_ptr<Core::Particle>> particlesSerial;
    std::vector<Core::Particle*>particlePtrsSerial;
    std::vector<std::unique_ptr<Core::Particle>> particlesParallelNew;
    std::vector<Core::Particle*>particlePtrsParallelNew;
    std::vector<std::unique_ptr<Core::Particle>> particlesFullSum;
    std::vector<Core::Particle*>particlePtrsFullSum;


    prepareIons(particlesSerial, particlePtrsSerial, nIons);
    prepareIons(particlesParallelNew, particlePtrsParallelNew, nIons);
    prepareIons(particlesFullSum, particlePtrsFullSum, nIons);


    // simulate ===============================================================================================
    Integration::VerletIntegrator verletIntegratorSerial(
            particlePtrsSerial, accelerationFunction);

    Integration::ParallelVerletIntegrator verletIntegratorParallelNew(
            particlePtrsParallelNew, accelerationFunction);

    Integration::FullSumVerletIntegrator verletIntegratorFullSum(
            particlePtrsFullSum, accelerationFunction);


    verletIntegratorSerial.run(timeSteps, dt);
    verletIntegratorParallelNew.run(timeSteps, dt);
    verletIntegratorFullSum.run(timeSteps, dt);

    std::vector<double> diffMags_s_p;
    std::vector<double> diffMags_s_fs;
    std::vector<double> diffMags_p_fs;

    for (unsigned int i=0; i<nIons; ++i){
        diffMags_s_p.push_back( (particlesSerial[i]->getLocation() - particlesParallelNew[i]->getLocation()).magnitude() );
        diffMags_s_fs.push_back( (particlesSerial[i]->getLocation() - particlesFullSum[i]->getLocation()).magnitude() );
        diffMags_p_fs.push_back( (particlesParallelNew[i]->getLocation() - particlesFullSum[i]->getLocation()).magnitude() );
    }
    double sum_s_p = std::accumulate(diffMags_s_p.begin(), diffMags_s_p.end(), 0.0);
    double sum_s_fs = std::accumulate(diffMags_s_fs.begin(), diffMags_s_fs.end(), 0.0);
    double sum_p_fs = std::accumulate(diffMags_p_fs.begin(), diffMags_p_fs.end(), 0.0);
    double maximumDiff_s_p = *std::max_element(diffMags_s_p.begin(), diffMags_s_p.end());
    double maximumDiff_s_fs = *std::max_element(diffMags_s_fs.begin(), diffMags_s_fs.end());
    double maximumDiff_p_fs = *std::max_element(diffMags_p_fs.begin(), diffMags_p_fs.end());

    /*for (unsigned int i=0; i<nIons; ++i){
        std::cout <<
            particlesSerial[i]->getLocation()<< " | " <<
            particlesParallelNew[i]->getLocation()<< " | " <<
            (particlesSerial[i]->getLocation() - particlesParallelNew[i]->getLocation()).magnitude()
        << std::endl;
    }*/

    CHECK(sum_s_p <= 1e-12);
    CHECK(sum_s_fs <= 1e-12);
    CHECK(sum_p_fs <= 1e-12);
    CHECK(maximumDiff_s_p <= 1e-14);
    CHECK(maximumDiff_s_fs <= 1e-14);
    CHECK(maximumDiff_p_fs <= 1e-14);
}