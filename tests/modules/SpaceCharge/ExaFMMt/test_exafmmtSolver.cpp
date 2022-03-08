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
 test_exafmmtSolver.cpp

 Basic testing of ExaFMM-t based FMM solver

 ****************************/
#include "ExaFMMt_fmmSolver.hpp"
#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include "PSim_boxStartZone.hpp"
#include "SC_fullSumSolver.hpp"
#include "catch.hpp"
#include "test_util.hpp"
#include "test_particleInit.hpp"

#include <iostream>

TEST_CASE( "Test particle/particle interaction calculation with ExaFMMt", "[ExaFMMt]"){

    ExaFMMt::FMMSolver exafmmSolver;
    SpaceCharge::FullSumSolver fullSumSolver;

    SECTION("Force in random box should be approximately the same as with full sum"){
        std::size_t nIons = 5000;
        Core::Vector boxSize(0.002, 0.002, 0.002);
        ParticleSimulation::BoxStartZone startZone(boxSize);
        std::vector<std::unique_ptr<Core::Particle>> ions= startZone.getRandomParticlesInStartZone(nIons, 1);

        for (std::size_t i=0; i<nIons; i++){
            exafmmSolver.insertParticle((*ions[i]), i+1);
            fullSumSolver.insertParticle((*ions[i]), i+1);
        }

        CHECK(exafmmSolver.getNumberOfParticles() == nIons);
        CHECK(fullSumSolver.getNumberOfParticles() == nIons);

        exafmmSolver.computeChargeDistribution();

        Core::Vector force1 = exafmmSolver.getEFieldFromSpaceCharge(*ions[1]);
        Core::Vector fullSumForce1 = fullSumSolver.getEFieldFromSpaceCharge(*ions[1]);
        CHECK(vectorApproxCompare(force1, fullSumForce1) == vectorsApproxEqual);

        Core::Vector force10 = exafmmSolver.getEFieldFromSpaceCharge(*ions[10]);
        Core::Vector fullSumForce10 = fullSumSolver.getEFieldFromSpaceCharge(*ions[10]);
        CHECK(vectorApproxCompare(force10, fullSumForce10) == vectorsApproxEqual);
    }

    SECTION( "Test force calculation with a large number of particles in a latticed cube"){
        std::size_t nPerDirection = 15;
        auto ions = getIonsInLattice(nPerDirection);

        std::size_t i = 0;
        for (auto& ion: ions){
            fullSumSolver.insertParticle(*ion, i);
            exafmmSolver.insertParticle(*ion, i);
            ++i;
        }

        exafmmSolver.computeChargeDistribution();
        CHECK(exafmmSolver.getNumberOfParticles() == nPerDirection * nPerDirection * nPerDirection);

        Core::Vector force1 = exafmmSolver.getEFieldFromSpaceCharge(*ions[1]);
        Core::Vector fullSumForce1 = fullSumSolver.getEFieldFromSpaceCharge(*ions[1]);
        CHECK(vectorApproxCompare(force1, fullSumForce1) == vectorsApproxEqual);

        Core::Vector force10 = exafmmSolver.getEFieldFromSpaceCharge(*ions[10]);
        Core::Vector fullSumForce10 = fullSumSolver.getEFieldFromSpaceCharge(*ions[10]);
        CHECK(vectorApproxCompare(force10, fullSumForce10) == vectorsApproxEqual);

        Core::Vector force100 = exafmmSolver.getEFieldFromSpaceCharge(*ions[100]);
        Core::Vector fullSumForce100 = fullSumSolver.getEFieldFromSpaceCharge(*ions[100]);
        CHECK(vectorApproxCompare(force100, fullSumForce100) == vectorsApproxEqual);
    }

    SECTION( "Test force calculation with small cube of charges"){

        std::vector<Core::Particle> particles = {
            Core::Particle(Core::Vector( 0.1, 0.1, 0.1), 1.0),
            Core::Particle(Core::Vector( 0.1,-0.1, 0.1), 1.0),
            Core::Particle(Core::Vector(-0.1, 0.1, 0.1), 1.0),
            Core::Particle(Core::Vector(-0.1,-0.1, 0.1), 1.0),
            Core::Particle(Core::Vector( 0.1, 0.1,-0.1), 1.0),
            Core::Particle(Core::Vector( 0.1,-0.1,-0.1), 1.0),
            Core::Particle(Core::Vector(-0.1, 0.1,-0.1), 1.0),
            Core::Particle(Core::Vector(-0.1,-0.1,-0.1), 1.0),
            Core::Particle(Core::Vector( 0.0, 0.0, 0.2), 1.0)
        };

        std::size_t i=0;
        for (auto& particle: particles){
            exafmmSolver.insertParticle(particle, i);
            fullSumSolver.insertParticle(particle, i);
            ++i;
        }

        exafmmSolver.computeChargeDistribution();

        Core::Vector force1 = exafmmSolver.getEFieldFromSpaceCharge(particles[8]);
        Core::Vector fullSumForce1 = fullSumSolver.getEFieldFromSpaceCharge(particles[8]);
        CHECK(vectorApproxCompare(force1, fullSumForce1) == vectorsApproxEqual);
    }
}