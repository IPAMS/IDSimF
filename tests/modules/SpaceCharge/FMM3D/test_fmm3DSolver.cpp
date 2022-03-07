/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2022 - Physical and Theoretical Chemistry /
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
 test_fmmSolver.cpp

 Basic testing of FMM3D based FMM solver

 ****************************/
#include "FMM3D_fmmSolver.hpp"
#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include "PSim_boxStartZone.hpp"
#include "SC_fullSumSolver.hpp"
#include "catch.hpp"
#include "test_util.hpp"

TEST_CASE( "Test basic particle/particle interaction calculation with FMM3D", "[FMM3D]"){

    FMM3D::FMMSolver fmm3dSolver;
    SpaceCharge::FullSumSolver fullSumSolver;

    SECTION("Force in random box should be approximately the same as with full sum"){
        std::size_t nIons = 10000;
        Core::Vector boxSize(0.002, 0.002, 0.002);
        ParticleSimulation::BoxStartZone startZone(boxSize);
        std::vector<std::unique_ptr<Core::Particle>> ions= startZone.getRandomParticlesInStartZone(nIons, 1);

        for (std::size_t i=0; i<nIons; i++){
            fmm3dSolver.insertParticle((*ions[i]), i+1);
            fullSumSolver.insertParticle((*ions[i]), i+1);
        }

        CHECK(fmm3dSolver.getNumberOfParticles() == nIons);
        CHECK(fullSumSolver.getNumberOfParticles() == nIons);

        fmm3dSolver.computeChargeDistribution();

        Core::Vector forceExa0 = fmm3dSolver.getEFieldFromSpaceCharge(*ions[0].get());
        Core::Vector forceFullSum0 = fullSumSolver.getEFieldFromSpaceCharge(*ions[0].get());

        CHECK( ((forceExa0-forceFullSum0).magnitude() / forceExa0.magnitude()) < 1e-6);
    }
}