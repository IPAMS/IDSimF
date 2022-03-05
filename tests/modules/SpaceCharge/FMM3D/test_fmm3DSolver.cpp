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
#include "catch.hpp"
#include "test_util.hpp"

TEST_CASE( "Test basic particle/particle interaction calculation with FMM3D", "[FMM3D]"){
    std::size_t nions = 10000;
    Core::Vector boxSize(0.002, 0.002, 0.002);
    ParticleSimulation::BoxStartZone startZone(boxSize);
    std::vector<std::unique_ptr<Core::Particle>> ions= startZone.getRandomParticlesInStartZone(nions, 1);

    FMM3D::FMMSolver fmmSolver;

    for (std::size_t i=0; i<nions; i++){
        fmmSolver.insertParticle((*ions[i]),i+1);
    }

    fmmSolver.computeChargeDistribution();

    Core::Vector spaceChargeForce = fmmSolver.getEFieldFromSpaceCharge(*ions[0].get());
    CHECK(spaceChargeForce.x() != 0.0);
    CHECK(spaceChargeForce.y() != 0.0);
    CHECK(spaceChargeForce.z() != 0.0);
}