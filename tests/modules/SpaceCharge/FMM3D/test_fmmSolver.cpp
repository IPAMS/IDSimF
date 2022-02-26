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
 test_fmmSolver.cpp

 Basic testing of FMM3D based FMM solver

 ****************************/
#include "FMM3D_fmmSolver.hpp"
#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include "catch.hpp"
#include "test_util.hpp"

TEST_CASE( "Test basic particle/particle interaction calculation with FMM3D", "[FMM3D]"){
    FMM3D::FMMSolver fmmSolver;

    Core::Particle testParticle;
    Core::Vector spaceChargeForce = fmmSolver.computeEFieldFromSpaceCharge(testParticle);
    CHECK(spaceChargeForce.x() == 1.0);
}