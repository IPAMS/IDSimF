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
 test_treeParticle.cpp

 Testing of simulated particles wrapped for BTree

 ****************************/

#include "Core_particle.hpp"
#include "BTree_treeParticle.hpp"
#include "catch.hpp"
#include "test_util.hpp"

TEST_CASE( "BTree Particle tests", "[Particle]") {

    Core::Particle baseParticle = Core::Particle(Core::Vector(1.0, 2.0, 3.0), 2.0);
    BTree::TreeParticle testParticle(&baseParticle);

    SECTION("Test for Host Node Management") {
        //Test ion location:
        CHECK(testParticle.getHostNode()==nullptr);
    }
}