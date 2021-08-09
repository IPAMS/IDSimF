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
 test_startZones.cpp

 Testing of generalized particle start zones

 ****************************/

#include "Core_vector.hpp"
#include "PSim_cylinderStartZone.hpp"
#include "PSim_boxStartZone.hpp"
#include "catch.hpp"
#include "test_util.hpp"
#include <vector>


TEST_CASE( "Test box particle start zone",
        "[ParticleSimulation][ParticleStartZone][particle definitions]") {

    int nSamples = 1000;

    SECTION("Test positions in shifted box start zone"){

        Core::Vector center(4.0, 3.0, 1.0);
        Core::Vector boxSize(3.0, 1.0, 2.0);
        Core::Vector lower = center - (boxSize*0.5);
        Core::Vector upper = center + (boxSize*0.5);

        ParticleSimulation::BoxStartZone boxStartZone(boxSize, center);

        bool outOfBox = false;
        for (int i=0; i<nSamples; ++i){
            Core::Vector startPos = boxStartZone.getRandomParticlePosition();
            if (
                startPos.x() < lower.x() || startPos.x() > upper.x() ||
                startPos.y() < lower.y() || startPos.y() > upper.y() ||
                startPos.z() < lower.z() || startPos.z() > upper.z()
            ){
                outOfBox = true;
            }
        }

        REQUIRE_FALSE(outOfBox);
    }
}


TEST_CASE( "Test cylinder particle start zone",
        "[ParticleSimulation][ParticleStartZone][particle definitions]") {

    // define cylinder start zone in x direction
    double radius = 1.0;
    double length = 4.0;
    unsigned int nSamples = 1000;

    SECTION("Check if position samples are out of cylindrical start zone in x-direction"){
        ParticleSimulation::CylinderStartZone cylinderStartZone_x(radius, length);

        bool outOfCylinder = false;
        for (unsigned int i=0; i<nSamples; ++i){
            Core::Vector startPos = cylinderStartZone_x.getRandomParticlePosition();
            double r = std::sqrt(startPos.y()*startPos.y() + startPos.z()*startPos.z());

            if (startPos.x() > length || startPos.x() < 0.0 || r > radius){
                outOfCylinder = true;
            }
        }
        REQUIRE_FALSE(outOfCylinder);
    }

    SECTION("Check if position samples are out of shifted cylindrical start zone in z-direction"){
        ParticleSimulation::CylinderStartZone cylinderStartZone_z(radius, length, {0.0, 0.0, 3.0});

        bool outOfCylinder = false;
        for (unsigned int i=0; i<nSamples; ++i){
            Core::Vector startPos = cylinderStartZone_z.getRandomParticlePosition();
            double r = std::sqrt(startPos.x()*startPos.x() + startPos.y()*startPos.y());

            if (startPos.z() > length || startPos.z() < 0.0 || r > radius){
                outOfCylinder = true;
            }
        }
        REQUIRE_FALSE(outOfCylinder);
    }

    SECTION("Check if position samples are out of shifted cylindrical start zone in y-direction with shift"){

        Core::Vector shift(1.0, 2.0, 0.0);
        ParticleSimulation::CylinderStartZone cylinderStartZone_y(radius, length, {0.0, 2.0, 0.0}, shift);

        bool outOfCylinder = false;
        for (unsigned int i=0; i<nSamples; ++i){
            Core::Vector startPos = cylinderStartZone_y.getRandomParticlePosition();
            startPos = startPos - shift;

            double r = std::sqrt(startPos.x()*startPos.x() + startPos.z()*startPos.z());

            if (startPos.y() > length || startPos.y() < 0.0 || r > radius){
                outOfCylinder = true;
            }
        }
        REQUIRE_FALSE(outOfCylinder);
    }

    SECTION("Check if two random position samples are different"){
        ParticleSimulation::CylinderStartZone cylinderStartZone(radius, length);
        CHECK(cylinderStartZone.getRandomParticlePosition() != cylinderStartZone.getRandomParticlePosition());
    }

    SECTION("Test generation of particle cloud in cylindrical start zone"){

        Core::Vector shift(2.0, 1.0, 0.0);
        ParticleSimulation::CylinderStartZone cylinderStartZone(radius, length, {0.0, 2.0, 0.0}, shift);

        std::vector<std::unique_ptr<BTree::Particle>> ions =
                cylinderStartZone.getRandomParticlesInStartZone(nSamples, 2.0, 5.0);


        double chargeExpected = 2.0 * Core::ELEMENTARY_CHARGE;
        bool incorrectIonFound = false;
        for (std::size_t i=0; i<nSamples; ++i){
            Core::Vector pos = ions[i]->getLocation();
            pos = pos - shift;
            double r = std::sqrt(pos.x()*pos.x() + pos.z()*pos.z());

            double charge = ions[i]->getCharge();
            double tob = ions[i]->getTimeOfBirth();

            if (pos.y() > length || pos.y() < 0.0 || r > radius ||
                !isExactDoubleEqual(charge, chargeExpected) || tob < 0.0 || tob > 5.0
            ){
                incorrectIonFound = true;
            }
        }
        REQUIRE_FALSE(incorrectIonFound);

        CHECK(ions[0]->getLocation() != ions[1]->getLocation());
    }

}
