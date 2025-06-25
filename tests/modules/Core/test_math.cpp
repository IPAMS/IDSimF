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
 test_math.cpp

 Testing of math functions in core module

 ****************************/


#include "Core_math.hpp"
#include "catch.hpp"
#include "CollisionModel_Atom.hpp"
#include "CollisionModel_MDForceField_Buckingham.hpp"

#include <cmath>
#include <iostream>


TEST_CASE("Test trigonometic functions implementation", "[CollisionModels][Math]") {

    SECTION("Radians to degree implementation should be correct") {

        REQUIRE( (Core::radToDeg(0.0 * M_PI) == Approx(0.0)));
        REQUIRE( (Core::radToDeg(1.0 * M_PI) == Approx(180.0)));
        REQUIRE( (Core::radToDeg(0.5 * M_PI) == Approx(90.0)));
        REQUIRE( (Core::radToDeg(-0.5 * M_PI) == Approx(-90.0)));

        REQUIRE( (Core::degToRad(0.0) == Approx(0.0 * M_PI)));
        REQUIRE( (Core::degToRad(180.0) == Approx(1.0 * M_PI)));
        REQUIRE( (Core::degToRad(90.0) == Approx(0.5 * M_PI)));
        REQUIRE( (Core::degToRad(-90.0) == Approx(-0.5 * M_PI)));
    }

    SECTION("Cartesian to polar should be correct") {

        SECTION("Test with {1,1,0}"){
            Core::Vector cart = Core::Vector(1.0,1.0,0);
            Core::Vector polar = Core::cartesianToPolar(cart);
            REQUIRE(polar.x() == Approx(sqrt(2.0)));
            REQUIRE(polar.y() == Approx(0.0));
            REQUIRE(polar.z() == Approx(M_PI_4));
        }

        SECTION("Test with {1,-1,-1}"){
            Core::Vector cart = Core::Vector(1.0,-1.0,-1.0);
            Core::Vector polar  = Core::cartesianToPolar(cart);
            REQUIRE(polar.x() == Approx(sqrt(3.0)));
            REQUIRE(polar.y() == Approx(M_PI_4));
        }
    }
}

double testFunc(const CollisionModel::Atom& , const CollisionModel::Atom& , double x){ return sin(x);}


TEST_CASE("Test golden section search", "[CollisionModels][Math]") {

    SECTION("Test base case"){

        
        CollisionModel::Atom dummyA;
        dummyA.setSigma(2.5e-10); 
        dummyA.setEpsilon(2.5e-22); 
        CollisionModel::Atom dummyB;
        dummyB.setSigma(2.5e-10);
        dummyB.setEpsilon(2.5e-22); 


        double result = Core::goldenSectionSearch<const CollisionModel::Atom&>(dummyA, dummyB, 
                                                 testFunc, 0, 3, 1e-10);
        double result2 = Core::goldenSectionSearch<const CollisionModel::Atom&>(dummyA, dummyB, 
                                               CollisionModel::MDForceField_Buckingham::calculateVDW, 0.5e-10, 1.5e-10, 1e-22);
                                            
        REQUIRE(result == Approx(1.5708));
        REQUIRE(result2 == Approx(7.55087e-11));

        
    }
}
