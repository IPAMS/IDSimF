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
 test_MathFunctions.cpp

 Testing of mathematical functions for collision modeling

 ****************************/

#include "catch.hpp"
#include "CollisionModel_MathFunctions.hpp"
#include <array>
#include <cmath>

TEST_CASE("Test trigonometic functions implementation", "[CollisionModels][Math]") {

    SECTION("Radians to degree implementation should be correct") {

        REQUIRE( (CollisionModel::radToDeg(0.0 * M_PI) == Approx(0.0)));
        REQUIRE( (CollisionModel::radToDeg(1.0 * M_PI) == Approx(180.0)));
        REQUIRE( (CollisionModel::radToDeg(0.5 * M_PI) == Approx(90.0)));
        REQUIRE( (CollisionModel::radToDeg(-0.5 * M_PI) == Approx(-90.0)));

        REQUIRE( (CollisionModel::degToRad(0.0) == Approx(0.0 * M_PI)));
        REQUIRE( (CollisionModel::degToRad(180.0) == Approx(1.0 * M_PI)));
        REQUIRE( (CollisionModel::degToRad(90.0) == Approx(0.5 * M_PI)));
        REQUIRE( (CollisionModel::degToRad(-90.0) == Approx(-0.5 * M_PI)));

    }

    SECTION("Cartesian to polar should be correct") {

        SECTION("Test with {1,1,0}"){
            Core::Vector cart = Core::Vector(1.0,1.0,0);
            Core::Vector polar = CollisionModel::cartesianToPolar(cart);
            REQUIRE(polar.x() == Approx(sqrt(2.0)));
            REQUIRE(polar.y() == Approx(0.0));
            REQUIRE(polar.z() == Approx(M_PI_4));
        }

        SECTION("Test with {1,-1,-1}"){
            Core::Vector cart = Core::Vector(1.0,-1.0,-1.0);
            Core::Vector polar  = CollisionModel::cartesianToPolar(cart);
            REQUIRE(polar.x() == Approx(sqrt(3.0)));
            REQUIRE(polar.y() == Approx(M_PI_4));
        }
    }
}

TEST_CASE("Test random functions implementation", "[CollisionModels][Math]") {
    SECTION("Spherical random vectors should have the correct length") {
        std::array<double,5> dists = {1.0,2.0,10.0,100.0,1000.0};

        for(const auto& dist: dists){
            Core::Vector vec = CollisionModel::sphereRand(dist);
            REQUIRE(Approx(vec.magnitude()) == dist);
        }
    }
}
