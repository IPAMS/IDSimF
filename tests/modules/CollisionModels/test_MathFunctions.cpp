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

#include "CollisionModel_MathFunctions.hpp"
#include "catch.hpp"
#include <array>


TEST_CASE("Test random functions implementation", "[CollisionModels][Math]") {
    SECTION("Spherical random vectors should have the correct length") {
        std::array<double,5> dists = {1.0,2.0,10.0,100.0,1000.0};

        for(const auto& dist: dists){
            Core::Vector vec = CollisionModel::sphereRand(dist);
            REQUIRE(Approx(vec.magnitude()) == dist);
        }
    }
}
