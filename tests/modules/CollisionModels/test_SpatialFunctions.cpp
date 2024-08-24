/***************************
Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2024 - Physical and Theoretical Chemistry /
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
 test_SpatialFunctions.cpp

 Testing of spatial functions / fields

 ****************************/

#include "CollisionModel_SpatialFieldFunctions.hpp"
#include "catch.hpp"

TEST_CASE("Test spatial functions", "[CollisionModels][SpatialFunctions]") {
    SECTION("Constant spatial functions should produce constant values") {
        auto constScalarFunction = CollisionModel::getConstantDoubleFunction(10.0);
        REQUIRE(constScalarFunction({1.0, 2.0, 3.0}) == Approx(10.0));
        REQUIRE(constScalarFunction({-2.0, -4.0, -6.0}) == Approx(10.0));
    }
}
