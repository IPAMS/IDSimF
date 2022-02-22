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
 test_vector.cpp

 Testing of basic vector class

 ****************************/

#include "Core_vector.hpp"
#include <cmath>
#include "catch.hpp"
#include "test_util.hpp"

TEST_CASE("Test Vector construction","[Core][Vector]"){
    SECTION("Vector construction with explicit parameters is working"){
        Core::Vector vec(1.0,2.0,3.0);
        REQUIRE(isExactDoubleEqual(vec.x(), 1.0));
        REQUIRE(isExactDoubleEqual(vec.y(), 2.0));
        REQUIRE(isExactDoubleEqual(vec.z(), 3.0));
    }

    SECTION("Vector construction with array as initializer is working"){
        double arr[3] ={5.0,6.0,7.0};
        Core::Vector vec(arr);
        REQUIRE(isExactDoubleEqual(vec.x(), 5.0));
        REQUIRE(isExactDoubleEqual(vec.y(), 6.0));
        REQUIRE(isExactDoubleEqual(vec.z(), 7.0));
    }
}

TEST_CASE("Test Vector operations","[Core][Vector]") {

    Core::Vector a(1.0, 2.5, 5.0);
    Core::Vector b(0.5, 1.25, 2.5);
    Core::Vector d(1.0, 1.0, 1.0);
    Core::Vector e(1.0, 0, 0);
    double c = 10.0;

    SECTION("Basic Vector operations are working") {
        REQUIRE((a+b)==Core::Vector(1.5, 3.75, 7.5));
        REQUIRE((a-b)==Core::Vector(0.5, 1.25, 2.5));
        REQUIRE(isExactDoubleEqual((a*a), 32.25));
        REQUIRE((a*c)==Core::Vector(10.0, 25.0, 50.0));

        REQUIRE(d.magnitude()== Approx(sqrt(3.0)));
        REQUIRE(e.magnitude()== Approx(1.0));
    }

    SECTION("Vector accessors are working") {
        REQUIRE(isExactDoubleEqual(a.x(), 1.0));
        REQUIRE(isExactDoubleEqual(a.y(), 2.5));
        REQUIRE(isExactDoubleEqual(a.z(), 5.0));
    }

    SECTION("Vector accessors return copies of the values") {
        Core::Vector d = Core::Vector(1.0, 2.5, 5.0);
        double x = d.x();
        REQUIRE(isExactDoubleEqual(x, 1.0));
        x = 2.0;
        REQUIRE(isExactDoubleEqual(a.x(), 1.0));
    }

    SECTION("Vector setter is working") {
        Core::Vector f = Core::Vector(1.0, 2.5, 5.0);
        f.set(2.0, 3.0, 4.0);
        REQUIRE(f==Core::Vector(2.0, 3.0, 4.0));
    }
}

