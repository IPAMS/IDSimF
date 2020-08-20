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

 Testing of math functions

 ****************************/

#include "PSim_math.hpp"
#include <iostream>
#include <limits>
#include "vector"
#include "catch.hpp"

bool compareVectors(std::vector<double> a, std::vector<double> b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] != Approx(b[i])) {
            std::cout << a[i] << " Should == " << b[i] << std::endl;
            return false;
        }
    }
    return true;
}

TEST_CASE("Test linspace function", "[ParticleSimulation][math]") {

    SECTION("Test linspace gives correct length") {
        std::vector<double> a = ParticleSimulation::linspace(1.0, 2.0, 10);
        REQUIRE(a.size() == 10);
    }

    SECTION("Test linspace gives correct vector with upper > lower") {
        std::vector<double> a = ParticleSimulation::linspace(2.0, 4.0, 5);
        std::vector<double> b = {2.0, 2.5, 3.0, 3.5, 4.0 };
        REQUIRE(compareVectors(a,b));
    }

    SECTION( "Test linspace gives correct vector with upper < lower") {
        std::vector<double> a = ParticleSimulation::linspace(4.0, 2.0, 5);
        std::vector<double> b = {4.0, 3.5, 3.0, 2.5, 2.0 };
        REQUIRE(compareVectors(a,b));
    }

    SECTION( "Test linspace gives correct vector if boundaries are equal") {
        std::vector<double> a = ParticleSimulation::linspace(2.0, 2.0, 5);
        std::vector<double> b = {2.0, 2.0, 2.0, 2.0, 2.0 };
        REQUIRE(compareVectors(a,b));
    }

}

TEST_CASE ("Check max_int value","[ParticleSimulation][math]"){
    std::cout <<"INT_MAX:"<< std::numeric_limits<int>::max() <<std::endl;
    std::cout <<"LONG_MAX:"<< std::numeric_limits<long>::max() <<std::endl;
}
