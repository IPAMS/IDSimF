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
 test_sampledFunction.cpp

 Test of sampled, discrete function

 ****************************/

#include "PSim_sampledFunction.hpp"
#include "catch.hpp"

TEST_CASE("Test sampled function", "[ParticleSimulation][SampledFunction][file readers]") {

    SECTION("Function data files should be readable and correct") {
        ParticleSimulation::SampledFunction sf_linear("function_linear.csv");
        ParticleSimulation::SampledFunction sf_quadratic("function_quadratic.csv");

        SECTION("State should be good and raw data should be correct") {
            CHECK(sf_linear.good());
            CHECK(sf_linear[2].first== Approx(10.0));
            CHECK(sf_linear[2].second== Approx(1.0));
            CHECK(sf_linear.getIndependentValue(3)==Approx(50.0));
            CHECK(sf_linear.getFunctionValue(3)==Approx(1000.0));

            CHECK(sf_quadratic.good());
            CHECK(sf_quadratic[5].first== Approx(3.3557046979866e-01));
            CHECK(sf_quadratic[5].second== Approx(1.1260754020089e-01));
            CHECK(sf_quadratic.getIndependentValue(6)==Approx(4.02684563758389e-01));
            CHECK(sf_quadratic.getFunctionValue(6)==Approx(1.621548578893e-01));

            //CHECK(sf_quadratic.getFunctionValue(3)==Approx(1000.0));
        }

        SECTION("Read from an non existing index should throw an exception") {
            CHECK(sf_linear.good());
            CHECK_THROWS(sf_linear[100000000]);
        }

        SECTION("Interpolated values should be correct") {
            CHECK(sf_linear.getInterpolatedValue(0.0) == Approx(1.0));
            CHECK(sf_linear.getInterpolatedValue(2.0) == Approx(1.5));
            CHECK(sf_linear.getInterpolatedValue(7.0) == Approx(1.5));
            CHECK(sf_linear.getInterpolatedValue(40.0) == Approx(750.25));

            CHECK(sf_quadratic.getInterpolatedValue(0.0) == Approx(0.0));
            CHECK(sf_quadratic.getInterpolatedValue(5.5) == Approx(5.5*5.5));
            CHECK(sf_quadratic.getInterpolatedValue(8.4) == Approx(8.4*8.4));
        }
    }
}
