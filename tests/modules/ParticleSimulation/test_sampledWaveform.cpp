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
 test_sampledWaveform.cpp

 Testing of sampled waveform / SWIFT waveform reader

 ****************************/

#include "PSim_sampledWaveform.hpp"
#include "catch.hpp"
#include <iostream>


TEST_CASE("Test sampled waveform / SWIFT waveform reader", "[ParticleSimulation][SampledWaveform][file readers]") {

    SECTION("SWIFT waveform data should be readable and correct") {
        ParticleSimulation::SampledWaveform sw("swift_test_transient.csv");

        SECTION("State should be good and correct") {
            CHECK(sw.good());
            CHECK(sw[15500] == Approx(0.76383));
            CHECK(sw.getValue(16500) == Approx(-0.8761));
        }

        SECTION("Read from an non existing index should throw an exception"){
            CHECK(sw.good() == true);
            CHECK_THROWS(sw[100000000]);
        }
    }

    SECTION("Synthetically numpy generated SWIFT waveform should be readable and correct") {
        ParticleSimulation::SampledWaveform sw("swift_test_sin.csv");
        CHECK(sw.good());
        CHECK(sw[24] == Approx(1.253332335643040918e-01) );
        CHECK(sw.getValue(24) == Approx(1.253332335643040918e-01) );

        CHECK(sw.size() == 200);
        CHECK(sw.getValue(4) == Approx(sw.getValueLooped(204)) );
        CHECK(sw.getValue(4) == Approx(sw.getValueLooped(404)) );
        CHECK(sw.getValue(0) == Approx(sw.getValueLooped(200)) );
    }

    SECTION("Waveform with low sample number should be interpolatable") {
        ParticleSimulation::SampledWaveform sw("low_sample_waveform.csv");
        CHECK(sw.good());
        CHECK(sw.size() == 20);

        CHECK(sw.getInterpolatedValue(0.2) == Approx(1.0));
        CHECK(sw.getInterpolatedValue(0.205) == Approx(1.0));
        CHECK(sw.getInterpolatedValue(0.35) == Approx(1.0));
        CHECK(sw.getInterpolatedValue(0.375) == Approx(1.5));
        CHECK(sw.getInterpolatedValue(0.4) == Approx(2.0));
        CHECK(sw.getInterpolatedValue(0.45) == Approx(2.0));
        CHECK(sw.getInterpolatedValue(0.475) == Approx(3.0));
        CHECK(sw.getInterpolatedValue(0.50) == Approx(4.0));
        CHECK(sw.getInterpolatedValue(0.50) == Approx(4.0));
        CHECK(sw.getInterpolatedValue(0.975) == Approx(0.5));
        CHECK_THROWS(sw.getInterpolatedValue(1.1));
    }

    SECTION("Triangle Waveform should work") {
        ParticleSimulation::SampledWaveform sw("triangle_waveform.csv");
        CHECK(sw.good());
        CHECK(sw.size() == 2);
        CHECK(sw.getInterpolatedValue(0.0) == Approx(0.0));
        CHECK(sw.getInterpolatedValue(0.25) == Approx(0.5));
        CHECK(sw.getInterpolatedValue(0.5) == Approx(1.0));
        CHECK(sw.getInterpolatedValue(0.75) == Approx(0.5));
        CHECK(sw.getInterpolatedValue(0.99) == Approx(0.0));
    }

    SECTION( "Non existing input file should lead to good()==false") {
        ParticleSimulation::SampledWaveform sw("not_a_file");
        CHECK(! sw.good());
    }
}