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
#include <iostream>
#include "catch.hpp"


TEST_CASE("Test sampled waveform / SWIFT waveform reader", "[ParticleSimulation][SampledWaveform][file readers]") {

    SECTION("SWIFT waveform data should be readable and correct") {
        ParticleSimulation::SampledWaveform sw("swift_test_transient.csv");

        SECTION("State should be good and correct") {
            REQUIRE(sw.good());
            REQUIRE(sw[15500] == 0.76383);
            REQUIRE(sw.getValue(16500) == -0.8761);
        }

        SECTION("Read from an non existing index should throw an exception"){
            REQUIRE(sw.good() == true);
            REQUIRE_THROWS(sw[100000000]);
        }
    }

    SECTION("Synthetically numpy generated SWIFT waveform should be readable and correct") {
        ParticleSimulation::SampledWaveform sw("swift_test_sin.csv");
        REQUIRE(sw.good());
        REQUIRE(sw[24] == 1.253332335643040918e-01 );
        REQUIRE(sw.getValue(24) == 1.253332335643040918e-01 );

        REQUIRE(sw.size() == 200);
        REQUIRE(sw.getValue(4) == sw.getValueLooped(204) );
        REQUIRE(sw.getValue(4) == sw.getValueLooped(404) );
        REQUIRE(sw.getValue(0) == sw.getValueLooped(200) );
    }

    SECTION( "Non existing input file should lead to good()==false") {
        ParticleSimulation::SampledWaveform sw("not_a_file");
        REQUIRE(! sw.good());
    }
}
