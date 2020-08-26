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
 test_simionPotentialArray.cpp

 Testing of Simion Potential array reader

 ****************************/

#include "Core_debug.hpp"
#include "PSim_simionPotentialArray.hpp"
#include "test_util.hpp"
#include "catch.hpp"

const double scale_mm_to_m = 1000;

TEST_CASE("SIMION PA basic tests","[SimionPotentialArray]"){

    SECTION("Test with unscaled mirrored planar 2d PA") {
        std::string paFilename = "simion_test_planar_2d.pa";
        ParticleSimulation::SimionPotentialArray simPa(paFilename);

        REQUIRE(simPa.getBounds() == std::array<double,6>{-0.099, 0.099, -0.049, 0.049, 0.0, 0.0});

        REQUIRE(simPa.getNumberOfGridPoints() == std::array<ParticleSimulation::index_t, 3>{100,50,1});

        REQUIRE_THROWS_MATCHES(simPa.getInterpolatedPotential(1.0,2.0,0.0),
                ParticleSimulation::PotentialArrayException,
                Catch::Matchers::Contains("is not in the planar potential array"));

        REQUIRE_THROWS_MATCHES(simPa.getInterpolatedPotential(0.1,0.04,0.0),
                ParticleSimulation::PotentialArrayException,
                Catch::Matchers::Contains("is not in the planar potential array"));

        REQUIRE(simPa.getInterpolatedPotential(0.099,0.04,0.0) == Approx(19.79632));

        // Bounds checks are only performed with safety guards on:
        Core::safetyGuards = false;
        REQUIRE_NOTHROW(simPa.getInterpolatedPotential(0.1,0.04,0.0));
        Core::safetyGuards = true;
    }

    SECTION("Test with planar 3d PA") {
        std::string paFilename = "simion_test_planar_3d.pa";
        ParticleSimulation::SimionPotentialArray simPa(paFilename);

        REQUIRE(simPa.getBounds() == std::array<double,6>{0.0, 0.039, 0.0, 0.019, 0.0, 0.039});

        REQUIRE_THROWS_MATCHES(simPa.isElectrode(0.038, 0.018, 0.0390001),
                ParticleSimulation::PotentialArrayException,
                Catch::Matchers::Contains("is not in the planar potential array"));

        REQUIRE (simPa.isElectrode(0.02, 0.015, 0.03) == false);
        REQUIRE (simPa.isElectrode(0.038, 0.018, 0.038) == true);
        REQUIRE (simPa.isElectrode(0.038, 0.018, 0.0389999) == true);

        // uppper and lower boundaries:
        REQUIRE (simPa.isElectrode(0.038, 0.018, 0.039) == false);
        REQUIRE (simPa.isElectrode(0.0, 0.0001, 0.0) == true);

        REQUIRE_THROWS_MATCHES(simPa.getInterpolatedPotential(1.0,2.0,1.0),
                ParticleSimulation::PotentialArrayException,
                Catch::Matchers::Contains("not in the planar potential array"));

        REQUIRE(simPa.getInterpolatedPotential(0.02, 0.015, 0.03) == Approx(-85.37619));
        REQUIRE(simPa.getInterpolatedPotential(0.038, 0.018, 0.038) == Approx(10.0));

        REQUIRE(
                vectorApproxCompare(
                        simPa.getField(0.02, 0.015, 0.03),
                        Core::Vector(-113.889,-3781.03,-8.519e-5))
                        ==  vectorsApproxEqual);
    }

    SECTION("Test with planar 3d PA with spatial scaling and translation") {
        std::string paFilename = "simion_test_planar_3d.pa";
        ParticleSimulation::SimionPotentialArray simPaScaled(paFilename,{10.0,10.0,0.0}, 1.0, 100.0);
        ParticleSimulation::SimionPotentialArray simPa(paFilename);

        REQUIRE(simPa.getBounds() == std::array<double,6>{0.0, 0.039, 0.0, 0.019, 0.0, 0.039});
        REQUIRE(simPaScaled.getBounds() == std::array<double,6>{10.0, 49.0, 10.0, 29.0, 0.0, 39.0});


        REQUIRE(simPaScaled.getInterpolatedPotential(30.0, 25.0, 30.0) == Approx(-8537.619));

        REQUIRE(
            vectorApproxCompare(
                    simPaScaled.getField(30.0, 25.0, 30.0),
                    simPa.getField(0.02, 0.015, 0.03)*0.1)
                    ==  vectorsApproxEqual);
    }

    SECTION("Test with planar 3d PA with potential down scaling") {
        std::string paFilename = "simion_test_planar_3d.pa";

        double potScaleFactor = 0.001;
        ParticleSimulation::SimionPotentialArray simPaScaled(paFilename, 0.001);
        ParticleSimulation::SimionPotentialArray simPa(paFilename);

        //check if scaled potential is correct on electrode
        REQUIRE(
            simPaScaled.getInterpolatedPotential(0.038, 0.018, 0.038) / potScaleFactor ==
            Approx(simPa.getInterpolatedPotential(0.038, 0.018, 0.038)));

        //... and on a arbitraty position which is not an electrode:
        REQUIRE(
                simPaScaled.getInterpolatedPotential(0.02, 0.015, 0.03) / potScaleFactor ==
                Approx(simPa.getInterpolatedPotential(0.02, 0.015, 0.03)));


        REQUIRE(
                vectorApproxCompare(
                        simPa.getField(0.0168, 0.012, 0.0)*potScaleFactor,
                        simPaScaled.getField(0.0168, 0.012, 0.0))
                        ==  vectorsApproxEqual);
    }

    SECTION("Test with cylindrical PA") {
        std::string paFilename = "simion_test_cylindrical.pa";
        ParticleSimulation::SimionPotentialArray simPa(paFilename);

        std::array<double,6> correctBounds{0.0, 0.149, -0.059, 0.059, -0.059, 0.059};
        std::array<double,6> bounds = simPa.getBounds();
        for (int i=0; i<6; ++i){
            REQUIRE( Approx(bounds[i]) == correctBounds[i]);
        }


        REQUIRE( simPa.isInside(0.07, 0.02, 0.03) );
        REQUIRE( ! simPa.isInside(0.07, 0.7, 0.03) );

        REQUIRE_THROWS_MATCHES(simPa.getInterpolatedPotential(1.0,2.0,0.0),
                ParticleSimulation::PotentialArrayException,
                Catch::Matchers::Contains("is not in the cylindrical potential array"));

        REQUIRE(simPa.getInterpolatedPotential(0.07, 0.02, 0.03) == Approx(72.57148));
        //check for circular mirroring / circular symmetry:
        REQUIRE(simPa.getInterpolatedPotential(0.07, 0.01, 0.02) == simPa.getInterpolatedPotential(0.07, -0.01, -0.02));

        REQUIRE_THROWS_MATCHES(simPa.getField(1.0,2.0,0.0),
                ParticleSimulation::PotentialArrayException,
                Catch::Matchers::Contains("is not in the cylindrical potential array"));

        REQUIRE_THROWS_MATCHES(simPa.getField(-0.0001,0.0001,0.0),
                ParticleSimulation::PotentialArrayException,
                Catch::Matchers::Contains("is not in the cylindrical potential array"));

        REQUIRE(simPa.getInterpolatedPotential(0.14899, 0.002, 0.002) ==
            simPa.getInterpolatedPotential(0.14899, 0.002, 0.002));

        REQUIRE(
            vectorApproxCompare(
                    simPa.getField(0.14899, 0.002, 0.002),
                    Core::Vector(0.601785,1.59432,1.59432)*scale_mm_to_m)
                    ==  vectorsApproxEqual);

        REQUIRE(
            vectorApproxCompare(
                    simPa.getField(0.109, 0.01, 0.0),
                    Core::Vector(5.3374,-3.4043,0)*scale_mm_to_m)
                    ==  vectorsApproxEqual);
    }
}
