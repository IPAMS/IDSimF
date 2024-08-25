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
 test_simulation_configuration.cpp

 Tests of simulation configuration class / object

 ****************************/

#include "appUtils_simulationConfiguration.hpp"
#include "catch.hpp"
#include "test_util.hpp"
#include "Core_utils.hpp"
#include <string>

TEST_CASE( "Test simulation configuration ", "[ApplicationUtils]") {

    SECTION("Sim Conf reading: Types should be readable") {

        AppUtils::SimulationConfiguration simConf("simulationConfiguration_typesTest.json");

        CHECK(simConf.intParameter("integer_parameter") == 1500);
        CHECK(simConf.intVectorParameter("integer_vector") == std::vector<int>({8000, 9000}));

        CHECK(Core::isDoubleEqual(simConf.doubleParameter("double_parameter"), 500.5));
        CHECK(simConf.doubleVectorParameter("double_vector") == std::vector<double>({35.5, 100.1, 45.5}));

        CHECK(simConf.stringParameter("string_parameter") == "cylinder");
        CHECK(simConf.stringVectorParameter("string_vector") == std::vector<std::string>{"string_1", "string_2", "string_3"});

        CHECK(simConf.isParameter("integer_parameter"));
        CHECK( !simConf.isParameter("not_a_parameter"));

        CHECK(simConf.isVectorParameter("integer_vector"));
        CHECK( !simConf.isVectorParameter("integer_parameter"));

        CHECK( !simConf.isVectorParameter("not_a_parameter"));
    }

    SECTION("Potential Array reading convienence function should be functional") {

        AppUtils::SimulationConfiguration simConf("simulationConfiguration_fileReadingTest.json");

        std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> potentialArrays =
            simConf.readPotentialArrays("potential_arrays", 1.0);

        CHECK(potentialArrays[0]->getBounds()[0] == Approx(0.0));
        CHECK(potentialArrays[0]->getBounds()[1] == Approx(99.0));

        CHECK(potentialArrays[1]->getBounds()[3] == Approx(34.0));
        CHECK(potentialArrays[2]->getBounds()[3] == Approx(34.0));

        CHECK(potentialArrays[0]->getNumberOfGridPoints()[0] == 100);
        CHECK(potentialArrays[1]->getNumberOfGridPoints()[0] == 100);
        CHECK(potentialArrays[2]->getNumberOfGridPoints()[0] == 100);

        CHECK(potentialArrays[0]->getInterpolatedPotential(80.1, 10.1, 0.0) == Approx(2.1801610625));
        CHECK(potentialArrays[1]->getInterpolatedPotential(80.1, 10.1, 0.0) == Approx(4950.0400784919));
        CHECK(potentialArrays[2]->getInterpolatedPotential(80.1, 10.1, 0.0) == Approx(103.8599748651));
    }
}