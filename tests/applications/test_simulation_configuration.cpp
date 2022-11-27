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
#include <string>

TEST_CASE( "Test simulation configuration ", "[ApplicationUtils]") {

    AppUtils::SimulationConfiguration simConf("simulationConfiguration_typesTest.json");

    SECTION("Sim Conf reading: Types should be readable") {
        CHECK(simConf.intParameter("integer_parameter") == 1500);
        CHECK(simConf.intVectorParameter("integer_vector") == std::vector<int>({8000, 9000}));

        CHECK(simConf.doubleParameter("double_parameter") == 500.5);
        CHECK(simConf.doubleVectorParameter("double_vector") == std::vector<double>({35.5, 100.1, 45.5}));

        CHECK(simConf.stringParameter("string_parameter") == "cylinder");
        CHECK(simConf.stringVectorParameter("string_vector") == std::vector<std::string>{"string_1", "string_2", "string_3"});

        CHECK(simConf.isParameter("integer_parameter"));
        CHECK( !simConf.isParameter("not_a_parameter"));

        CHECK(simConf.isVectorParameter("integer_vector"));
        CHECK( !simConf.isVectorParameter("integer_parameter"));

        CHECK( !simConf.isVectorParameter("not_a_parameter"));
    }

}