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
 test_simulationConfiguration.cpp

 Testing of reaction simulation configuration in RS

 ****************************/

#include "RS_SimulationConfiguration.hpp"
#include "catch.hpp"

TEST_CASE("Test RS simulation configuration", "[RS][SimulationConfiguration]") {

    RS::SimulationConfiguration simConf = RS::SimulationConfiguration();

    std::unique_ptr<RS::Substance> Cluster_1 = std::unique_ptr<RS::Substance>(
            new RS::Substance("[H3O]+",RS::Substance::substanceType::discrete));

    std::unique_ptr<RS::Substance> Cluster_2 = std::unique_ptr<RS::Substance>(
            new RS::Substance("Cl_2",RS::Substance::substanceType::discrete));

    std::unique_ptr<RS::Substance> N2 = std::unique_ptr<RS::Substance>(
            new RS::Substance("N2",RS::Substance::substanceType::isotropic));


    SECTION("Substances can be added to SimulationConfiguration and retrieved by index") {
        simConf.addSubstance(Cluster_1);
        simConf.addSubstance(N2);

        CHECK(simConf.substance(0)->name() == "[H3O]+");
        CHECK(simConf.substance(0)->type() == RS::Substance::substanceType::discrete);
        CHECK(simConf.substance(1)->name() == "N2");
        CHECK(simConf.substance(1)->type() == RS::Substance::substanceType::isotropic);
    }

    SECTION("Substances can be added to SimulationConfiguration and retrieved by name") {
        simConf.addSubstance(Cluster_1);
        simConf.addSubstance(N2);

        CHECK(simConf.substanceByName("[H3O]+")->name() == "[H3O]+");
        CHECK(simConf.substanceByName("[H3O]+")->type() == RS::Substance::substanceType::discrete);
        CHECK(simConf.substanceByName("N2")->name() == "N2");
        CHECK(simConf.substanceByName("N2")->type() == RS::Substance::substanceType::isotropic);
    }

    SECTION("Discrete Substances can be retrieved by name") {
        simConf.addSubstance(Cluster_1);
        simConf.addSubstance(N2);
        simConf.addSubstance(Cluster_2);

        CHECK(simConf.getAllDiscreteSubstances().size() == 2);
        CHECK(simConf.getAllDiscreteSubstances()[0]->name() == "[H3O]+");
        CHECK(simConf.getAllDiscreteSubstances()[1]->name() == "Cl_2");
    }

    SECTION("Substances with the same name cannot be added") {
        std::unique_ptr<RS::Substance> Cluster1_duplicate = std::unique_ptr<RS::Substance>(
                new RS::Substance("[H3O]+",RS::Substance::substanceType::isotropic));

        CHECK(simConf.addSubstance(Cluster_1) == true);
        CHECK(simConf.addSubstance(Cluster1_duplicate) == false);

        CHECK(simConf.substanceByName("[H3O]+")->name() == "[H3O]+");
        CHECK(simConf.substanceByName("[H3O]+")->type() == RS::Substance::substanceType::discrete);
        REQUIRE(simConf.getAllSubstances().size() == 1);
    }

    SECTION( "Retrieval of non existing substances is not possible") {
        REQUIRE(simConf.addSubstance(Cluster_1) == true);
        REQUIRE_THROWS_WITH(simConf.substanceByName("notexisting"), "Substance notexisting is not existing");
        REQUIRE_THROWS_AS(simConf.substance(10), std::out_of_range);
    }
}