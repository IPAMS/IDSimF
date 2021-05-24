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
 test_configFileParser.cpp

 Testing of RS config file parser

 ****************************/

#include "catch.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_Substance.hpp"
#include "RS_AbstractReaction.hpp"
#include "RS_VantHoffReaction.hpp"
#include "RS_SimpleCollisionStepReaction.hpp"
#include <vector>
#include <utility>

using sMap = std::map<RS::Substance*,int>;
using sPair= sMap::value_type;


TEST_CASE("Test basic RS config file parser file reading", "[RS][ConfigFileParser][file readers]") {

    RS::ConfigFileParser parser;

    SECTION("RS Config file parser: existing file should parse without exception") {
        REQUIRE_NOTHROW(parser.parseFile("RS_waterCluster_test.conf"));
    }

    SECTION("RS Config file parser: non existing file should throw an exception") {
        REQUIRE_THROWS(parser.parseFile("i_do_not_exist.conf"));
    }
}

TEST_CASE("Test RS config file parser methods", "[RS][ConfigFileParser]") {

    RS::ConfigFileParser parser;

    SECTION( "RS Config file parser: test string splitting") {

        std::string testString = "ABC[token1]123[token2]DEF";
        std::string pattern = "\\[\\w+\\]";
        std::pair<std::vector<std::string>,std::vector<std::string>> splitResult = parser.splitString(testString,pattern);
        std::vector<std::string> result1 = {"ABC","123","DEF"};
        REQUIRE(splitResult.first == result1);


        testString = "[token1]XYZ[token2]ABC";
        splitResult = parser.splitString(testString,pattern);
        std::vector<std::string> result2 = {"XYZ","ABC"};
        REQUIRE(splitResult.first == result2);


        testString = "XYZ[token2]ABC[token3]DEF[token4]";
        splitResult = parser.splitString(testString,pattern);
        std::vector<std::string> result3 = {"XYZ","ABC","DEF"};
        REQUIRE(splitResult.first == result3);
    }
}


TEST_CASE("Test parsing of chemical systems with RS config file parser", "[RS][ConfigFileParser][file readers]") {

    RS::ConfigFileParser parser;
    RS::Substance ed_1("educt_1",RS::Substance::substanceType::discrete);
    RS::Substance pro_1("product_1",RS::Substance::substanceType::discrete);
    RS::Substance N2("N2",RS::Substance::substanceType::isotropic);

    sMap educts;
    sMap products;

    educts.insert(sPair(&ed_1,0));
    educts.insert(sPair(&N2,1));
    products.insert(sPair(&pro_1,0));

    RS::ReactionConditions reactionConditions;
    RS::ReactiveParticle dummyParticle(&ed_1);

    SECTION("Config file with all reaction types should be parsed correctly") {
        N2.staticConcentration(3.58e16);
        std::unique_ptr<RS::SimulationConfiguration> simConf = parser.parseFile("RS_minimal_test.conf");

        RS::AbstractReaction* r0 = simConf->reaction(0);
        RS::AbstractReaction* r1 = simConf->reaction(1);
        RS::AbstractReaction* r3 = simConf->reaction(3);
        RS::AbstractReaction* r4 = simConf->reaction(4);
        RS::AbstractReaction* r5 = simConf->reaction(5);
        REQUIRE(r0->getTypeLabel() == "static");
        REQUIRE(r1->getTypeLabel() == "vanthoff");
        REQUIRE(r3->getTypeLabel() == "vanthoff_field");
        REQUIRE(r4->getTypeLabel() == "simple_step");
        REQUIRE(r5->getTypeLabel() == "static_thermalizing");

        //Test if the rate constants are calculated correctly

        RS::VantHoffReaction reac_compare_1(
                educts,products,
                132.0*1000.0,6.8408e17,2e-9,
                "a test reaction");
        reac_compare_1.updateStaticReactionConcentration();

        reactionConditions.temperature = 598;
        double k_test = reac_compare_1.attemptReaction(reactionConditions, &dummyParticle, 1.0).reactionProbability;
        double k_1 = r1->attemptReaction(reactionConditions, &dummyParticle, 1.0).reactionProbability;

        REQUIRE(Approx(k_1).epsilon(1e-1) == k_test);

        RS::FieldDependentVantHoffReaction reac_compare_3(
                educts,products,
                82000,
                3.3806e+09,
                2e-9,
                2.35e-4,
                28,
                "a field dependent vant hoff test reaction");

        reactionConditions.electricField= 10000;


        k_test = reac_compare_3.attemptReaction(reactionConditions, &dummyParticle, 1.0).reactionProbability;
        double k_3 = r3->attemptReaction(reactionConditions, &dummyParticle, 1.0).reactionProbability;

        REQUIRE(Approx(k_3) == k_test);

        RS::SimpleCollisionStepReaction reac_compare_4(
                educts,products,
                10,
                "a simple step test reaction");

        RS::CollisionConditions collisionConditions = {.totalCollisionEnergy = 1.0e-5};
        k_test = reac_compare_4.attemptReaction(collisionConditions, &dummyParticle).reactionProbability;
        double k_4 = r4->attemptReaction(collisionConditions, &dummyParticle).reactionProbability;

        REQUIRE(Approx(k_4) == k_test);

        dummyParticle.setMassAMU(200);
        RS::StaticThermalizingReaction reac_compare_5(
                educts,products,
                2.5e-28,
                "simple thermalizing reaction");
        k_test = reac_compare_5.attemptReaction(reactionConditions, &dummyParticle, 1e-5).reactionProbability;
        double k_5 = r5->attemptReaction(reactionConditions, &dummyParticle, 1e-5).reactionProbability;

        REQUIRE(Approx(k_5) == k_test);
    }

    SECTION ("Existing water cluster RS config file should parse correctly") {
        std::unique_ptr<RS::SimulationConfiguration> simConf = parser.parseFile("RS_waterCluster_test.conf");
        REQUIRE(simConf->substance(6)->type() == RS::Substance::substanceType::isotropic);
        REQUIRE(simConf->substance(6)->name() == "N2");
        REQUIRE(Approx(simConf->substance(6)->staticConcentration()) == 3.58e16);

        REQUIRE(simConf->substance(2)->name() == "Cl_3");
        REQUIRE(simConf->substance(2)->type() == RS::Substance::substanceType::discrete);
        REQUIRE(simConf->substance(2)->mass() == 55);
        REQUIRE(simConf->substance(2)->charge() == 1.0);
        REQUIRE(simConf->substance(2)->mobility() == Approx(0.000235).epsilon(1e-6));
        REQUIRE(simConf->substance(2)->collisionDiameter() == Approx(7.0e-10).epsilon(1e-5));

        REQUIRE(simConf->substance(3)->name() == "Cl_4");
        REQUIRE(simConf->substance(3)->type() == RS::Substance::substanceType::discrete);
        REQUIRE(simConf->substance(3)->mass() == 73);
        REQUIRE(simConf->substance(3)->charge() == 1.0);
        REQUIRE(simConf->substance(3)->mobility() == Approx(2.09e-4).epsilon(1e-6));
        REQUIRE(simConf->substance(3)->collisionDiameter() == Approx(9.0e-10).epsilon(1e-5));

        std::vector<RS::Substance*> discreteSubstances = simConf->getAllDiscreteSubstances();
        REQUIRE(discreteSubstances[2]->name() == "Cl_3");
        REQUIRE(discreteSubstances[3]->name() == "Cl_4");

        RS::AbstractReaction* r0 = simConf->reaction(0);
        REQUIRE(r0->getTypeLabel() == "static");
        REQUIRE(r0->getLabel() == "cl1_forward");
        REQUIRE(r0->discreteEducts()->at(simConf->substanceByName("Cl_1"))== 1);
        REQUIRE(r0->educts()->at(simConf->substanceByName("H2O"))== 1);

        REQUIRE(r0->discreteProducts()->at(simConf->substanceByName("Cl_2"))== 1);
        REQUIRE(r0->products()->at(simConf->substanceByName("N2"))== 1);
    }

    SECTION( "Existing field dependent water cluster RS config file should parse correctly") {
        std::unique_ptr<RS::SimulationConfiguration> simConf = parser.parseFile("RS_waterCluster_test_fieldDependent.conf");

        REQUIRE(simConf->substance(4)->name() == "Cl_5");
        REQUIRE(simConf->substance(4)->type() == RS::Substance::substanceType::discrete);
        REQUIRE(simConf->substance(4)->mass() == 91);
        REQUIRE(simConf->substance(4)->charge() == 1.0);
        REQUIRE(simConf->substance(4)->mobility() == Approx(1.90e-4).epsilon(1e-6));
        REQUIRE(simConf->substance(4)->collisionDiameter() == Approx(10.0e-10).epsilon(1e-5));

        std::vector<RS::Substance*> discreteSubstances = simConf->getAllDiscreteSubstances();
        REQUIRE(discreteSubstances[2]->name() == "Cl_3");

        RS::AbstractReaction* r8 = simConf->reaction(8);
        REQUIRE(r8->getTypeLabel() == "vanthoff_field");
        REQUIRE(r8->getLabel() == "cl2_backward");
        REQUIRE(r8->discreteEducts()->at(simConf->substanceByName("Cl_2")) == 1);
        REQUIRE_THROWS(r8->educts()->at(simConf->substanceByName("H2O")));
        REQUIRE(r8->educts()->at(simConf->substanceByName("N2"))== 1);

        RS::FieldDependentVantHoffReaction reac_compare(
                educts, products,
                132000,
                6.8408e+17,
                2e-9,
                3.56629e-4,
                28,
                "a field dependent vant hoff test reaction");

        reactionConditions.electricField= 10000;
    }

    SECTION("Exception conditions in substance parsing should throw exceptions"){

        SECTION("Missing concentration value for isotropic substance should throw exception"){
            std::string teststring= "[SUBSTANCES]\n"
                                    "N2 isotropic 3.58e16\n"
                                    "H2O isotropic\n "
                                    "[REACTIONS]\n";

            REQUIRE_THROWS_WITH(parser.parseText(teststring), "Isotropic substance without concentration value found");
        }

        SECTION("Missing values for discrete substance should throw exception") {
            std::string teststring = "[SUBSTANCES]\n"
                         "N2 discrete 10\n"
                         "H2O discrete 20 10\n"
                         "[REACTIONS]\n";
            REQUIRE_THROWS_WITH(parser.parseText(teststring), "Discrete substance: Mass, charge, mobility "
                                                              "or collision diameter value missing");
        }
    }

    SECTION("Exception conditions in reaction parsing should throw exceptions"){

        SECTION("Reaction with an missing educt or product should throw an exception"){
            std::string teststring= "[SUBSTANCES]\n"
                                    "N2 isotropic 3.58e16\n"
                                    "H2O isotropic 3.58e14\n"
                                    "Cl_1 discrete 19 1 3.57e-4 3.0e-10\n"
                                    "[REACTIONS]\n"
                                    "Cl_3 + H2O + N2 => Cl_1 + N2 | static ; 6.98E-29 #cl2_forward";

            REQUIRE_THROWS_WITH(parser.parseText(teststring), "Substance Cl_3 is not existing");

            teststring= "[SUBSTANCES]\n"
                        "N2 isotropic 3.58e16\n"
                        "H2O isotropic 3.58e14\n"
                        "Cl_1 discrete 19 1 3.57e-4 3.0e-10\n"
                        "[REACTIONS]\n"
                        "Cl_1 + H2O + N2 => Cl_2 + N2 | static ; 6.98E-29 #cl1_forward";

            REQUIRE_THROWS_WITH(parser.parseText(teststring), "Substance Cl_2 is not existing");
        }
    }
}

