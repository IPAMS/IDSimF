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
 test_substance.cpp

 Testing of chemical substance in RS

 ****************************/

#include "RS_Substance.hpp"
#include "catch.hpp"
#include "test_util.hpp"

TEST_CASE("RS Substance: Test basic instantiation", "[RS][Substance]") {

    RS::Substance sub("Testsubstance", RS::Substance::substanceType::discrete);

    SECTION("RS Substance: Test basic instantiation") {
        REQUIRE(sub.name() == "Testsubstance");
        REQUIRE(sub.type() == RS::Substance::substanceType::discrete);

        sub= RS::Substance("Testsubstance 2", "isotropic");
        REQUIRE(sub.name() == "Testsubstance 2");
        REQUIRE(sub.type() == RS::Substance::substanceType::isotropic);

        sub= RS::Substance("Testsubstance 3", "field");
        REQUIRE(sub.name() == "Testsubstance 3");
        REQUIRE(sub.type() == RS::Substance::substanceType::field);

        REQUIRE_THROWS(sub= RS::Substance("Testsubstance 2", "not_existing"));
    }

    SECTION("RS Substance: Test basic getter and setters") {
        sub.charge(10.0);
        REQUIRE(isExactDoubleEqual(sub.charge(), 10.0));

        sub.mass(100.0);
        REQUIRE(isExactDoubleEqual(sub.mass(), 100.0));

        sub.staticConcentration(1.0e10);
        REQUIRE(isExactDoubleEqual(sub.staticConcentration(),  1.0e10));

        sub.lowFieldMobility(5.0);
        REQUIRE(isExactDoubleEqual(sub.lowFieldMobility(), 5.0));
    }
}
