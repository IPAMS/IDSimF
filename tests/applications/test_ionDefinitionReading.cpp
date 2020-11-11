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
 test_ionDefinitionReading.cpp

 Testing of application utilites for ion definition reading

 ****************************/

#include "ionDefinitionReading.hpp"
#include "parameterParsing.hpp"
#include "catch.hpp"

TEST_CASE( "Test basic ion definition reading", "[ApplicationUtils]"){

    Json::Value conf_box = readConfigurationJson("ionBox.json");
    Json::Value conf_ionCloud = readConfigurationJson("ionCloudFile.json");

    REQUIRE_FALSE(AppUtils::isIonCloudDefinitionPresent(conf_box));
    REQUIRE(AppUtils::isIonCloudDefinitionPresent(conf_ionCloud));
}
