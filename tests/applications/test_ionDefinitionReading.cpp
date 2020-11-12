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
#include "test_util.hpp"
#include "catch.hpp"

bool testParticleBox(std::vector<std::unique_ptr<BTree::Particle>>& particles){

    bool particleOutOfBoxFound = false;

    for(const auto& part: particles){
        Core::Vector pLocation = part->getLocation();
        if (pLocation.x() < 0.001 || pLocation.x() > 0.003 ||
                pLocation.y() < 0.001 || pLocation.y() > 0.003 ||
                pLocation.z() < 0.001 || pLocation.z() > 0.003 )
        {
            particleOutOfBoxFound = true;
        }
    }

    return (particleOutOfBoxFound);
}

TEST_CASE( "Test basic ion definition reading", "[ApplicationUtils]"){

    Json::Value conf_box = readConfigurationJson("ionBox.json");
    Json::Value conf_ionCloud = readConfigurationJson("ionCloudFile.json");

    std::vector<std::unique_ptr<BTree::Particle>>particles;
    std::vector<BTree::Particle*>particlePtrs;

    SECTION("Ion definition reading: If ion cloud file is present should be recognized"){
        REQUIRE_FALSE(AppUtils::isIonCloudDefinitionPresent(conf_box));
        REQUIRE(AppUtils::isIonCloudDefinitionPresent(conf_ionCloud));
    }

    SECTION("Ion definition reading: Ion cloud file should be readable"){
        AppUtils::readIonDefinitionFromIonCloudFile(particles, particlePtrs, conf_ionCloud);

        REQUIRE(particles[1]->getCharge() == Approx(1.0 * Core::ELEMENTARY_CHARGE));
        REQUIRE(vectorApproxCompare(particles[1]->getLocation(), Core::Vector(-0.001, 0.0, 0.0)) == "Vectors approximately equal");
    }

    SECTION("Ion definition reading: Random ion box definition should be readable"){
        AppUtils::readRandomIonDefinition(particles, particlePtrs, conf_box);
        bool particleOutOfBoxFound = testParticleBox(particles);
        REQUIRE( !particleOutOfBoxFound );
    }

    SECTION("Ion definition reading: Full ion definition reading with ion cloud file should work"){
        AppUtils::readIonDefinition(particles, particlePtrs, conf_ionCloud);
        REQUIRE(vectorApproxCompare(particles[1]->getLocation(), Core::Vector(-0.001, 0.0, 0.0)) == "Vectors approximately equal");
    }

    SECTION("Ion definition reading: Full ion definition reading with random box should work"){
        AppUtils::readIonDefinition(particles, particlePtrs, conf_box);
        bool particleOutOfBoxFound = testParticleBox(particles);
        REQUIRE( !particleOutOfBoxFound );
    }
}
