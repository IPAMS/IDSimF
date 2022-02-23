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

#include "appUtils_ionDefinitionReading.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "catch.hpp"
#include "test_util.hpp"

bool testParticleBox(std::vector<std::unique_ptr<Core::Particle>>& particles){

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

bool testParticleCylinder(
        std::vector<std::unique_ptr<Core::Particle>>& particles,
        double expectedLength = 0.006,
        double expectedRadius = 0.004
        ){

    Core::Vector expectedShift(0.002, 0.005, -0.002);

    bool incorrectIonFound = false;

    for(const auto& part: particles){
        Core::Vector pos = part->getLocation();
        pos = pos - expectedShift;
        double r = std::sqrt(pos.x()*pos.x() + pos.y()*pos.y());

        if (pos.z() > expectedLength || pos.z() < 0.0 || r > expectedRadius){
            incorrectIonFound = true;
        }
    }

    return (incorrectIonFound);
}

TEST_CASE( "Test ion definition reading", "[ApplicationUtils]") {

    AppUtils::SimulationConfiguration simConf_box("ionBox.json");
    AppUtils::SimulationConfiguration simConf_cylinder_minimal("ionCylinder_minimal.json");
    AppUtils::SimulationConfiguration simConf_cylinder_full("ionCylinder_full.json");
    AppUtils::SimulationConfiguration simConf_ionCloud("ionCloudFile.json");

    //std::string configurationPath = "./";

    std::vector<std::unique_ptr<Core::Particle>> particles;
    std::vector<Core::Particle*> particlePtrs;

    SECTION("Ion definition reading: Ion definition reading with invalid ion start geometry should throw") {
        AppUtils::SimulationConfiguration simConf_invalid("ionDefinition_invalid.json");
        REQUIRE_THROWS_AS(
                AppUtils::readIonDefinition(particles, particlePtrs, simConf_invalid),
                std::invalid_argument);
    }

    SECTION("Ion definition reading: Ion definition reading with invalid ion number should throw") {
        AppUtils::SimulationConfiguration simConf_invalid("ionDefinition_invalid_2.json");
        REQUIRE_THROWS_WITH(
                AppUtils::readIonDefinition(particles, particlePtrs, simConf_invalid),
                "LargestInt out of UInt range");
    }

    SECTION("Test ion definitions with ion cloud files") {

        SECTION("Ion definition reading: If ion cloud file is present should be recognized") {
            REQUIRE_FALSE(AppUtils::isIonCloudDefinitionPresent(simConf_box));
            CHECK(AppUtils::isIonCloudDefinitionPresent(simConf_ionCloud));
        }

        SECTION("Ion definition reading: Ion cloud file should be readable") {
            AppUtils::readIonDefinitionFromIonCloudFile(particles, particlePtrs, simConf_ionCloud);

            CHECK(particles[1]->getCharge()==Approx(-1.0*Core::ELEMENTARY_CHARGE));
            CHECK(vectorApproxCompare(particles[1]->getLocation(), Core::Vector(1.0, 2.0, 1.0))
                    =="Vectors approximately equal");
        }

        SECTION("Ion definition reading: Full ion definition reading with ion cloud file should work") {
            AppUtils::readIonDefinition(particles, particlePtrs, simConf_ionCloud);
            CHECK(vectorApproxCompare(particles[1]->getLocation(), Core::Vector(1.0, 2.0, 1.0))
                    =="Vectors approximately equal");
        }
    }

    SECTION("Test ion definitions with box start zone") {

        SECTION("Ion definition reading: Random ion box definition should be readable") {
            AppUtils::readRandomIonDefinition(particles, particlePtrs, simConf_box);
            bool particleOutOfBoxFound = testParticleBox(particles);
            CHECK(!particleOutOfBoxFound);
        }

        SECTION("Ion definition reading: Full ion definition reading with random box should work") {
            AppUtils::readIonDefinition(particles, particlePtrs, simConf_box);
            bool particleOutOfBoxFound = testParticleBox(particles);
            CHECK(!particleOutOfBoxFound);


        }
    }

    SECTION("Test ion definitions with minimal cylinder start zone") {

        SECTION("Ion definition reading: Minimal random ion cylinder definition should be readable") {
            AppUtils::readRandomIonDefinition(particles, particlePtrs, simConf_cylinder_minimal);
            bool particleOutOfCylinderFound = testParticleCylinder(particles);
            CHECK(!particleOutOfCylinderFound);
        }

        SECTION("Ion definition reading: Full ion definition reading with minimal cylinder definition should work") {
            AppUtils::readIonDefinition(particles, particlePtrs, simConf_cylinder_minimal);
            bool particleOutOfCylinderFound = testParticleCylinder(particles);
            CHECK(!particleOutOfCylinderFound);
        }
    }

    SECTION("Test ion definitions with full cylinder start zone") {

        SECTION("Ion definition reading: Full ion definition reading with full cylinder definition should work") {
            AppUtils::readIonDefinition(particles, particlePtrs, simConf_cylinder_full);
            bool particleOutOfCylinderFound = testParticleCylinder(particles, 0.08, 0.01);
            REQUIRE(!particleOutOfCylinderFound);

            auto& p_35 = particles[1];
            auto& p_100 = particles[8001];
            REQUIRE(p_35->getMass() == Approx(35*Core::AMU_TO_KG));
            REQUIRE(p_100->getMass() == Approx(100*Core::AMU_TO_KG));

            double kE_35 = 1/2.0 * p_35->getMass() * p_35->getVelocity().magnitudeSquared() * Core::JOULE_TO_EV;
            double kE_100 = 1/2.0 * p_100->getMass() * p_100->getVelocity().magnitudeSquared() * Core::JOULE_TO_EV;
            REQUIRE(kE_35 == Approx(100.0));
            REQUIRE(kE_100 == Approx(100.0));

            Core::Vector v_100 = p_100->getVelocity();
            REQUIRE(isExactDoubleEqual(v_100.z(), 0.0));
            REQUIRE(2.0 * v_100.x() == Approx(v_100.y()));
        }
    }
}
