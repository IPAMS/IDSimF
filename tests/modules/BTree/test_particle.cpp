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
 test_particle.cpp

 Testing of basic simulated particles

 ****************************/

#include "BTree_particle.hpp"
#include "BTree_node.hpp"
#include "Core_constants.hpp"
#include "RS_Substance.hpp"
#include "catch.hpp"

TEST_CASE( "Basic Particle semantics tests", "[Particle]") {

    BTree::Particle testIon = BTree::Particle(Core::Vector(1.0, 2.0, 3.0), 2.0);

    SECTION("Basic Particle creation tests") {
        //Test ion location:
        REQUIRE(testIon.getLocation()==Core::Vector(1.0, 2.0, 3.0));
        REQUIRE(testIon.getCharge()==2.0*Core::ELEMENTARY_CHARGE);
        REQUIRE(testIon.isActive());
    }

    SECTION("Test particle constructors") {

        SECTION("Test constructor with time of birth"){
            //Test ion location:
            BTree::Particle testIon2 = BTree::Particle(
                    Core::Vector(2.0, 3.0, 4.0),
                    Core::Vector(1.0, 1.0, 1.0),
                    1.0,
                    100.0,
                    1e-5);
            REQUIRE(testIon2.getLocation()==Core::Vector(2.0, 3.0, 4.0));
            REQUIRE(testIon2.getCharge()==1.0*Core::ELEMENTARY_CHARGE);
            REQUIRE(testIon2.getMass()==100.0*Core::AMU_TO_KG);
            REQUIRE(testIon2.getTimeOfBirth()==1e-5);
            REQUIRE(testIon2.isActive());
        }

        SECTION("Test constructor with time of birth and collision diameter"){
            //Test ion location:
            BTree::Particle testIon2 = BTree::Particle(
                    Core::Vector(2.0, 3.0, 4.0),
                    Core::Vector(1.0, 1.0, 1.0),
                    1.0,
                    100.0,
                    3.2e-10,
                    1e-5);
            REQUIRE(testIon2.getLocation()==Core::Vector(2.0, 3.0, 4.0));
            REQUIRE(testIon2.getCharge()==1.0*Core::ELEMENTARY_CHARGE);
            REQUIRE(testIon2.getMass()==100.0*Core::AMU_TO_KG);
            REQUIRE(testIon2.getDiameter()==3.2e-10);
            REQUIRE(testIon2.getTimeOfBirth()==1e-5);
            REQUIRE(testIon2.isActive());
        }

    }

    SECTION("Basic Particle setter tests") {

        testIon.setActive(false);
        REQUIRE(!testIon.isActive());

        testIon.setChargeElementary(1.0);
        REQUIRE(testIon.getCharge()==1.0*Core::ELEMENTARY_CHARGE);

        testIon.setLocation(Core::Vector(2.0, 2.0, 2.0));
        REQUIRE(testIon.getLocation()==Core::Vector(2.0, 2.0, 2.0));

        REQUIRE(testIon.getHostNode()==nullptr);

        BTree::Node testNode = BTree::Node(
                Core::Vector(1.0, 1.0, 1.0),
                Core::Vector(2.0, 2.0, 2.0),
                nullptr
        );

        testIon.setHostNode(&testNode);
        REQUIRE(testIon.getHostNode()==&testNode);
    }

    SECTION("Particle location modification possible, location is a reference") {
        testIon.getLocation().set(10.0, 20.0, 30.0);
        REQUIRE(testIon.getLocation()==Core::Vector(10.0, 20.0, 30.0));
    }

    SECTION("Auxiliary parameters are accessible as parameters") {

        std::string testKey = "a testkey";

        testIon.setAuxScalarParam(testKey, 100.0);
        REQUIRE(testIon.getAuxScalarParam(testKey)==100.0);

        testIon.setAuxScalarParam(testKey, 200.0);
        REQUIRE(testIon.getAuxScalarParam(testKey)==200.0);
    }
}
