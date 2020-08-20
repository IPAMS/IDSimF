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
 test_ionCloudReader.cpp

 Testing of ion cloud reader implementations

 ****************************/

#include "PSim_ionCloudReader.hpp"
#include "Core_vector.hpp"
#include "BTree_particle.hpp"
#include <vector>
#include "catch.hpp"

TEST_CASE("Test ion cloud reader", "[ParticleSimulation][IonCloudReader][file readers]") {

    SECTION( "Ion cloud reader: existing file should open without exception") {
        ParticleSimulation::IonCloudReader reader = ParticleSimulation::IonCloudReader();
        REQUIRE_NOTHROW(reader.readIonCloud("test_ion_cloud_01.csv"));
    }

    SECTION( "Ion cloud reader: non existing file should throw an exception") {
        ParticleSimulation::IonCloudReader reader = ParticleSimulation::IonCloudReader();
        REQUIRE_THROWS(reader.readIonCloud("i_do_not_exist.csv"));
    }

    SECTION( "Ion cloud reader: ions defined in file are parsed correctly") {
        ParticleSimulation::IonCloudReader reader = ParticleSimulation::IonCloudReader();
        std::vector<std::unique_ptr<BTree::Particle>> iCl = reader.readIonCloud("test_ion_cloud_01.csv");
        std::vector<BTree::Particle> iClRef = std::vector<BTree::Particle>();
        iClRef.push_back(BTree::Particle(Core::Vector(1.0,1.0,1.0),Core::Vector(1.00,1.00,1.00),1.0,100.0));
        iClRef.push_back(BTree::Particle(Core::Vector(1.0,2.0,1.0),Core::Vector(10.0,10.0,10.0),-1.0,200.0));
        iClRef.push_back(BTree::Particle(Core::Vector(-10,-20,-10.0),Core::Vector(-10.00,10.0,-10.0),2.0,300.0));
        iClRef.push_back(BTree::Particle(Core::Vector(1.0,2.0,1.0),Core::Vector(10.0,10.0,10.0),-10.5,200.0));

        REQUIRE(iCl.size() == 4);

        for (int i=0; i<iCl.size(); i++){
            REQUIRE( (*iCl[i]).getLocation() == iClRef[i].getLocation() );
            REQUIRE( (*iCl[i]).getVelocity() == iClRef[i].getVelocity() );
            REQUIRE( (*iCl[i]).getCharge() == iClRef[i].getCharge() );
            REQUIRE( (*iCl[i]).getMass() == iClRef[i].getMass() );
        }
    }
}