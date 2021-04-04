/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2021 - Physical and Theoretical Chemistry /
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
 test_jsonFileWriter.cpp

 Testing of JSON trajectory file writer

 ****************************/


#include "catch.hpp"
#include "PSim_particleStartSplatTracker.hpp"
#include "BTree_particle.hpp"
#include "Core_vector.hpp"
#include <iostream>


TEST_CASE("TestParticleStartSplatTracker", "[ParticleSimulation][ParticleStartSplatTracker][particle trackers]") {
    ParticleSimulation::ParticleStartSplatTracker tracker;

    BTree::Particle particle_0({0.5, 1.0, 1.5}, 1.0);
    tracker.particleStart(&particle_0, 1.0);
    particle_0.setLocation({1.0, 1.0, 1.0});
    tracker.particleSplat(&particle_0, 2.0);

    SECTION("Test of valid use with one particle"){
        ParticleSimulation::pMapEntry entry_0 = tracker.get(&particle_0);

        CHECK(entry_0.startTime == 1.0);
        CHECK(entry_0.startLocation == Core::Vector(0.5, 1.0, 1.5));
        CHECK(entry_0.splatTime == 2.0);
        CHECK(entry_0.splatLocation == Core::Vector(1.0, 1.0, 1.0));
        CHECK(entry_0.globalIndex == 0);
        CHECK(particle_0.getIntegerAttribute("global index") == 0);
    }

    SECTION("Test of valid usage with two particles"){
        BTree::Particle particle_1({1.5, 2.0, 2.5}, 1.0);
        tracker.particleStart(&particle_1, 1.1);
        particle_1.setLocation({2.0, 2.0, 2.0});
        tracker.particleSplat(&particle_1, 2.1);

        ParticleSimulation::pMapEntry entry_1 = tracker.get(&particle_1);

        CHECK(entry_1.startTime == 1.1);
        CHECK(entry_1.startLocation == Core::Vector(1.5, 2.0, 2.5));
        CHECK(entry_1.splatTime == 2.1);
        CHECK(entry_1.splatLocation == Core::Vector(2.0, 2.0, 2.0));
        CHECK(entry_1.globalIndex == 1);
        CHECK(particle_1.getIntegerAttribute("global index") == 1);
    }

    SECTION("Test double insert"){
        REQUIRE_THROWS_AS(
                tracker.particleStart(&particle_0, 1.2),
                std::invalid_argument);
    }
}