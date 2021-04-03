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

    SECTION("Test of valid usage with single particle"){
        Core::Vector startLocation(0.5, 1.0, 1.5);
        BTree::Particle particle(startLocation, 1.0);

        tracker.particleStart(&particle, 1.0);

        Core::Vector stopLocation(2.0, 2.0, 2.0);
        particle.setLocation(stopLocation);

        tracker.particleSplat(&particle, 2.0);

        ParticleSimulation::pMapEntry entry = tracker.get(&particle);

        CHECK(entry.startTime == 1.0);
        CHECK(entry.startLocation == startLocation);
        CHECK(entry.splatTime == 2.0);
        CHECK(entry.splatLocation == stopLocation);
    }
}