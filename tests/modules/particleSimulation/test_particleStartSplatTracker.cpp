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
#include <memory>


TEST_CASE("TestParticleStartSplatTracker", "[ParticleSimulation][ParticleStartSplatTracker][particle trackers]") {

    SECTION("Basic StartSplatTracker tests") {
        ParticleSimulation::ParticleStartSplatTracker tracker;

        BTree::Particle particle_0({0.5, 1.0, 1.5}, 1.0);
        tracker.particleStart(&particle_0, 1.0);
        particle_0.setLocation({1.0, 1.0, 1.0});
        tracker.particleSplat(&particle_0, 2.0);

        SECTION("Test of valid use with one particle") {
            ParticleSimulation::ParticleStartSplatTracker::pMapEntry entry_0 = tracker.get(&particle_0);

            CHECK(entry_0.startTime==1.0);
            CHECK(entry_0.startLocation==Core::Vector(0.5, 1.0, 1.5));
            CHECK(entry_0.splatTime==2.0);
            CHECK(entry_0.splatLocation==Core::Vector(1.0, 1.0, 1.0));
            CHECK(entry_0.globalIndex==0);
            CHECK(particle_0.getIntegerAttribute("global index")==0);
        }

        SECTION("Test of valid usage with two particles") {
            BTree::Particle particle_1({1.5, 2.0, 2.5}, 1.0);
            tracker.particleStart(&particle_1, 1.1);
            particle_1.setLocation({2.0, 2.0, 2.0});
            tracker.particleSplat(&particle_1, 2.1);

            ParticleSimulation::ParticleStartSplatTracker::pMapEntry entry_1 = tracker.get(&particle_1);

            CHECK(entry_1.startTime==1.1);
            CHECK(entry_1.startLocation==Core::Vector(1.5, 2.0, 2.5));
            CHECK(entry_1.splatTime==2.1);
            CHECK(entry_1.splatLocation==Core::Vector(2.0, 2.0, 2.0));
            CHECK(entry_1.globalIndex==1);
            CHECK(particle_1.getIntegerAttribute("global index")==1);
        }

        SECTION("Test double insert") {
            REQUIRE_THROWS_AS(
                    tracker.particleStart(&particle_0, 1.2),
                    std::invalid_argument);
        }

        SECTION("Test splat of not started particle") {
            BTree::Particle particle_3({3.0, 3.0, 3-0}, 1.0);
            REQUIRE_THROWS_AS(
                    tracker.particleSplat(&particle_3, 1.1),
                    std::invalid_argument);
        }
    }

    SECTION("Tests with large particle cloud"){

        int nParticles = 100;
        ParticleSimulation::ParticleStartSplatTracker tracker;
        std::vector<BTree::uniquePartPtr>particles;
        std::vector<BTree::Particle*>particlesPtrs;

        for (int i=0; i<nParticles; ++i){
            double timeOfBirth = i*0.01;
            BTree::uniquePartPtr particle = std::make_unique<BTree::Particle>(
                    Core::Vector(i*1.0, 1.0, 1.0),
                    Core::Vector(0.0, 0.0, 0.0),
                    1.0,
                    100.0,
                    timeOfBirth);
            particlesPtrs.emplace_back(particle.get());
            particles.emplace_back(std::move(particle));
            tracker.particleStart(particlesPtrs.at(i), timeOfBirth);
        }

        for (int i=nParticles-10; i>=0; --i){ //keep deliberately the last 10 particles non splatted...
            double splatTime = (2*nParticles)*0.01-(i*0.01);
            particles.at(i)->setLocation({1.0, i*1.0, 1.0});
            particles.at(i)->setSplatTime(splatTime);
            tracker.particleSplat(particlesPtrs.at(i), splatTime);
        }

        tracker.sortStartSplatData();
        std::vector<ParticleSimulation::ParticleStartSplatTracker::pMapEntry> pData = tracker.getStartSplatData();
        CHECK(pData[10].globalIndex == 10);
        CHECK(pData[10].startLocation == Core::Vector(10*1.0, 1.0, 1.0));
        CHECK(pData[10].splatLocation == particles[10]->getLocation());
        CHECK(pData[10].startTime == particles[10]->getTimeOfBirth());
        CHECK(pData[10].splatTime == particles[10]->getSplatTime());
        CHECK(pData[10].splatLocation == particles[10]->getLocation());
        CHECK(pData[10].state == ParticleSimulation::ParticleStartSplatTracker::SPLATTED);

        CHECK(pData[25].globalIndex == 25);

        std::vector<double> startTimes = tracker.getStartTimes();
        std::vector<double> splatTimes = tracker.getSplatTimes();
        std::vector<Core::Vector> startLocations = tracker.getStartLocations();
        std::vector<Core::Vector> splatLocations = tracker.getSplatLocations();

        CHECK(startTimes[25] == particles[25]->getTimeOfBirth());
        CHECK(splatTimes[25] == particles[25]->getSplatTime());
        CHECK(startLocations[25] == Core::Vector(25*1.0, 1.0, 1.0));
        CHECK(splatLocations[25] == particles[25]->getLocation());

        CHECK(pData[nParticles-1].state == ParticleSimulation::ParticleStartSplatTracker::STARTED);
        CHECK(pData[nParticles-1].splatTime == 0.0); //since the particle is not splatted
    }
}