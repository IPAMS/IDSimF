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
#include "test_util.hpp"
#include "PSim_particleStartSplatTracker.hpp"
#include "PSim_parallelVerletIntegrator.hpp"
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

        SECTION("Particles should be restartable"){
            BTree::Particle particle_1({1.5, 2.0, 2.5}, 1.0);
            tracker.particleStart(&particle_1, 1.1);
            particle_1.setLocation({2.0, 2.0, 2.0});
            tracker.particleSplat(&particle_1, 2.1);

            BTree::Particle particle_2({-1.5, -2.0, -2.5}, 1.0);
            tracker.particleStart(&particle_2, 1.1);

            CHECK(particle_2.getIntegerAttribute("global index")==2);

            Core::Vector newLocation({-1.5, -2.0, -3.0});
            tracker.particleRestart(&particle_2, particle_2.getLocation(), newLocation, 1.5);
            particle_2.setLocation(newLocation);
            CHECK(particle_2.getIntegerAttribute("global index")==3);

            Core::Vector newLocation2({-1.5, -2.0, -3.5});
            tracker.particleRestart(&particle_2, particle_2.getLocation(), newLocation2, 1.6);
            particle_2.setLocation(newLocation2);
            CHECK(particle_2.getIntegerAttribute("global index")==4);

            tracker.sortStartSplatData();

            std::vector<double> startTimes = tracker.getStartTimes();
            std::vector<double> splatTimes = tracker.getSplatTimes();
            std::vector<Core::Vector> startLocations = tracker.getStartLocations();
            std::vector<Core::Vector> splatLocations = tracker.getSplatLocations();
            std::vector<int> status = tracker.getSplatState();

            CHECK(status[0] == ParticleSimulation::ParticleStartSplatTracker::SPLATTED);
            CHECK(status[1] == ParticleSimulation::ParticleStartSplatTracker::SPLATTED);
            CHECK(status[2] == ParticleSimulation::ParticleStartSplatTracker::SPLATTED_AND_RESTARTED);
            CHECK(status[3] == ParticleSimulation::ParticleStartSplatTracker::SPLATTED_AND_RESTARTED);
            CHECK(status[4] == ParticleSimulation::ParticleStartSplatTracker::RESTARTED);

            CHECK(startTimes == std::vector<double>{1.0, 1.1, 1.1, 1.5, 1.6});
            CHECK(splatTimes == std::vector<double>{2.0, 2.1, 1.5, 1.6, 0.0});

            CHECK(startLocations[3] == Core::Vector{-1.5, -2.0, -3.0});
            CHECK(splatLocations[3] == Core::Vector{-1.5, -2.0, -3.0});
            CHECK(startLocations[4] == Core::Vector{-1.5, -2.0, -3.5});
        }


        SECTION("Double insert should throw") {
            REQUIRE_THROWS_AS(
                    tracker.particleStart(&particle_0, 1.2),
                    std::invalid_argument);
        }

        SECTION("Splat of not started particle should throw") {
            BTree::Particle particle_3({3.0, 3.0, 3-0}, 1.0);
            REQUIRE_THROWS_AS(
                    tracker.particleSplat(&particle_3, 1.1),
                    std::invalid_argument);
        }
    }

    SECTION("Tests with large particle cloud"){

        unsigned int nParticles = 100;
        ParticleSimulation::ParticleStartSplatTracker tracker;
        std::vector<BTree::uniquePartPtr>particles;
        std::vector<BTree::Particle*>particlesPtrs;

        for (unsigned int i=0; i<nParticles; ++i){
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

        for (std::size_t i=nParticles-10; i>0; --i){ //keep deliberately the last 10 particles non splatted...
            double i_d = static_cast<double>(i);
            double splatTime = (2*nParticles)*0.01-(i_d*0.01);
            particles.at(i-1)->setLocation({1.0, (i_d-1)*1.0, 1.0});
            particles.at(i-1)->setSplatTime(splatTime);
            tracker.particleSplat(particlesPtrs.at(i-1), splatTime);
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


    SECTION("Tests with actual particle trajectory integration"){

        ParticleSimulation::ParticleStartSplatTracker tracker;

        // init variables / functions
        double ionAcceleration = 5.0; //((1000V / 100mm) * elementary charge) / 100 amu = 9.64e9 m/s^2
        double dt = 1e-4;


        auto accelerationFct = [ionAcceleration](BTree::Particle* /*particle*/, int /*particleIndex*/, BTree::ParallelTree& /*tree*/,
                                                 double /*time*/, int /*timestep*/){
            Core::Vector result(ionAcceleration, 0, ionAcceleration * 0.5);
            return (result);
        };

        double nParticles = 10;
        unsigned int timeSteps = 100;

        std::vector<BTree::uniquePartPtr>particles;
        std::vector<BTree::Particle*>particlesPtrs;

        double timeOfBirth = 0.0;
        double yPos = 0.0;
        for (int i=0; i<nParticles; ++i){
            BTree::uniquePartPtr particle = std::make_unique<BTree::Particle>(
                    Core::Vector(0.0, yPos, 0.0),
                    Core::Vector(0.0,0.0,0.0),
                    1.0,
                    100.0,
                    timeOfBirth);
            particlesPtrs.push_back(particle.get());
            particles.push_back(std::move(particle));
            yPos += 0.01;
            timeOfBirth += dt*0.5;
        }


        int nParticlesStartMonitored = 0;
        auto particleStartMonitoringFct = [&nParticlesStartMonitored, &tracker] (BTree::Particle* particle, double time){
            nParticlesStartMonitored++;
            tracker.particleStart(particle, time);
        };

        double zMax = 1e-4;
        auto otherActionsFct = [zMax, &tracker] (
                Core::Vector& /*newPartPos*/, BTree::Particle* particle,
                int /*particleIndex*/, BTree::ParallelTree& /*tree*/, double time, int /*timestep*/){
            if (particle->getLocation().z() > zMax){
                particle->setActive(false);
                particle->setSplatTime(time);
                tracker.particleSplat(particle, time);
            }
        };

        ParticleSimulation::ParallelVerletIntegrator verletIntegrator(
                particlesPtrs, accelerationFct, nullptr, otherActionsFct, particleStartMonitoringFct);

        verletIntegrator.run(timeSteps, dt);

        tracker.sortStartSplatData();

        std::vector<double> startTimes = tracker.getStartTimes();
        std::vector<double> splatTimes = tracker.getSplatTimes();
        std::vector<Core::Vector> startLocations = tracker.getStartLocations();
        std::vector<Core::Vector> splatLocations = tracker.getSplatLocations();
        std::vector<int> states = tracker.getSplatState();

        CHECK(states[0] == 2);
        CHECK(states[9] == 2);
        CHECK(startTimes[9] == 10*dt*0.5);
        CHECK(splatTimes[9] == Approx(0.0095));
        CHECK(vectorApproxCompare(splatLocations[9], Core::Vector(0.00020025, 0.09, 0.000100125)) == vectorsApproxEqual);
        CHECK(vectorApproxCompare(startLocations[9], Core::Vector(0.0, 0.09, 0.0)) == vectorsApproxEqual);
    }
}