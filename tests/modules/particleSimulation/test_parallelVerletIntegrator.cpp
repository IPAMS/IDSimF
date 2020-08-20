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
 test_parallelVerletIntegrator.cpp

 Testing of parallelized verlet trajectory integrator

 ****************************/

#include "PSim_parallelVerletIntegrator.hpp"
#include "CollisionModel_EmptyCollisionModel.hpp"
#include "Core_vector.hpp"
#include "BTree_particle.hpp"
#include "Core_randomGenerators.hpp"
#include "RS_Substance.hpp"
#include "RS_Simulation.hpp"
#include "catch.hpp"

using sMap = std::map<RS::Substance,int>;
using sPair= sMap::value_type;

TEST_CASE( "Test parallel verlet integrator", "[ParticleSimulation][ParallelVerletIntegrator][trajectory integration]") {

    SECTION( "Parallel Verlet integrator should be working with deferred particle addition") {

        // init empty verlet integrator
        double ionAcceleration = 10.0; //((1000V / 100mm) * elementary charge) / 100 amu = 9.64e9 m/s^2

        auto accelerationFct = [ionAcceleration](BTree::Particle *particle, int particleIndex, BTree::ParallelTree &tree,
                                                 double time, int timestep){
            Core::Vector result(ionAcceleration, 0, ionAcceleration * 0.5);
            return (result);
        };

        auto timestepWriteFct = [](std::vector<BTree::Particle*>& particles, BTree::ParallelTree& tree, double time, int timestep,
                                   bool lastTimestep){};

        auto otherActionsFct = [] (
                Core::Vector& newPartPos,BTree::Particle* particle,
                int particleIndex, BTree::ParallelTree& tree, double time,int timestep){
        };

        CollisionModel::EmptyCollisionModel collisionModel;

        ParticleSimulation::ParallelVerletIntegrator verletIntegrator(
                accelerationFct, timestepWriteFct, otherActionsFct,
                collisionModel);

        //should not crash / throw without particles
        double dt = 1e-4;
        REQUIRE_NOTHROW(verletIntegrator.run(1,dt));

        //particles should be addable and integrator should be able to run:
        int nSteps = 100;
        BTree::Particle testParticle1(Core::Vector(0.0, 0.0, 0.0),
                Core::Vector(0.0, 0.0, 0.0),
                1.0, 100.0);

        REQUIRE_NOTHROW(verletIntegrator.addParticle(&testParticle1));
        REQUIRE_NOTHROW(verletIntegrator.run(nSteps,dt));

        std::cout<<"p 1: "<<testParticle1.getLocation()<<std::endl;

        BTree::Particle testParticle2(Core::Vector(0.0, 0.01, 0.0),
                Core::Vector(0.0, 0.0, 0.0),
                1.0, 100.0);

        REQUIRE_NOTHROW(verletIntegrator.addParticle(&testParticle2));
        REQUIRE_NOTHROW(verletIntegrator.run(nSteps,dt));

        verletIntegrator.run(nSteps,dt);

        Core::Vector ionPos = testParticle2.getLocation();
        REQUIRE(Approx(ionPos.x()).epsilon(1e-6) == 0.00199);
        REQUIRE(Approx(ionPos.y()).epsilon(1e-2) == 0.01);
        REQUIRE(Approx(ionPos.z()).epsilon(1e-7) == 0.000995);
    }

    SECTION( "Parallel Verlet integrator should be able to integrate correctly non reactive particles") {
        //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
        Core::globalRandomGenerator = std::make_unique<Core::TestRandomGenerator>();

        //Test with verlet integration:
        double nParticles = 10;
        double dt = 1e-4;
        double timeSteps = 60;

        double ionAcceleration = 10.0; //((1000V / 100mm) * elementary charge) / 100 amu = 9.64e9 m/s^2

        auto accelerationFct = [ionAcceleration](BTree::Particle *particle, int particleIndex, BTree::ParallelTree &tree,
                                                 double time, int timestep){
            Core::Vector result(ionAcceleration, 0, ionAcceleration * 0.5);
            return (result);
        };

        auto timestepWriteFct = [](std::vector<BTree::Particle*>& particles, BTree::ParallelTree& tree, double time, int timestep,
                                   bool lastTimestep){};

        auto otherActionsFct = [] (
                Core::Vector& newPartPos,BTree::Particle* particle,
                int particleIndex, BTree::ParallelTree& tree, double time,int timestep){
        };


        std::vector<BTree::uniquePartPtr>particles;
        std::vector<BTree::Particle*>particlesPtrs;

        //prepare particles with a distributed time of birth (TOB):
        double yPos = 0;
        double lastTime = timeSteps * dt - 4*dt;
        double timeOfBirth = lastTime;
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
            timeOfBirth -= dt*0.5;
        }


        CollisionModel::EmptyCollisionModel collisionModel;

        ParticleSimulation::ParallelVerletIntegrator verletIntegrator(
                particlesPtrs,
                accelerationFct, timestepWriteFct, otherActionsFct,
                collisionModel);

        verletIntegrator.run(timeSteps,dt);

        double endTime = timeSteps * dt;
        for (int i=0; i<nParticles; ++i){
            Core::Vector ionPos = particles[i]-> getLocation();

            //calculate approximate position according to a pure linear uniform acceleration
            // according to the real time the particles were present in the simulation:
            double diffTime = endTime - (0.5*dt)  - particles[i]-> getTimeOfBirth();
            double xCalculated= 0.5 * ionAcceleration * diffTime * diffTime;
            double zCalculated= 0.5 * xCalculated;

            REQUIRE(Approx(ionPos.x()).epsilon(0.05) == xCalculated);
            REQUIRE(Approx(ionPos.y()).epsilon(1e-7) == i*0.01);
            REQUIRE(Approx(ionPos.z()).epsilon(0.05) == zCalculated);
        }
    }

    SECTION( "Parallel Verlet integrator should be able to integrate correctly particles born at start") {
        //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
        Core::globalRandomGenerator = std::make_unique<Core::TestRandomGenerator>();

        //Test with verlet integration:
        double nParticles = 10;
        double dt = 1e-4;
        double timeSteps = 60;

        double ionAcceleration = 10.0; //((1000V / 100mm) * elementary charge) / 100 amu = 9.64e9 m/s^2

        auto accelerationFct = [ionAcceleration](BTree::Particle* particle, int particleIndex,
                                                 BTree::ParallelTree& tree,
                                                 double time, int timestep) {
            Core::Vector result(ionAcceleration, 0, ionAcceleration*0.5);
            return (result);
        };

        auto timestepWriteFct = [](std::vector<BTree::Particle*>& particles, BTree::ParallelTree& tree, double time,
                                   int timestep,
                                   bool lastTimestep) { };

        auto otherActionsFct = [](
                Core::Vector& newPartPos, BTree::Particle* particle,
                int particleIndex, BTree::ParallelTree& tree, double time, int timestep) {
        };

        std::vector<BTree::uniquePartPtr> particles;
        std::vector<BTree::Particle*> particlesPtrs;

        //prepare particles with a distributed time of birth (TOB):
        double yPos = 0;
        double lastTime = timeSteps*dt-4*dt;
        double timeOfBirth = 0.0;
        for (int i = 0; i<nParticles; ++i) {
            BTree::uniquePartPtr particle = std::make_unique<BTree::Particle>(
                    Core::Vector(0.0, yPos, 0.0),
                    Core::Vector(0.0, 0.0, 0.0),
                    1.0,
                    100.0,
                    timeOfBirth);
            particlesPtrs.push_back(particle.get());
            particles.push_back(std::move(particle));
            yPos += 0.01;
            timeOfBirth -= dt*0.5;
        }

        CollisionModel::EmptyCollisionModel collisionModel;

        ParticleSimulation::ParallelVerletIntegrator verletIntegrator(
                particlesPtrs,
                accelerationFct, timestepWriteFct, otherActionsFct,
                collisionModel);

        verletIntegrator.run(timeSteps, dt);

        double endTime = timeSteps*dt;
        for (int i = 0; i<nParticles; ++i) {
            Core::Vector ionPos = particles[i]->getLocation();

            //calculate approximate position according to a pure linear uniform acceleration
            // according to the real time the particles were present in the simulation:
            double diffTime = endTime-(0.5*dt)-particles[i]->getTimeOfBirth();
            double xCalculated = 0.5*ionAcceleration*diffTime*diffTime;
            double zCalculated = 0.5*xCalculated;

            REQUIRE(Approx(ionPos.x()).epsilon(0.05)==xCalculated);
            REQUIRE(Approx(ionPos.y()).epsilon(1e-7)==i*0.01);
            REQUIRE(Approx(ionPos.z()).epsilon(0.05)==zCalculated);
        }
    }
}
