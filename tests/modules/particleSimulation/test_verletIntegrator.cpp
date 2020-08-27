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
 test_verletIntegrator.cpp

 Testing of serial version of the verlet integrator

 ****************************/

#include "PSim_verletIntegrator.hpp"
#include "CollisionModel_EmptyCollisionModel.hpp"
#include "Core_vector.hpp"
#include "BTree_particle.hpp"
#include "Core_randomGenerators.hpp"
#include "RS_Substance.hpp"
#include "RS_Simulation.hpp"
#include "catch.hpp"

using sMap = std::map<RS::Substance,int>;
using sPair= sMap::value_type;

TEST_CASE("Test serial verlet integrator", "[ParticleSimulation][VerletIntegrator][trajectory integration]") {

    //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
    Core::globalRandomGenerator = std::make_unique<Core::TestRandomGenerator>();

    double ionAcceleration = 10.0; //((1000V / 100mm) * elementary charge) / 100 amu = 9.64e9 m/s^2

    auto accelerationFct = [ionAcceleration](BTree::Particle *particle, int particleIndex, BTree::Tree &tree,
                                             double time, int timestep){
        Core::Vector result(ionAcceleration, 0, ionAcceleration * 0.5);
        return (result);
    };

    auto timestepWriteFct = [](std::vector<BTree::Particle*>& particles, BTree::Tree& tree, double time, int timestep,
                               bool lastTimestep){};

    auto otherActionsFct = [] (
            Core::Vector& newPartPos,BTree::Particle* particle,
            int particleIndex, BTree::Tree& tree, double time,int timestep){
    };

    SECTION( "Verlet integrator should be working with deferred particle addition") {

        CollisionModel::EmptyCollisionModel collisionModel;

        ParticleSimulation::VerletIntegrator verletIntegrator(
                accelerationFct, timestepWriteFct, otherActionsFct,
                collisionModel);

        //should not crash / throw without particles
        double dt = 1e-4;
        REQUIRE_NOTHROW(verletIntegrator.runSingleStep(dt));

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


    SECTION( "Verlet integrator should be able to integrate correctly non reactive particles") {
        //Test with verlet integration:
        double nParticles = 10;
        double dt = 1e-4;
        double timeSteps = 50;

        std::vector<BTree::uniquePartPtr>particles;
        std::vector<BTree::Particle*>particlesPtrs;

        double yPos = 0;
        for (int i=0; i<nParticles; i++){
            BTree::uniquePartPtr particle = std::make_unique<BTree::Particle>(
                    Core::Vector(0.0,yPos,0.0),
                    Core::Vector(0.0,0.0,0.0),
                    1.0,
                    100.0);
            particlesPtrs.push_back(particle.get());
            particles.push_back(std::move(particle));
            yPos = yPos+0.01;
        }

        CollisionModel::EmptyCollisionModel collisionModel;

        ParticleSimulation::VerletIntegrator verletIntegrator(
                particlesPtrs,
                accelerationFct, timestepWriteFct, otherActionsFct,
                collisionModel);

        verletIntegrator.runSingleStep(dt);
        verletIntegrator.run(timeSteps,dt);

        for (int i=0; i<nParticles; i++){
            Core::Vector ionPos = particles[i]-> getLocation();
            REQUIRE(Approx(ionPos.x()).epsilon(1e-6) == 0.0001275);
            REQUIRE(Approx(ionPos.y()).epsilon(1e-7) == i*0.01);
            REQUIRE(Approx(ionPos.z()).epsilon(1e-7) == 6.375e-05);
        }
    }

    SECTION( "Verlet integrator should work with reactive particles and reactions") {

        //prepare data structures:
        double nParticles = 2;
        double dt = 1e-4;
        int timeSteps = 50;

        auto accelerationFctReactive = [ionAcceleration] (BTree::Particle* particle, int particleIndex,
                BTree::Tree& tree, double time,int timestep){

            Core::Vector result(
                    ionAcceleration * particle->getCharge() / particle->getMass(),
                    0,0);
            return(result);
        };

        //prepare reaction system and particles:
        RS::ConfigFileParser parser = RS::ConfigFileParser();
        RS::Simulation rsSim = RS::Simulation(parser.parseFile("RS_verlet_test.conf"));
        RS::SimulationConfiguration* simConf = rsSim.simulationConfiguration();
        RS::Substance* ed_1 = simConf->substance(0);
        RS::Substance* ed_2 = simConf->substance(1);

        std::vector<uniqueReactivePartPtr>particles;
        std::vector<BTree::Particle*>particlesPtrs; // pointers on the raw particles for trajectory integrator
        double yPos = 0;
        for (int i=0; i<nParticles; i++){
            uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(
                    ed_1,
                    Core::Vector(0.1,yPos,0.0));
            particlesPtrs.push_back(particle.get());
            particles.push_back(std::move(particle));

            particle = std::make_unique<RS::ReactiveParticle>(
                    ed_2,
                    Core::Vector(-0.1,yPos,0.0));
            particlesPtrs.push_back(particle.get());
            particles.push_back(std::move(particle));

            yPos = yPos+0.01;
        }
        for (int i=0; i<particles.size(); i++){
            rsSim.addParticle(particles[i].get(),i);
        }

        // run trajectory integration without reactions:
        CollisionModel::EmptyCollisionModel collisionModel;

        ParticleSimulation::VerletIntegrator verletIntegrator(
                particlesPtrs,
                accelerationFctReactive, timestepWriteFct, otherActionsFct,
                collisionModel);

        verletIntegrator.run(timeSteps,dt);

        Core::Vector ionPos = particles[0]-> getLocation();
        REQUIRE(Approx(ionPos.x()) == 118.298);
        REQUIRE(Approx(ionPos.y()) == 0.0);
        REQUIRE(Approx(ionPos.z()) == 0.0);

        ionPos = particles[3]-> getLocation();
        REQUIRE(Approx(ionPos.x()) == -59.1991);
        REQUIRE(Approx(ionPos.y()) == 0.01);
        REQUIRE(Approx(ionPos.z()) == 0.0);

        //provoke chemical reaction:
        RS::ReactionConditions reactionConditions = RS::ReactionConditions();
        reactionConditions.temperature = 298;
        reactionConditions.pressure = 100000.0;

        for (int ts=0; ts<timeSteps; ts++){
            for (int i=0; i<particles.size(); i++){
                rsSim.react(i, reactionConditions, dt);
            }
            rsSim.advanceTimestep(dt);
            verletIntegrator.runSingleStep(dt);
        }

        ionPos = particles[0]-> getLocation();
        REQUIRE(Approx(ionPos.x()) == 241.321);
        REQUIRE(Approx(ionPos.y()) == 0.0);
        REQUIRE(Approx(ionPos.z()) == 0.0);
        Core::Vector ionVelo = particles[0]-> getVelocity();
        REQUIRE(Approx(ionVelo.x()) == 482.442);
        REQUIRE(Approx(ionVelo.y()) == 0.0);
        REQUIRE(Approx(ionVelo.z()) == 0.0);


        ionPos = particles[3]-> getLocation();
        REQUIRE(Approx(ionPos.x()) == -61.6113);
        REQUIRE(Approx(ionPos.y()) == 0.01);
        REQUIRE(Approx(ionPos.z()) == 0.0);
        ionVelo = particles[3]-> getVelocity();
        REQUIRE(Approx(ionVelo.x()) == 23639.6);
        REQUIRE(Approx(ionVelo.y()) == 0.0);
        REQUIRE(Approx(ionVelo.z()) == 0.0);

        REQUIRE(particles[0]->getSpecies() == simConf->substance(2));
        REQUIRE(particles[3]->getSpecies() == simConf->substance(3));

        rsSim.printConcentrations();
    }
}
