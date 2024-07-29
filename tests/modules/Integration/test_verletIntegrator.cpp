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

#include "Integration_verletIntegrator.hpp"
#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include "Core_randomGenerators.hpp"
#include "RS_Substance.hpp"
#include "RS_Simulation.hpp"
#include "catch.hpp"

using sMap = std::map<RS::Substance,int>;
using sPair= sMap::value_type;

TEST_CASE("Test serial verlet integrator", "[ParticleSimulation][VerletIntegrator][trajectory integration]") {

    //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
    Core::globalRandomGeneratorPool = std::make_unique<Core::TestRandomGeneratorPool>();

    double ionAcceleration = 10.0; //((1000V / 100mm) * elementary charge) / 100 amu = 9.64e9 m/s^2

    auto accelerationFct = [ionAcceleration](Core::Particle* /*particle*/, int /*particleIndex*/, SpaceCharge::FieldCalculator& /*tree*/,
                                             double /*time*/, int /*timestep*/){
        Core::Vector result(ionAcceleration, 0, ionAcceleration * 0.5);
        return (result);
    };

    SECTION( "Verlet integrator should be working with deferred particle addition") {

        Integration::VerletIntegrator verletIntegrator(accelerationFct);

        //should not crash / throw without particles
        double dt = 1e-4;
        CHECK_NOTHROW(verletIntegrator.runSingleStep(dt));

        //particles should be addable and integrator should be able to run:
        unsigned int nSteps = 100;
        Core::Particle testParticle1(Core::Vector(0.0, 0.0, 0.0),
                Core::Vector(0.0, 0.0, 0.0),
                1.0, 100.0);

        CHECK_NOTHROW(verletIntegrator.addParticle(&testParticle1));
        CHECK_NOTHROW(verletIntegrator.run(nSteps,dt));

        Core::Vector ionPos1 = testParticle1.getLocation();
        CHECK(Approx(ionPos1.x()).epsilon(1e-6) == 0.000495);
        CHECK(Approx(ionPos1.y()).epsilon(1e-2) == 0.00);
        CHECK(Approx(ionPos1.z()).epsilon(1e-7) == 0.0002475);


        std::cout<<"p 1: "<<testParticle1.getLocation()<<std::endl;

        Core::Particle testParticle2(Core::Vector(0.0, 0.01, 0.0),
                Core::Vector(0.0, 0.0, 0.0),
                1.0, 100.0);

        CHECK_NOTHROW(verletIntegrator.addParticle(&testParticle2));
        CHECK_NOTHROW(verletIntegrator.run(nSteps,dt));

        //verletIntegrator.run(nSteps,dt);

        Core::Vector ionPos2 = testParticle2.getLocation();
        CHECK(Approx(ionPos2.x()).epsilon(1e-6) == 0.000495);
        CHECK(Approx(ionPos2.y()).epsilon(1e-2) == 0.01);
        CHECK(Approx(ionPos2.z()).epsilon(1e-7) == 0.0002475);
    }

    SECTION( "Verlet integrator should be able to integrate correctly non reactive particles with TOB distribution") {
        unsigned int nParticles = 10;
        unsigned int timeSteps = 60;
        double dt = 1e-4;

        //prepare particles with a distributed time of birth (TOB):
        std::vector<Core::uniquePartPtr>particles;
        std::vector<Core::Particle*>particlesPtrs;

        double yPos = 0;
        double lastTime = timeSteps * dt - 4*dt;
        double timeOfBirth = lastTime;
        for (unsigned int i=0; i<nParticles; ++i){
            Core::uniquePartPtr particle = std::make_unique<Core::Particle>(
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

        SECTION("Verlet integrator should run through and should call functions"){

            unsigned int nTimeStepsRecorded = 0;
            auto postTimestepFct = [&nTimeStepsRecorded](
                    Integration::AbstractTimeIntegrator* /*integrator*/, std::vector<Core::Particle*>& /*particles*/, double /*time*/, int /*timestep*/,
                                       bool /*lastTimestep*/){
                nTimeStepsRecorded++;
            };

            int nParticlesTouched = 0;
            auto otherActionsFct = [&nParticlesTouched] (
                    Core::Vector& /*newPartPos*/,Core::Particle* /*particle*/,
                    int /*particleIndex*/, double /*time*/, int /*timestep*/){
                nParticlesTouched++;
            };

            unsigned int nParticlesStartMonitored = 0;
            auto particleStartMonitoringFct = [&nParticlesStartMonitored] (Core::Particle* /*particle*/, double /*time*/){
                nParticlesStartMonitored++;
            };

            Integration::VerletIntegrator verletIntegrator(
                    particlesPtrs,
                    accelerationFct, postTimestepFct, otherActionsFct, particleStartMonitoringFct);

            verletIntegrator.run(timeSteps, dt);

            CHECK(nTimeStepsRecorded == timeSteps + 2);
            CHECK(nParticlesTouched == 58);
            CHECK(nParticlesStartMonitored == nParticles);

            CHECK(verletIntegrator.timeStep() == timeSteps);
            CHECK(verletIntegrator.time() == Approx(timeSteps * dt));

            double endTime = timeSteps * dt;
            for (std::size_t i=0; i<nParticles; ++i){
                Core::Vector ionPos = particles[i]-> getLocation();

                //calculate approximate position according to a pure linear uniform acceleration
                // according to the real time the particles were present in the simulation:
                double diffTime = endTime - (0.5*dt)  - particles[i]-> getTimeOfBirth();
                double xCalculated= 0.5 * ionAcceleration * diffTime * diffTime;
                double zCalculated= 0.5 * xCalculated;

                CHECK(Approx(ionPos.x()).epsilon(0.05) == xCalculated);
                CHECK(Approx(ionPos.y()).epsilon(1e-7) == i*0.01);
                CHECK(Approx(ionPos.z()).epsilon(0.05) == zCalculated);
            }
        }

        SECTION("Verlet integrator should be stoppable"){
            unsigned int terminationTimeStep = 50;

            auto timestepStopFct = [terminationTimeStep](
                    Integration::AbstractTimeIntegrator* integrator, std::vector<Core::Particle*>& /*particles*/, double /*time*/, unsigned int timestep, bool /*lastTimestep*/){
                if (timestep >= terminationTimeStep){
                    integrator->setTerminationState();
                }
            };

            Integration::VerletIntegrator verletIntegrator(
                    particlesPtrs,
                    accelerationFct, timestepStopFct);
            verletIntegrator.run(timeSteps, dt);

            CHECK(verletIntegrator.timeStep() == terminationTimeStep);
        }
    }

    SECTION( "Verlet integrator should work with reactive particles and reactions") {

        //prepare data structures:
        double nParticles = 2;
        double dt = 1e-4;
        unsigned int timeSteps = 50;

        auto accelerationFctReactive = [ionAcceleration] (
                Core::Particle* particle, int /*particleIndex*/,
                SpaceCharge::FieldCalculator& /*fc*/, double /*time*/, int /*timestep*/){
            Core::Vector result(ionAcceleration * particle->getCharge() / particle->getMass(), 0,0);
            return(result);
        };

        //prepare reaction system and particles:
        RS::ConfigFileParser parser = RS::ConfigFileParser();
        RS::Simulation rsSim = RS::Simulation(parser.parseFile("RS_verlet_test.conf"));
        RS::SimulationConfiguration* simConf = rsSim.simulationConfiguration();
        RS::Substance* ed_1 = simConf->substance(0);
        RS::Substance* ed_2 = simConf->substance(1);

        std::vector<uniqueReactivePartPtr>particles;
        std::vector<Core::Particle*>particlesPtrs; // pointers on the raw particles for trajectory integrator
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
        for (std::size_t i=0; i<particles.size(); i++){
            rsSim.addParticle(particles[i].get(), i);
        }

        // run trajectory integration without reactions:

        Integration::VerletIntegrator verletIntegrator(
                particlesPtrs, accelerationFctReactive);

        verletIntegrator.run(timeSteps,dt);

        Core::Vector ionPos = particles[0]-> getLocation();
        CHECK(Approx(ionPos.x()) == 118.298);
        CHECK(Approx(ionPos.y()) == 0.0);
        CHECK(Approx(ionPos.z()) == 0.0);

        ionPos = particles[3]-> getLocation();
        CHECK(Approx(ionPos.x()) == -59.1991);
        CHECK(Approx(ionPos.y()) == 0.01);
        CHECK(Approx(ionPos.z()) == 0.0);

        //provoke chemical reaction:
        RS::ReactionConditions reactionConditions = RS::ReactionConditions();
        reactionConditions.temperature = 298;
        reactionConditions.pressure = 100000.0;

        for (unsigned int ts=0; ts<timeSteps; ts++){
            for (std::size_t i=0; i<particles.size(); i++){
                rsSim.react(i, reactionConditions, dt);
            }
            rsSim.advanceTimestep(dt);
            verletIntegrator.runSingleStep(dt);
        }

        ionPos = particles[0]-> getLocation();
        CHECK(Approx(ionPos.x()) == 241.321);
        CHECK(Approx(ionPos.y()) == 0.0);
        CHECK(Approx(ionPos.z()) == 0.0);
        Core::Vector ionVelo = particles[0]-> getVelocity();
        CHECK(Approx(ionVelo.x()) == 482.442);
        CHECK(Approx(ionVelo.y()) == 0.0);
        CHECK(Approx(ionVelo.z()) == 0.0);

        ionPos = particles[3]-> getLocation();
        CHECK(Approx(ionPos.x()) == -61.6113);
        CHECK(Approx(ionPos.y()) == 0.01);
        CHECK(Approx(ionPos.z()) == 0.0);
        ionVelo = particles[3]-> getVelocity();
        CHECK(Approx(ionVelo.x()) == 23639.6);
        CHECK(Approx(ionVelo.y()) == 0.0);
        CHECK(Approx(ionVelo.z()) == 0.0);

        CHECK(particles[0]->getSpecies() == simConf->substance(2));
        CHECK(particles[3]->getSpecies() == simConf->substance(3));

        rsSim.printConcentrations();
    }
}
