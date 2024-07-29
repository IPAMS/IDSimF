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

#include "Integration_fullSumRK4Integrator.hpp"
#include "SC_generic.hpp"
#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include "Core_randomGenerators.hpp"
#include "RS_Substance.hpp"
#include "RS_Simulation.hpp"
#include "catch.hpp"

using sMap = std::map<RS::Substance,int>;
using sPair= sMap::value_type;

TEST_CASE( "Test full sum runge kutta 4 integrator", "[ParticleSimulation][ParallelRK4Integrator][trajectory integration]") {
    //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
    Core::globalRandomGeneratorPool = std::make_unique<Core::TestRandomGeneratorPool>();

    // init variables / functions
    double ionAcceleration = 10.0; //((1000V / 100mm) * elementary charge) / 100 amu = 9.64e9 m/s^2
    double spaceChargeAcceleration = 1.0;
    double dt = 1e-4;

    auto spaceChargeAccelerationFct = [spaceChargeAcceleration](Core::Particle* /*particle*/, int /*particleIndex*/, SpaceCharge::FieldCalculator& /*fieldSolver*/,
                                             double /*time*/, int /*timestep*/){
        Core::Vector result(0, spaceChargeAcceleration, 0);
        return (result);
    };

    auto accelerationFct = [ionAcceleration](
            Core::Particle* /*particle*/, Core::Vector /*position*/, Core::Vector /*velocity*/, double /*time*/, unsigned int /*timestep*/){
        Core::Vector result(ionAcceleration, 0, ionAcceleration * 0.5);
        return (result);
    };

    Core::Particle testParticle1(Core::Vector(0.0, 0.0, 0.0),
            Core::Vector(0.0, 0.0, 0.0),
            1.0, 100.0);

    Core::Particle testParticle2(Core::Vector(0.0, 0.01, 0.0),
            Core::Vector(0.0, 0.0, 0.0),
            1.0, 100.0);


    SECTION( "Parallel full sum RK integrator should be working with deferred particle addition") {

        // bare integrator without time step or other actions function
        Integration::FullSumRK4Integrator RK4Integrator(accelerationFct, spaceChargeAccelerationFct);

        //should not crash / throw without particles
        REQUIRE_NOTHROW(RK4Integrator.run(1,dt));

        //particles should be addable and integrator should be able to run:
        unsigned int nSteps = 100;
        unsigned int nStepsLong = 1000;
        REQUIRE_NOTHROW(RK4Integrator.addParticle(&testParticle1));
        REQUIRE_NOTHROW(RK4Integrator.run(nSteps,dt));

        REQUIRE_NOTHROW(RK4Integrator.addParticle(&testParticle2));
        REQUIRE_NOTHROW(RK4Integrator.run(nSteps*2,dt));

        Core::Vector ionPos = testParticle2.getLocation();
        CHECK(Approx(ionPos.x()).epsilon(1e-6) == 0.002);
        CHECK(Approx(ionPos.y()).epsilon(1e-2) == 0.0102);
        CHECK(Approx(ionPos.z()).epsilon(1e-7) == 0.001);

        REQUIRE_NOTHROW(RK4Integrator.run(nStepsLong,dt));

        ionPos = testParticle2.getLocation();
        CHECK(Approx(ionPos.x()).epsilon(1e-6) == 0.072);
        CHECK(Approx(ionPos.y()).epsilon(1e-2) == 0.0172);
        CHECK(Approx(ionPos.z()).epsilon(1e-7) == 0.036);
    }

    SECTION("Tests with particle lists"){

        unsigned int nParticles = 10;
        unsigned int timeSteps = 60;

        std::vector<Core::uniquePartPtr>particles;
        std::vector<Core::Particle*>particlesPtrs;

        SECTION( "Full sum RK4 integrator should be able to integrate correctly non reactive particles with ToB") {

            //prepare particles with a distributed time of birth (TOB):
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

            Integration::FullSumRK4Integrator RK4Integrator(
                    particlesPtrs, accelerationFct, spaceChargeAccelerationFct, nullptr, nullptr, nullptr);

            RK4Integrator.run(timeSteps, dt);

            double endTime = timeSteps * dt;
            for (std::size_t i=0; i<nParticles; ++i){
                Core::Vector ionPos = particles[i]-> getLocation();

                //calculate approximate position according to a pure linear uniform acceleration
                // according to the real time the particles were present in the simulation:
                double diffTime = endTime - particles[i]-> getTimeOfBirth();
                double xCalculated= 0.5 * ionAcceleration * diffTime * diffTime;
                double yDriftDistance = 0.5 * spaceChargeAcceleration * diffTime * diffTime;
                double zCalculated= 0.5 * xCalculated;

                CHECK(Approx(ionPos.x()).epsilon(0.05) == xCalculated);
                CHECK(Approx(ionPos.y()).epsilon(1e-7) == i*0.01 + yDriftDistance);
                CHECK(Approx(ionPos.z()).epsilon(0.05) == zCalculated);
            }
        }

        SECTION( "Full sum RK4 integrator should be able to integrate correctly particles born at start") {

            //prepare particles all born at t=0:
            double yPos = 0;
            double timeOfBirth = 0.0;
            for (unsigned int i = 0; i<nParticles; ++i) {
                Core::uniquePartPtr particle = std::make_unique<Core::Particle>(
                        Core::Vector(0.0, yPos, 0.0),
                        Core::Vector(0.0, 0.0, 0.0),
                        1.0,
                        100.0,
                        timeOfBirth);
                particlesPtrs.push_back(particle.get());
                particles.push_back(std::move(particle));
                yPos += 0.01;
            }

            SECTION("Integration should run through and functions should be called") {

                unsigned int nTimestepsRecorded = 0;
                auto postTimestepFct = [&nTimestepsRecorded](
                        Integration::AbstractTimeIntegrator* /*integrator*/, std::vector<Core::Particle*>& /*particles*/,
                        double /*time*/, int /*timestep*/, bool /*lastTimestep*/){
                    nTimestepsRecorded++;
                };

                unsigned int nParticlesTouched = 0;
                auto otherActionsFct = [&nParticlesTouched] (
                        Core::Vector& /*newPartPos*/, Core::Particle* /*particle*/,
                        int /*particleIndex*/, double /*time*/, int /*timestep*/){
                    nParticlesTouched++;
                };

                unsigned int nParticlesStartMonitored = 0;
                auto particleStartMonitoringFct = [&nParticlesStartMonitored] (Core::Particle* /*particle*/, double /*time*/){
                    nParticlesStartMonitored++;
                };

                Integration::FullSumRK4Integrator RK4Integrator(
                        particlesPtrs, accelerationFct, spaceChargeAccelerationFct, postTimestepFct, otherActionsFct, particleStartMonitoringFct);

                RK4Integrator.run(timeSteps, dt);

                CHECK(RK4Integrator.time() == Approx(timeSteps*dt));
                CHECK(RK4Integrator.timeStep() == timeSteps);

                double endTime = timeSteps*dt;
                for (std::size_t i = 0; i<nParticles; ++i) {
                    Core::Vector ionPos = particles[i]->getLocation();

                    //calculate approximate position according to a pure linear uniform acceleration
                    // according to the real time the particles were present in the simulation:
                    double diffTime = endTime-(0.5*dt)-particles[i]->getTimeOfBirth();
                    double xCalculated = 0.5*ionAcceleration*diffTime*diffTime;
                    double zCalculated = 0.5*xCalculated;

                    CHECK(Approx(ionPos.x()).epsilon(0.05)==xCalculated);
                    //CHECK(Approx(ionPos.y()).epsilon(1e-7)==i*0.01);
                    CHECK(Approx(ionPos.z()).epsilon(0.05)==zCalculated);
                }

                CHECK(nTimestepsRecorded == timeSteps + 2);
                CHECK(nParticlesTouched == timeSteps * nParticles);
                CHECK(nParticlesStartMonitored == nParticles);
            }

            SECTION("Integration should be stoppable") {
                unsigned int terminationTimeStep = 40;

                unsigned int nTimestepsRecorded = 0;
                auto postTimestepFct = [&nTimestepsRecorded, terminationTimeStep](
                        Integration::AbstractTimeIntegrator* integrator,
                        std::vector<Core::Particle*>& /*particles*/,
                        double /*time*/, unsigned int timestep, bool /*lastTimestep*/)
                {
                    nTimestepsRecorded++;
                    if (timestep >= terminationTimeStep) {
                        integrator->setTerminationState();
                    }
                };

                Integration::FullSumRK4Integrator verletIntegrator(
                        particlesPtrs,
                        accelerationFct, spaceChargeAccelerationFct,
                        postTimestepFct);

                verletIntegrator.run(timeSteps, dt);
                CHECK(nTimestepsRecorded == terminationTimeStep+2);
                CHECK(verletIntegrator.timeStep() == terminationTimeStep);
                CHECK(verletIntegrator.time() == Approx(dt*(terminationTimeStep)));
            }
        }
    }
}
