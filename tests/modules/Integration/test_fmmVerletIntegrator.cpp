/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2022 - Physical and Theoretical Chemistry /
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
 test_fmmVerletIntegrator.cpp

 Testing of fmm based verlet trajectory integrator

 ****************************/

#include "Integration_fmmIntegrator.hpp"
#ifdef WITH_FMM_3d
    #include "FMM3D_fmmSolver.hpp"
#endif

#ifdef WITH_EXAFMMT
    #include "ExaFMMt_fmmSolver.hpp"
#endif

#include "SC_generic.hpp"
#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include "catch.hpp"

template<class solverType>
        void runTestSimulation(std::vector<Core::uniquePartPtr> &particles,
                               std::vector<Core::Particle*> &particlesPtrs,
                               std::size_t nParticles, int timeSteps, double dt
        ){
            double ionAcceleration = 10.0; //((1000V / 100mm) * elementary charge) / 100 amu = 9.64e9 m/s^2

            auto accelerationFct = [ionAcceleration](Core::Particle* /*particle*/, int /*particleIndex*/, SpaceCharge::FieldCalculator& /*tree*/,
                    double /*time*/, int /*timestep*/){
                Core::Vector result(ionAcceleration, 0, ionAcceleration * 0.5);
                return (result);
            };

            unsigned int nTimestepsRecorded = 0;
            auto timestepWriteFct = [&nTimestepsRecorded](std::vector<Core::Particle*>& /*particles*/,
                    double /*time*/, int /*timestep*/,
                    bool /*lastTimestep*/) {
                nTimestepsRecorded++;
            };

            unsigned int nParticlesTouched = 0;
            auto otherActionsFct = [&nParticlesTouched](
                    Core::Vector& /*newPartPos*/, Core::Particle* /*particle*/,
                    int /*particleIndex*/, double /*time*/, int /*timestep*/) {

                #pragma omp atomic
                nParticlesTouched++;
            };

            unsigned int nParticlesStartMonitored = 0;
            auto particleStartMonitoringFct = [&nParticlesStartMonitored](Core::Particle* /*particle*/,
                    double /*time*/) {
                nParticlesStartMonitored++;
            };

            Integration::FMMVerletIntegrator<solverType> verletIntegrator(
                    particlesPtrs, accelerationFct, timestepWriteFct, otherActionsFct, particleStartMonitoringFct);

#ifdef WITH_EXAFMMT
            if constexpr (std::is_same_v<solverType, ExaFMMt::FMMSolver>) {
                verletIntegrator.getFMMSolver()->setExpansionOrder(7);
            }
#endif

#ifdef WITH_FMM_3d
            if constexpr (std::is_same_v<solverType, FMM3D::FMMSolver>) {
                verletIntegrator.getFMMSolver()->setRequestedPrecision(0.3e-6);
            }
#endif
            verletIntegrator.run(timeSteps, dt);

            CHECK(verletIntegrator.time()==Approx(timeSteps*dt));
            CHECK(verletIntegrator.timeStep()==timeSteps);

            double endTime = timeSteps*dt;
            for (std::size_t i = 0; i<nParticles; ++i) {
                Core::Vector ionPos = particles[i]->getLocation();

                //calculate approximate position according to a pure linear uniform acceleration
                // according to the real time the particles were present in the simulation:
                double diffTime = endTime-(0.5*dt)-particles[i]->getTimeOfBirth();
                double xCalculated = 0.5*ionAcceleration*diffTime*diffTime;
                double zCalculated = 0.5*xCalculated;

                CHECK(Approx(ionPos.x()).epsilon(0.05)==xCalculated);
                CHECK(Approx(ionPos.y()).epsilon(1e-7)==i*0.01);
                CHECK(Approx(ionPos.z()).epsilon(0.05)==zCalculated);
            }

            CHECK(nTimestepsRecorded==timeSteps+2);
            CHECK(nParticlesTouched==timeSteps*nParticles);
            CHECK(nParticlesStartMonitored==nParticles);
        }



TEST_CASE( "Test fmm verlet integrator", "[Integration][FMMVerletIntegrator][trajectory integration]") {
    // init variables / functions

    double dt = 1e-4;
    unsigned int nParticles = 10;
    unsigned int timeSteps = 60;

    std::vector<Core::uniquePartPtr> particles;
    std::vector<Core::Particle*> particlesPtrs;

    SECTION("FMM Verlet integrator should be able to integrate correctly particles born at start") {

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

#ifdef WITH_FMM_3d
        SECTION("FMM3D FMM integrator should integrate correctly") {
            runTestSimulation<FMM3D::FMMSolver>(particles, particlesPtrs, nParticles, timeSteps, dt);
        }
#endif

#ifdef WITH_EXAFMMT
        SECTION("EXAFMMT integrator should integrate correctly") {
            runTestSimulation<ExaFMMt::FMMSolver>(particles, particlesPtrs, nParticles, timeSteps, dt);
        }
#endif
    }
}
