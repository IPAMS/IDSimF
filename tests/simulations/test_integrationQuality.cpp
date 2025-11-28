/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2025 - Physical and Theoretical Chemistry /
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
 test_integrationQuality.cpp

 Test integration quality of numerical integrators with simple parabolic potential

 ****************************/

#include "Core_particle.hpp"
#include "Integration_abstractTimeIntegrator.hpp"
#include "Integration_verletIntegrator.hpp"
#include "Integration_parallelVerletIntegrator.hpp"
#include "Integration_fullSumVerletIntegrator.hpp"
#include "FileIO_scalar_writer.hpp"
#include "PSim_util.hpp"
#include "catch.hpp"
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <omp.h>

#include "Integration_fmmIntegrator.hpp"
#include "Integration_fullSumRK4Integrator.hpp"
#ifdef WITH_FMM_3d
    #include "FMM3D_fmmSolver.hpp"
#endif

#ifdef WITH_EXAFMMT
    #include "ExaFMMt_fmmSolver.hpp"
#else
    namespace ExaFMMt{
        class FMMSolver;
    }
#endif

#include "Integration_parallelRK4Integrator.hpp"

template<class INTEGRATOR>
void runIntegrator(std::size_t timeSteps, double dt, double forceConstant, double maxAcceptedPos, std::size_t writePeriod, std::string fileName) {

    std::vector<double> xPositions;
    FileIO::Scalar_writer fileWriter(fileName);
    auto postTimeStepFunction =
        [&xPositions, &fileWriter, writePeriod](
            Integration::AbstractTimeIntegrator* /*integrator*/,
            std::vector<Core::Particle*>& particles,
            double time,
            unsigned int timestep,
            bool /*lastTimestep*/) -> void
        {
            xPositions.push_back(particles[0]->getLocation().x());
            if (timestep % writePeriod == 0) {
                fileWriter.writeTimestep(
                    {particles[0]->getLocation().x(), particles[0]->getVelocity().x()},
                    time);
            }
        };
    std::vector<std::unique_ptr<Core::Particle>> particles;
    std::vector<Core::Particle*>particlePtrs;

    double posx = 0.1;
    std::unique_ptr<Core::Particle> newIon = std::make_unique<Core::Particle>(Core::Vector(posx, 0, 0), 1.0);
    newIon->setMassAMU(100);
    particlePtrs.push_back(newIon.get());
    particles.push_back(std::move(newIon));

    if constexpr (
            std::is_same_v<INTEGRATOR, Integration::ParallelRK4Integrator> ||
            std::is_same_v<INTEGRATOR, Integration::FullSumRK4Integrator>
            ) {
        auto multistepAccelerationFunction =
                    [forceConstant](
                            Core::Particle* particle, Core::Vector position, Core::Vector /*velocity*/,
                            double /*time*/, unsigned int /*timestep*/) -> Core::Vector{

                        Core::Vector force(-forceConstant*position.x(), 0.0, 0.0);
                        Core::Vector accel(force / particle->getMass());
                        return (accel);
        };

        auto spaceChargeAccelerationFct = [](Core::Particle* /*particle*/, int /*particleIndex*/, SpaceCharge::FieldCalculator& /*tree*/,
                                             double /*time*/, int /*timestep*/){
            Core::Vector result(0, 0, 0);
            return (result);
        };

        INTEGRATOR integrator(particlePtrs, multistepAccelerationFunction, spaceChargeAccelerationFct, postTimeStepFunction);
        integrator.run(timeSteps,dt);
    }
    else {
        auto accelerationFunction =
            [forceConstant](
                    Core::Particle *particle, int /*particleIndex*/,
                    SpaceCharge::FieldCalculator& /*fieldCalculator*/, double /*time*/, int /*timestep*/) -> Core::Vector{

                Core::Vector force(-forceConstant*particle->getLocation().x(), 0.0, 0.0);
                Core::Vector accel(force / particle->getMass());
                return (accel);
        };

        INTEGRATOR integrator(particlePtrs, accelerationFunction, postTimeStepFunction);
        integrator.run(timeSteps,dt);
    }

    double max = *std::max_element(xPositions.begin(), xPositions.end());
    double min = *std::min_element(xPositions.begin(), xPositions.end());
    CHECK(max < maxAcceptedPos);
    CHECK(min > -maxAcceptedPos);
    //std::cout << "Max value: " << max << std::endl;
    //std::cout << "Min value: " << min << std::endl;

}

TEST_CASE("Compare integration quality with parabolic profile", "[Simulation]") {
    unsigned int timeSteps = 100000;
    unsigned int writePeriod = 100;
    double dt = 4e-5;
    double forceConstant = 1.0e-20;
    double maxAccepted = 0.100002;


    omp_set_num_threads(1);

    SECTION("Test integration quality of serial verlet integrator") {
        runIntegrator<Integration::VerletIntegrator>(timeSteps, dt, forceConstant, maxAccepted, writePeriod,"integration_test_verlet_serial.txt");
    }
    SECTION("Test integration quality of parallel verlet integrator") {
        runIntegrator<Integration::ParallelVerletIntegrator>(timeSteps, dt, forceConstant, maxAccepted, writePeriod,"integration_test_verlet_parallel.txt");
    }
    SECTION("Test integration quality of full sum verlet integrator") {
        runIntegrator<Integration::FullSumVerletIntegrator>(timeSteps, dt, forceConstant, maxAccepted, writePeriod,"integration_test_verlet_full_sum.txt");
    }
    SECTION("Test integration quality of parallel runge kutta integrator") {
        runIntegrator<Integration::ParallelRK4Integrator>(timeSteps, dt, forceConstant, maxAccepted, writePeriod,"integration_test_RK4_parallel.txt");
    }
    SECTION("Test integration quality of full sum runge kutta integrator") {
        runIntegrator<Integration::FullSumRK4Integrator>(timeSteps, dt, forceConstant, maxAccepted, writePeriod,"integration_test_RK4_full_sum.txt");
    }

#ifdef WITH_EXAFMMT
    SECTION("Test integration quality of EXAFMM verlet integrator") {
        runIntegrator<Integration::FMMVerletIntegrator<ExaFMMt::FMMSolver>>(5000, dt, forceConstant, maxAccepted, writePeriod,"integration_test_verlet_EXAFMM.txt");
    }
#endif

#ifdef WITH_FMM_3d
    SECTION("Test integration quality of FMM3D verlet integrator") {
        runIntegrator<Integration::FMMVerletIntegrator<FMM3D::FMMSolver>>(5000, dt, forceConstant, maxAccepted, writePeriod,"integration_test_verlet_EXAFMM.txt");
    }
#endif

}
