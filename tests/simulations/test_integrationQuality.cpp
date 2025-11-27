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

template<class INTEGRATOR>
void runIntegrator(std::size_t timeSteps, double dt, double forceConstant, double maxAcceptedPos, std::size_t writePeriod, std::string fileName) {

    auto accelerationFunction =
        [forceConstant](
                Core::Particle *particle, int /*particleIndex*/,
                SpaceCharge::FieldCalculator& /*fieldCalculator*/, double /*time*/, int /*timestep*/) -> Core::Vector{

            Core::Vector force(-forceConstant*particle->getLocation().x(), 0.0, 0.0);
            Core::Vector accel(force / particle->getMass());
            //std::cout << "Acceleration: " << accel << std::endl;
            return (accel);
    };

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

    INTEGRATOR integrator(particlePtrs, accelerationFunction, postTimeStepFunction);
    integrator.run(timeSteps,dt);

    double max = *std::max_element(xPositions.begin(), xPositions.end());
    double min = *std::min_element(xPositions.begin(), xPositions.end());
    CHECK(max < maxAcceptedPos);
    CHECK(min > -maxAcceptedPos);
    std::cout << "Max value: " << max << std::endl;
    std::cout << "Min value: " << min << std::endl;

}

TEST_CASE("Compare integration quality with parabolic profile", "[Simulation]") {
    unsigned int timeSteps = 100000;
    unsigned int writePeriod = 100;
    double dt = 4e-5;
    double forceConstant = 1.0e-20;

    //runIntegrator<Integration::VerletIntegrator>(timeSteps, dt, forceConstant, 0.15, writePeriod,"integration_test_verlet_serial.txt");
    //runIntegrator<Integration::ParallelVerletIntegrator>(timeSteps, dt, forceConstant, 0.15, writePeriod,"integration_test_verlet_parallel.txt");
}