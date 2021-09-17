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
 test_util.cpp

 Testing of utility functions in RS

 ****************************/

#include "RS_util.hpp"
#include "Core_randomGenerators.hpp"
#include "catch.hpp"
#include <numeric>
#include <vector>
#include <RS_constants.hpp>


TEST_CASE("Random sampled Maxwell-Boltzmann velocity vectors should be correct", "[RS][utility][random]") {
    //This test is a statistical test: use a real random number generator:
    Core::globalRandomGeneratorPool = std::make_unique<Core::RandomGeneratorPool>();
    unsigned int nSamples = 20000;

    double temperature = 298;
    double mass_amu = 19;

    std::vector<double> velocities(nSamples);
    std::generate(velocities.begin(), velocities.end(),
    [temperature,mass_amu](){
        return RS::util::maxwellBoltzmannRandomVelocity(temperature, mass_amu).magnitude();
    } );

    double meanVelocity = std::accumulate(velocities.begin(), velocities.end(), 0.0) / velocities.size();
    double analyticalMeanVelocity = std::sqrt( (8.0*RS::K_BOLTZMANN*temperature)/(M_PI*mass_amu*RS::KG_PER_AMU));
    CHECK(meanVelocity == Approx(analyticalMeanVelocity).margin(5));
}