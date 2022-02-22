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

 Testing of utility functions for particle simulations

 ****************************/

#include "PSim_util.hpp"
#include "BTree_particle.hpp"
#include "Core_vector.hpp"
#include "catch.hpp"
#include <cmath>
TEST_CASE( "Test random particle generation","[ParticleSimulation][utilities][random]") {

    SECTION( "Randomly generated ions on a cylinder wall should be on that cylinder wall"){
        unsigned int nions = 100;
        double cylinderRadius = 10;
        double cylinderLength = 20;
        std::vector<std::unique_ptr<BTree::Particle>> ions =
                ParticleSimulation::util::prepareIonsOnCylinderWalls(nions, 1,cylinderRadius,cylinderLength);

        bool invalidIonFound = false;
        for (std::size_t i=0; i<nions; ++i){
            Core::Vector pos = ions[i]->getLocation();
            double ionRadius = std::sqrt(pos.x()*pos.x() + pos.y()*pos.y());

            if ( !(ionRadius == Approx(cylinderRadius)) ){
                invalidIonFound = true;
            }

            if (!(pos.z() >= -cylinderLength && pos.z() <= cylinderLength)){
                invalidIonFound = true;
            }
        }
        CHECK(! invalidIonFound);
    }
}

TEST_CASE("Test utility functions","[ParticleSimulation][utilities]") {
    SECTION("Probing electric force should at least produce the right spatial positions and result length"){
        std::vector<BTree::Particle> particles;
        int nIons = 10;
        for (int i=0; i<nIons; ++i){
            particles.push_back(BTree::Particle(Core::Vector(1.0,1.0,i*0.1),1.0));
        }

        std::vector<std::tuple<double,double,Core::Vector>> result =
                ParticleSimulation::util::probeForces(particles,ParticleSimulation::Plane::XZ,10,10,0.0,0.0,2.0,10.0,1.0);

        CHECK(result.size() == 100);
        auto lastrow = result.back();
        CHECK(std::get<0>(lastrow) == Approx(1.8));
        CHECK(std::get<1>(lastrow) == Approx(9.0));
    }
}