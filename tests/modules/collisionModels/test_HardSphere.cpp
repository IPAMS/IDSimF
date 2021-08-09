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
 test_HardSphere.cpp

 Testing of hard sphere collision model

 ****************************/

#include "catch.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "Core_constants.hpp"
#include "Core_randomGenerators.hpp"
#include "BTree_particle.hpp"

#include <memory>

TEST_CASE( "Basic test Hard Sphere model", "[CollisionModels][HardSphereModel]") {
    //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
    Core::globalRandomGenerator = std::make_unique<Core::TestRandomGenerator>();

    double diameterHe = CollisionModel::HardSphereModel::DIAMETER_HE;
    BTree::Particle ion = BTree::Particle();
    ion.setDiameter(CollisionModel::HardSphereModel::DIAMETER_HE);
    ion.setVelocity(Core::Vector(100,0,0));
    ion.setMassAMU(28.0);
    int n = 100;

    SECTION("Test without maxwellian approximation"){
        CollisionModel::HardSphereModel hs = CollisionModel::HardSphereModel(
                1.0,298,4.0,diameterHe,false);

        for (int i =0; i<n; i++){
            hs.modifyVelocity(ion, 2e-7);
        }
        Core::Vector ionVelo = ion.getVelocity();

        CHECK(Approx(ionVelo.x()) ==  -27.2881052);
        CHECK(Approx(ionVelo.y()) ==  32.5138876);
        CHECK(Approx(ionVelo.z()) ==  -128.48497);
    }

    SECTION("Test with maxwellian approximation"){
        CollisionModel::HardSphereModel hs = CollisionModel::HardSphereModel(
                1.0,298,4.0,diameterHe,true);

        for (int i =0; i<n; i++){
            hs.modifyVelocity(ion, 2e-7);
        }
        Core::Vector ionVelo = ion.getVelocity();

        CHECK(Approx(ionVelo.x()) ==  147.714697438);
        CHECK(Approx(ionVelo.y()) ==  58.6209442494);
        CHECK(Approx(ionVelo.z()) ==  -91.735975867);
    }
}

TEST_CASE( "Test Hard Sphere model after collision function", "[CollisionModels][HardSphereModel]") {
    //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
    Core::globalRandomGenerator = std::make_unique<Core::TestRandomGenerator>();

    //define a simple after collision function, which is called after a collision has taken place:
    int nTotalCollisions = 0;
    std::vector<double> collisionEnergies;
    auto afterCollisionFct = [&nTotalCollisions, &collisionEnergies](
            RS::CollisionConditions collisionConditions, BTree::Particle &){
        nTotalCollisions++;
        collisionEnergies.push_back(collisionConditions.totalCollisionEnergy * Core::JOULE_TO_EV);
    };
    double diameterHe = CollisionModel::HardSphereModel::DIAMETER_HE;
    CollisionModel::HardSphereModel hs = CollisionModel::HardSphereModel(1.0,298,4.0,diameterHe,afterCollisionFct);
    BTree::Particle ion = BTree::Particle();
    ion.setDiameter(CollisionModel::HardSphereModel::DIAMETER_HE);
    ion.setVelocity(Core::Vector(100,0,0));
    ion.setMassAMU(28.0);
    int n = 60;
    for (int i =0; i<n; i++){
        hs.modifyVelocity(ion, 1e-6);
    }

    CHECK(nTotalCollisions == 6);
    CHECK(collisionEnergies.size() == 6);
    CHECK(Approx(collisionEnergies[0]) == 0.007162);
    CHECK(Approx(collisionEnergies[1]) == 0.0224838);
}

