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
 test_MultiCollisionModel.cpp

 Testing of Multi Collision Model bundling multiple collision models

 ****************************/

#include "CollisionModel_AbstractCollisionModel.hpp"
#include "CollisionModel_MultiCollisionModel.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "CollisionModel_util.hpp"
#include "BTree_particle.hpp"
#include "Core_randomGenerators.hpp"
#include "catch.hpp"

#include <memory>
#include <vector>


TEST_CASE( "Test multi collision model with multiple hard sphere models ","[CollisionModels][MultiModel]") {
    //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
    Core::globalRandomGeneratorPool = std::make_unique<Core::TestRandomGeneratorPool>();

    double diameterHe = CollisionModel::HardSphereModel::DIAMETER_HE;
    double diameterN2 = CollisionModel::HardSphereModel::DIAMETER_N2;

    int collisionCounts[2] ={0,0};
    auto countFct1 = CollisionModel::util::getCollisionCountFunction(&collisionCounts[0]);
    auto countFct2 = CollisionModel::util::getCollisionCountFunction(&collisionCounts[1]);

    auto hs1 = std::make_unique<CollisionModel::HardSphereModel>(1.0,298,4.0,diameterHe,countFct1);
    auto hs2 = std::make_unique<CollisionModel::HardSphereModel>(0.05,298,28.0,diameterN2,countFct2);

    BTree::Particle ion = BTree::Particle();
    ion.setDiameter(diameterN2);
    ion.setVelocity(Core::Vector(400,0,0));
    ion.setMassAMU(28.0);
    int n = 31;

    std::vector<std::unique_ptr<CollisionModel::AbstractCollisionModel>> models;
    models.emplace_back(std::move(hs1));
    models.emplace_back(std::move(hs2));
    CollisionModel::MultiCollisionModel multiModel(std::move(models));

    for (int i =0; i<n; i++){
        multiModel.modifyVelocity(ion, 2e-5);
    }

    CHECK(collisionCounts[0] == 25);
    CHECK(collisionCounts[1] == 2);
}

