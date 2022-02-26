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
 test_tree.cpp

 Testing of serial version of the barnes-hut tree

 ****************************/

#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include "BTree_tree.hpp"
#include "PSim_boxStartZone.hpp"
#include "PSim_util.hpp"
#include "test_particleStarting.hpp"
#include "catch.hpp"
#include "test_util.hpp"
#include <iostream>

TEST_CASE( "Test serial tree semantics / particle management","[Tree]") {
    BTree::Tree testTree(
            Core::Vector(-2.0, -2.0, -2.0),
            Core::Vector(2.0, 2.0, 2.0));

    SECTION ("Test particle insert and delete"){
        //Test particle add increases the number of particles:
        Core::Particle testIon1(Core::Vector(-0.001,0.0,0.0),1.0);
        Core::Particle testIon2(Core::Vector(0.0,0.0,0.0),3.0);
        Core::Particle testIon3(Core::Vector(0.001,0.0,0.0),8.0);

        testTree.insertParticle(testIon1,1);
        testTree.insertParticle(testIon2,2);
        testTree.insertParticle(testIon3,3);
        testTree.computeChargeDistribution();

        BTree::TreeParticle* checkParticlePtr1 = testTree.getParticle(1);
        BTree::TreeParticle* checkParticlePtr2 = testTree.getParticle(2);
        BTree::TreeParticle* checkParticlePtr3 = testTree.getParticle(3);

        BTree::Node* root = testTree.getRoot();
        CHECK(testTree.getParticleList()->size() == 3);
        CHECK(root->getNumberOfParticles() == 3);

        //test particle removement by particle reference:
        testTree.removeParticle(1);
        CHECK(testTree.getParticleList()->size() == 2);
        CHECK(root->getNumberOfParticles() == 2);
        CHECK( (root->getCharge() - 11*Core::ELEMENTARY_CHARGE) < 1e-100);

        testTree.removeParticle(2);
        CHECK( root->getCenterOfCharge() == Core::Vector(0.001,0.0,0.0));
        CHECK( (root->getCharge() - 8*Core::ELEMENTARY_CHARGE) < 1e-100);
        CHECK( ! root->isParticleInSubtree(checkParticlePtr1,false));
        CHECK( ! root->isParticleInSubtree(checkParticlePtr2,false));
        CHECK( root->isParticleInSubtree(checkParticlePtr3,false));

        testTree.removeParticle(3);

        CHECK( root->getCenterOfCharge() == Core::Vector(0.0,0.0,0.0));
        CHECK( root->getCharge() < 1e-100);
        CHECK( ! root->isParticleInSubtree(checkParticlePtr3,false));
    }

    SECTION( "Test particle insertion edge cases"){
        Core::Particle testIon1(Core::Vector(1.0,1.0,1.0),1.0);
        Core::Particle testIon2(Core::Vector(1.0,1.0,1.0),1.0);

        testTree.insertParticle(testIon1,1);
        CHECK_THROWS(testTree.insertParticle(testIon2,2));

        CHECK( testTree.getNumberOfParticles() == 1);
        CHECK(testTree.getRoot()->getNumberOfParticles() == 1);

        CHECK(testTree.getParticle(1)->wrappedParticle == &testIon1);
    }

    SECTION( "Test external index particle access"){
        //Test particle add increases the number of particles:
        Core::Particle testIon1(Core::Vector(1.0,1.0,1.0), 1.0);
        Core::Particle testIon2(Core::Vector(1.1,1.0,1.0), 1.0);
        Core::Particle testIon3(Core::Vector(1.2,1.0,1.0), 1.0);
        Core::Particle testIon4(Core::Vector(1.3,1.0,1.0), 1.0);
        Core::Particle testIon5(Core::Vector(1.4,1.0,1.0), 1.0);

        testTree.insertParticle(testIon1,10);
        testTree.insertParticle(testIon2,20);
        testTree.insertParticle(testIon3,30);
        testTree.insertParticle(testIon4,40);

        CHECK( (testTree.getParticle(30)->wrappedParticle == &testIon3));
        CHECK( (testTree.getParticle(40)->wrappedParticle == &testIon4));

        testTree.insertParticle(testIon5,15);

        CHECK( (testTree.getParticle(15)->wrappedParticle == &testIon5));

        testTree.removeParticle(20);
        testTree.removeParticle(30);

        CHECK( (testTree.getParticle(15)->wrappedParticle == &testIon5));
        CHECK( (testTree.getParticle(10)->wrappedParticle == &testIon1));

        CHECK( testTree.getNumberOfParticles() == 3);
    }

    SECTION("Test particle update"){
        BTree::Tree testTree_2 = BTree::Tree(
                Core::Vector(0,0,0),
                Core::Vector(10,10,10));

        Core::Particle testIon1(Core::Vector(1.0,1.0,1.0), 1.0);
        testTree_2.insertParticle(testIon1,10);
        Core::Particle testIon2(Core::Vector(2.0,1.0,1.0), 1.0);
        testTree_2.insertParticle(testIon2,20);
        testTree_2.computeChargeDistribution();

        testTree_2.updateParticleLocation(10, Core::Vector(2.01,1.0,1.0));
        CHECK( testIon1.getLocation() == Core::Vector(2.01,1.0,1.0));
        CHECK( testTree_2.getNumberOfParticles() == 2);
        CHECK( testTree_2.getRoot()->getCenterOfCharge() == Core::Vector(2.005,1.0,1.0));

        testTree_2.updateParticleLocation(10,Core::Vector(2.01001,1.0,1.0));
        CHECK( testIon1.getLocation() == Core::Vector(2.01001,1.0,1.0));
        CHECK( testTree_2.getRoot()->getCenterOfCharge() == Core::Vector(2.005005,1.0,1.0));
    }

    SECTION( "Test tree integrity with large number of random particles"){
        std::size_t nions = 10000;
        Core::Vector boxSize(0.002, 0.002, 0.002);
        ParticleSimulation::BoxStartZone startZone(boxSize);
        std::vector<std::unique_ptr<Core::Particle>> ions= startZone.getRandomParticlesInStartZone(nions, 1);

        for (std::size_t i=0; i<nions; i++){
            testTree.insertParticle((*ions[i]), i+1);
        }
        testTree.computeChargeDistribution();

        //The integrity of the resulting large random tree should be valid:
        CHECK_NOTHROW(testTree.getRoot()->testSpatialTreeIntegrity());
        CHECK_NOTHROW(testTree.getRoot()->testNodeIntegrity(0));
        CHECK_NOTHROW(testTree.getRoot()->testNodeParticleIntegrity());
    }
}

TEST_CASE( "Test serial tree charge distribution calculation","[Tree]"){
    BTree::Tree testTree(
            Core::Vector(-2.0, -2.0, -2.0),
            Core::Vector(2.0, 2.0, 2.0));

    SECTION("Test basic particle insertion, deletion and charge distribution calculation"){
        //Test particle add increases the number of particles:
        Core::Particle testIon1(Core::Vector(1.0,1.0,1.0), 1.0);
        Core::Particle testIon2(Core::Vector(1.5,1.0,1.0), 1.0);
        Core::Particle testIon3(Core::Vector(0.5,1.0,1.0), 1.0);
        testTree.insertParticle(testIon1,1);
        testTree.insertParticle(testIon2,2);
        testTree.insertParticle(testIon3,3);

        testTree.computeChargeDistribution();

        CHECK(testTree.getNumberOfParticles() == 3);
        CHECK(isExactDoubleEqual(testTree.getRoot()->getCharge(), 3.0*Core::ELEMENTARY_CHARGE));

        CHECK(testTree.getEFieldFromSpaceCharge(testIon1) == Core::Vector(0.0,0.0,0.0));
        CHECK(testTree.getEFieldFromSpaceCharge(testIon2) != Core::Vector(0.0,0.0,0.0));

        BTree::treeParticlePtrList * particleList= testTree.getParticleList();
        CHECK(particleList->size() == 3);
        CHECK(particleList->front()->wrappedParticle == &testIon3);
    }

    SECTION("Test particle insertion and carge calculation with another particle configuration"){
        Core::Particle testIon1(Core::Vector(-0.001,0.0,0.0), 1.0);
        Core::Particle testIon2(Core::Vector(0.0,0.0,0.0), 1.0);
        Core::Particle testIon3(Core::Vector(0.001,0.0,0.0), 1.0);
        testTree.insertParticle(testIon1,1);
        testTree.insertParticle(testIon2,2);
        testTree.insertParticle(testIon3,3);

        testTree.computeChargeDistribution();
        CHECK(testTree.getEFieldFromSpaceCharge(testIon1) != Core::Vector(0.0,0.0,0.0));
        CHECK(
                vectorApproxCompare(
                        testTree.getEFieldFromSpaceCharge(testIon2),
                        Core::Vector(0.0,0.0,0.0))
                        ==  vectorsApproxEqual);
        CHECK(testTree.getEFieldFromSpaceCharge(testIon3) != Core::Vector(0.0,0.0,0.0));
    }

    SECTION( "Test force calculation with a large number of particles in a cube"){

        //Test particle add increases the number of particles:
        Core::Particle testIon1(Core::Vector(0.0,0.0,0.0), 1.0);
        Core::Particle testIon2(Core::Vector(-0.0008,0.0,0.0), 1.0);
        Core::Particle testIon3(Core::Vector(0.0008,0.0,0.0), 1.0);
        Core::Particle testIon4(Core::Vector(0.0,-0.0008,0.0), 1.0);
        Core::Particle testIon5(Core::Vector(0.0,0.0008,0.0), 1.0);

        std::size_t nions = 10000;
        Core::Vector corner1(-0.00051,-0.0005,-0.0005);
        Core::Vector boxSize1(0.00001,0.001,0.001);
        auto ions1= getRandomIonsInBox(nions,corner1,boxSize1);

        Core::Vector corner2(0.0005,-0.0005,-0.0005);
        Core::Vector boxSize2(0.00001,0.001,0.001);
        auto ions2= getRandomIonsInBox(nions,corner2,boxSize2);

        Core::Vector corner3(-0.0005,-0.00051,-0.0005);
        Core::Vector boxSize3(0.001,0.00001,0.001);
        auto ions3= getRandomIonsInBox(nions,corner3,boxSize3);

        Core::Vector corner4(-0.0005,0.0005,-0.0005);
        Core::Vector boxSize4(0.001,0.00001,0.001);
        auto ions4= getRandomIonsInBox(nions,corner4,boxSize4);

        Core::Vector corner5(-0.0005,-0.0005,-0.00051);
        Core::Vector boxSize5(0.001,0.001,0.00001);
        auto ions5= getRandomIonsInBox(nions,corner5,boxSize5);

        Core::Vector corner6(-0.0005,-0.0005,0.0005);
        Core::Vector boxSize6(0.001,0.001,0.00001);
        auto ions6= getRandomIonsInBox(nions,corner6,boxSize6);

        for (std::size_t i=0; i<nions; i++){
            testTree.insertParticle(*ions1[i],1);
            testTree.insertParticle(*ions2[i],2);
            testTree.insertParticle(*ions3[i],3);
            testTree.insertParticle(*ions4[i],4);
            testTree.insertParticle(*ions5[i],5);
            testTree.insertParticle(*ions6[i],6);
        }

        testTree.computeChargeDistribution();
        CHECK(testTree.getNumberOfParticles() == 6 * nions);
        Core::Vector centralForce = testTree.getEFieldFromSpaceCharge(testIon1);
        Core::Vector leftForce = testTree.getEFieldFromSpaceCharge(testIon2);
        Core::Vector topForce = testTree.getEFieldFromSpaceCharge(testIon4);
        CHECK(std::abs(centralForce.x()) < 1.0 );
        CHECK(std::abs(centralForce.y()) < 1.0 );
        CHECK(std::abs(centralForce.z()) < 1.0 );

        CHECK(std::abs(leftForce.x())-107 < 1.5 );
        CHECK(std::abs(leftForce.y()) < 1 );
        CHECK(std::abs(leftForce.z()) < 1 );

        CHECK(std::abs(topForce.y())-107 < 1.5 );
        CHECK(std::abs(topForce.x()) < 1 );
        CHECK(std::abs(topForce.z()) < 1 );
    }
}
