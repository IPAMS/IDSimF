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
 test_parallelTree.cpp

 Testing of parallel version of the barnes-hut tree

 ****************************/

#include "Core_vector.hpp"
#include "BTree_parallelTree.hpp"
#include "PSim_util.hpp"
#include <iostream>
#include "test_util.hpp"
#include "catch.hpp"

TEST_CASE( "Test parallel tree semantics / particle management","[Tree]") {
    BTree::ParallelTree testTree(
            Core::Vector(-2.0, -2.0, -2.0),
            Core::Vector(2.0, 2.0, 2.0));

    SECTION ("Test particle insert and delete"){
        //Test particle add increases the number of particles:
        BTree::Particle testIon1(Core::Vector(-0.001,0.0,0.0),1.0);
        BTree::Particle testIon2(Core::Vector(0.0,0.0,0.0),3.0);
        BTree::Particle testIon3(Core::Vector(0.001,0.0,0.0),8.0);
        testTree.insertParticle(testIon1,1);
        testTree.insertParticle(testIon2,2);
        testTree.insertParticle(testIon3,3);

        BTree::ParallelNode* root = testTree.getRoot();
        REQUIRE(testTree.getParticleList()->size() == 3);
        REQUIRE(root->getNumberOfParticles() == 3);

        //test particle removement by particle reference:
        testTree.removeParticle(1);
        REQUIRE(testTree.getParticleList()->size() == 2);
        REQUIRE(root->getNumberOfParticles() == 2);
        REQUIRE( (root->getCharge() - 11*Core::ELEMENTARY_CHARGE) < 1e-100);

        testTree.removeParticle(2);
        REQUIRE( root->getCenterOfCharge() == Core::Vector(0.001,0.0,0.0));
        REQUIRE( (root->getCharge() - 8*Core::ELEMENTARY_CHARGE) < 1e-100);

        testTree.removeParticle(3);

        REQUIRE( root->getCenterOfCharge() == Core::Vector(0.0,0.0,0.0));
        REQUIRE( root->getCharge() < 1e-100);
    }

    SECTION( "Test particle insertion edge cases"){
        BTree::Particle testIon1(Core::Vector(1.0,1.0,1.0),1.0);
        BTree::Particle testIon2(Core::Vector(1.0,1.0,1.0),1.0);

        testTree.insertParticle(testIon1,1);
        REQUIRE_THROWS(testTree.insertParticle(testIon2,2));

        REQUIRE( testTree.getNumberOfParticles() == 1);
        REQUIRE(testTree.getRoot()->getNumberOfParticles() == 1);

        REQUIRE(testTree.getParticle(1) == &testIon1);
    }

    SECTION( "Test external index particle access"){
        //Test particle add increases the number of particles:
        BTree::Particle testIon1(Core::Vector(1.0,1.0,1.0), 1.0);
        BTree::Particle testIon2(Core::Vector(1.1,1.0,1.0), 1.0);
        BTree::Particle testIon3(Core::Vector(1.2,1.0,1.0), 1.0);
        BTree::Particle testIon4(Core::Vector(1.3,1.0,1.0), 1.0);
        BTree::Particle testIon5(Core::Vector(1.4,1.0,1.0), 1.0);

        testTree.insertParticle(testIon1,10);
        testTree.insertParticle(testIon2,20);
        testTree.insertParticle(testIon3,30);
        testTree.insertParticle(testIon4,40);

        REQUIRE( (testTree.getParticle(30) == &testIon3));
        REQUIRE( (testTree.getParticle(40) == &testIon4));

        testTree.insertParticle(testIon5,15);

        REQUIRE( (testTree.getParticle(15) == &testIon5));

        testTree.removeParticle(20);
        testTree.removeParticle(30);

        REQUIRE( (testTree.getParticle(15) == &testIon5));
        REQUIRE( (testTree.getParticle(10) == &testIon1));

        REQUIRE( testTree.getNumberOfParticles() == 3);
    }

    SECTION ( "Test particle update"){
        BTree::ParallelTree testTree_2(
                Core::Vector(0,0,0),
                Core::Vector(10,10,10));

        BTree::Particle testIon1(Core::Vector(1.0,1.0,1.0), 1.0);
        testTree_2.insertParticle(testIon1,10);
        BTree::Particle testIon2(Core::Vector(2.0,1.0,1.0), 1.0);
        testTree_2.insertParticle(testIon2,20);

        int updated = 0;
        testTree_2.updateParticleLocation(10, Core::Vector(2.01, 1.0, 1.0), &updated);
        testTree_2.updateNodes(updated);

        REQUIRE( testIon1.getLocation() == Core::Vector(2.01,1.0,1.0));
        REQUIRE( testTree_2.getNumberOfParticles() == 2);
        REQUIRE( testTree_2.getRoot()->getCenterOfCharge() == Core::Vector(2.005,1.0,1.0));

        testTree_2.updateParticleLocation(10, Core::Vector(2.01001, 1.0, 1.0), &updated);
        testTree_2.updateNodes(updated);
        REQUIRE( testIon1.getLocation() == Core::Vector(2.01001,1.0,1.0));
        REQUIRE( testTree_2.getRoot()->getCenterOfCharge() == Core::Vector(2.005005,1.0,1.0));
    }

    SECTION( "Test tree integrity with large number of random particles"){
        int nions = 10000;
        Core::Vector boxSize(0.002, 0.002, 0.002);
        ParticleSimulation::BoxStartZone startZone(boxSize);
        std::vector<std::unique_ptr<BTree::Particle>> ions= startZone.getRandomParticlesInStartZone(nions, 1);

        for (int i=0; i<nions; i++){
            testTree.insertParticle((*ions[i]),i+1);
        }
        testTree.init();
        //The integrity of the resulting large random tree should be valid:
        REQUIRE_NOTHROW(testTree.getRoot()->testNodeIntegrity(0));
    }
}

TEST_CASE( "Test parallel tree charge distribution calculation","[Tree]"){
    BTree::ParallelTree testTree(
            Core::Vector(-2.0, -2.0, -2.0),
            Core::Vector(2.0, 2.0, 2.0));

    SECTION("Test basic particle insertion, deletion and charge distribution calculation"){
        //Test particle add increases the number of particles:
        BTree::Particle testIon1(Core::Vector(1.0,1.0,1.0), 1.0);
        BTree::Particle testIon2(Core::Vector(1.5,1.0,1.0), 1.0);
        BTree::Particle testIon3(Core::Vector(0.5,1.0,1.0), 1.0);
        testTree.insertParticle(testIon1,1);
        testTree.insertParticle(testIon2,2);
        testTree.insertParticle(testIon3,3);

        int numberOfNodesInTree = testTree.init();
        REQUIRE(testTree.getNumberOfParticles() == 3);
        REQUIRE(testTree.getRoot()->getCharge() == 3.0*Core::ELEMENTARY_CHARGE);

        REQUIRE(testTree.computeEFieldFromTree(testIon1) == Core::Vector(0.0,0.0,0.0));
        REQUIRE(testTree.computeEFieldFromTree(testIon2) != Core::Vector(0.0,0.0,0.0));

        std::list<BTree::Particle*>* particleList= testTree.getParticleList();
        REQUIRE(particleList->size() == 3);
        REQUIRE(particleList->front() == &testIon3);
    }

    SECTION("Test particle insertion and carge calculation with another particle configuration"){
        BTree::Particle testIon1(Core::Vector(-0.001,0.0,0.0), 1.0);
        BTree::Particle testIon2(Core::Vector(0.0,0.0,0.0), 1.0);
        BTree::Particle testIon3(Core::Vector(0.001,0.0,0.0), 1.0);
        testTree.insertParticle(testIon1,1);
        testTree.insertParticle(testIon2,2);
        testTree.insertParticle(testIon3,3);

        int numberOfNodesInTree = testTree.init();
        REQUIRE(testTree.computeEFieldFromTree(testIon1) != Core::Vector(0.0,0.0,0.0));
        REQUIRE(
                vectorApproxCompare(
                        testTree.computeEFieldFromTree(testIon2),
                        Core::Vector(0.0,0.0,0.0))
                        ==  vectorsApproxEqual);
        REQUIRE(testTree.computeEFieldFromTree(testIon3) != Core::Vector(0.0,0.0,0.0));
    }

    SECTION( "Test force calculation with a large number of particles in a cube"){

        //Test particle add increases the number of particles:
        BTree::Particle testIon1(Core::Vector(0.0,0.0,0.0), 1.0);
        BTree::Particle testIon2(Core::Vector(-0.0008,0.0,0.0), 1.0);
        BTree::Particle testIon3(Core::Vector(0.0008,0.0,0.0), 1.0);
        BTree::Particle testIon4(Core::Vector(0.0,-0.0008,0.0), 1.0);
        BTree::Particle testIon5(Core::Vector(0.0,0.0008,0.0), 1.0);

        int nions = 10000;
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

        for (int i=0; i<nions; i++){
            testTree.insertParticle(*ions1[i],1);
            testTree.insertParticle(*ions2[i],2);
            testTree.insertParticle(*ions3[i],3);
            testTree.insertParticle(*ions4[i],4);
            testTree.insertParticle(*ions5[i],5);
            testTree.insertParticle(*ions6[i],6);
        }

        //testTree.computeChargeDistributionRecursive();
        int numberOfNodesInTree = testTree.init();
        REQUIRE(testTree.getNumberOfParticles() == 6 * nions);
        Core::Vector centralForce = testTree.computeEFieldFromTree(testIon1);
        Core::Vector leftForce = testTree.computeEFieldFromTree(testIon2);
        Core::Vector topForce = testTree.computeEFieldFromTree(testIon4);
        REQUIRE(std::abs(centralForce.x()) < 1.0 );
        REQUIRE(std::abs(centralForce.y()) < 1.0 );
        REQUIRE(std::abs(centralForce.z()) < 1.0 );

        REQUIRE(std::abs(leftForce.x())-107 < 1.5 );
        REQUIRE(std::abs(leftForce.y()) < 1 );
        REQUIRE(std::abs(leftForce.z()) < 1 );

        REQUIRE(std::abs(topForce.y())-107 < 1.5 );
        REQUIRE(std::abs(topForce.x()) < 1 );
        REQUIRE(std::abs(topForce.z()) < 1 );
    }
}