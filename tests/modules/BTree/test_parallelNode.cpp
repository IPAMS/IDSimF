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
 test_parallelNode.cpp

 Testing of nodes of the parallel version of the Barnes-Hut tree

 ****************************/

#include "BTree_parallelNode.hpp"
#include "BTree_particle.hpp"
#include "Core_constants.hpp"
#include "test_util.hpp"
#include <iostream>
#include <cmath>
#include "catch.hpp"

TEST_CASE( "Test basic parallel node semantics", "[Node]"){

    SECTION( "Parallel nodes should initialize correctly", "[Node]") {
        //Test ion location:
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(1.0,1.0,1.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );

        //Test bare node is root:
        REQUIRE( testNode.isRoot() == true);

        //Test center is calculated correctly:
        REQUIRE( testNode.getCenter() == Core::Vector(1.5,1.5,1.5));

        //Test corners of the block are returned correctly:
        REQUIRE( testNode.getMin() == Core::Vector(1.0,1.0,1.0));
        REQUIRE( testNode.getMax() == Core::Vector(2.0,2.0,2.0));
    }

    SECTION( "Sub octants should have the correct geometry / locations", "[Node]") {
        //Test ion location:
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(-2.0,-2.0,-2.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );

        REQUIRE( testNode.getOctant(Core::Vector(1.0,1.0,1.0)) == BTree::ParallelNode::NET);
        REQUIRE( testNode.getOctant(Core::Vector(-1.0,1.0,1.0)) == BTree::ParallelNode::NWT);
        REQUIRE( testNode.getOctant(Core::Vector(1.0,-1.0,1.0)) == BTree::ParallelNode::SET);
        REQUIRE( testNode.getOctant(Core::Vector(-1.0,-1.0,1.0)) == BTree::ParallelNode::SWT);
        REQUIRE( testNode.getOctant(Core::Vector(1.0,1.0,-1.0)) == BTree::ParallelNode::NEB);
        REQUIRE( testNode.getOctant(Core::Vector(-1.0,1.0,-1.0)) == BTree::ParallelNode::NWB);
        REQUIRE( testNode.getOctant(Core::Vector(1.0,-1.0,-1.0)) == BTree::ParallelNode::SEB);
        REQUIRE( testNode.getOctant(Core::Vector(-1.0,-1.0,-1.0)) == BTree::ParallelNode::SWB);
    }

    SECTION( "All sub nodes / octants of a parallel node should be creatable and correct") {
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(1.0,1.0,1.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );

        //Explicitly test the generation of ALL subnodes:
        BTree::ParallelNode* subNode = testNode.createOctNode(BTree::ParallelNode::SWB);
        REQUIRE( subNode->getCenter() == Core::Vector(1.25,1.25,1.25));

        subNode = testNode.createOctNode(BTree::ParallelNode::SWT);
        REQUIRE( subNode->getCenter() == Core::Vector(1.25,1.25,1.75));

        subNode = testNode.createOctNode(BTree::ParallelNode::SEB);
        REQUIRE( subNode->getCenter() == Core::Vector(1.75,1.25,1.25));

        subNode = testNode.createOctNode(BTree::ParallelNode::SET);
        REQUIRE( subNode->getCenter() == Core::Vector(1.75,1.25,1.75));

        subNode = testNode.createOctNode(BTree::ParallelNode::NWB);
        REQUIRE( subNode->getCenter() == Core::Vector(1.25,1.75,1.25));

        subNode = testNode.createOctNode(BTree::ParallelNode::NWT);
        REQUIRE( subNode->getCenter() == Core::Vector(1.25,1.75,1.75));

        subNode = testNode.createOctNode(BTree::ParallelNode::NET);
        REQUIRE( subNode->getCenter() == Core::Vector(1.75,1.75,1.75));

        subNode = testNode.createOctNode(BTree::ParallelNode::NEB);
        REQUIRE( subNode->getCenter() == Core::Vector(1.75,1.75,1.25));
    }
}

TEST_CASE( "Test particle management in parallel node", "[Node]") {

    SECTION( "Particles should insert correctly into parallel node", "[Node]") {
        //Test ion location:
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(1.0,1.0,1.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );


        //Test bare node has no particles:
        REQUIRE( testNode.getNumberOfParticles() == 0);
        REQUIRE( testNode.getParticle() == nullptr);

        //Test particle add increases the number of particles:
        BTree::Particle testIon1 = BTree::Particle(Core::Vector(1.2,1.2,1.2),2.0);
        BTree::Particle testIon2 = BTree::Particle(Core::Vector(1.6,1.2,1.2),2.0);
        BTree::Particle testIon3 = BTree::Particle(Core::Vector(1.65,1.2,1.2),2.0);
        BTree::Particle testIon4 = BTree::Particle(Core::Vector(1.66,1.2,1.2),2.0);
        BTree::Particle testIon5 = BTree::Particle(Core::Vector(1.665,1.2,1.2),2.0);

        testNode.insertParticle(&testIon1);
        REQUIRE( testNode.getNumberOfParticles() == 1);

        //Test particle is added as a reference:
        REQUIRE( testNode.getParticle() == &testIon1);


        testNode.insertParticle(&testIon2);
        REQUIRE( testNode.getNumberOfParticles() == 2);

        testNode.insertParticle(&testIon3);
        testNode.insertParticle(&testIon4);
        testNode.insertParticle(&testIon5);
        REQUIRE( testNode.getNumberOfParticles() == 5);

        //Test if particles host node is set correctly and the host node has correct characteristics:
        BTree::AbstractNode* hostNode5 = testIon5.getHostNode();
        Core::Vector pos = testIon5.getLocation();
        Core::Vector min = hostNode5->getMin();
        Core::Vector max = hostNode5->getMax();

        REQUIRE(hostNode5->getNumberOfParticles() == 1);
        REQUIRE((min.x() < pos.x() && pos.x() < max.x()));
        REQUIRE((min.y() < pos.y() && pos.y() < max.y()));
        REQUIRE((min.z() < pos.z() && pos.z() < max.z()));
    }

    SECTION("Edge cases of particle insertion into parallel node should be handled correctly"){

        //test particle insertion with same particle position:
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(0.0,0.0,0.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );

        BTree::Particle testIon1 = BTree::Particle(Core::Vector(1.0,1.0,1.0),2.0);
        BTree::Particle testIon2 = BTree::Particle(Core::Vector(1.0,1.0,1.0),2.0);

        testNode.insertParticle(&testIon1);
        REQUIRE_THROWS(testNode.insertParticle(&testIon2));

    }

    SECTION( "Particles remove from parallel node correctly"){
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(0.0,0.0,0.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );


        //Test bare node has no particles:
        REQUIRE( testNode.getNumberOfParticles() == 0);
        REQUIRE( testNode.getParticle() == nullptr);

        //Test particle add increases the number of particles:
        BTree::Particle testIon1 = BTree::Particle(Core::Vector(1.2,1.2,1.2),2.0);
        BTree::Particle testIon2 = BTree::Particle(Core::Vector(1.6,1.2,1.2),2.0);
        BTree::Particle testIon3 = BTree::Particle(Core::Vector(1.65,1.2,1.2),2.0);
        BTree::Particle testIon4 = BTree::Particle(Core::Vector(1.66,1.2,1.2),2.0);
        BTree::Particle testIon5 = BTree::Particle(Core::Vector(1.5,1.2,1.2),6.0);
        BTree::Particle testIon6 = BTree::Particle(Core::Vector(1.6,1.2,1.2),2.0);
        testNode.insertParticle(&testIon5);
        testNode.insertParticle(&testIon6);

        BTree::AbstractNode* leafNode6 = testIon6.getHostNode();
        leafNode6->removeMyselfFromTree();
        REQUIRE(testNode.getNumberOfParticles() == 1);
        REQUIRE(testNode.getCharge() == 6.0*Core::ELEMENTARY_CHARGE);

        BTree::AbstractNode* leafNode5 = testIon5.getHostNode();
        leafNode5->removeMyselfFromTree();
        REQUIRE(testNode.getNumberOfParticles() == 0);
        REQUIRE( testNode.getCharge() == 0.0*Core::ELEMENTARY_CHARGE);
        testNode.insertParticle(&testIon1);
        testNode.insertParticle(&testIon2);
        testNode.insertParticle(&testIon3);
        testNode.insertParticle(&testIon4);
        testNode.insertParticle(&testIon5);
        //testNode.insertParticle(&testIon6);

        //testNode.computeChargeDistributionRecursive();
        REQUIRE(testNode.getNumberOfParticles() == 5);
        REQUIRE( (testNode.getCharge() - 16.0*Core::ELEMENTARY_CHARGE) < 1e-100);

        BTree::AbstractNode* leafNode1 = testIon1.getHostNode();
        leafNode1->removeMyselfFromTree();
        REQUIRE(testNode.getNumberOfParticles() == 4);
        REQUIRE( (testNode.getCharge() - 14.0*Core::ELEMENTARY_CHARGE) < 1e-100);

        BTree::AbstractNode* leafNode2 = testIon2.getHostNode();
        leafNode2->removeMyselfFromTree();
        REQUIRE(testNode.getNumberOfParticles() == 3);
        REQUIRE( (testNode.getCharge() - 12.0*Core::ELEMENTARY_CHARGE) < 1e-100);
    }

    SECTION("Particles should be removable from node, center of charge should update"){
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(0.0,0.0,0.0),
                Core::Vector(4.0,4.0,4.0),
                nullptr
        );

        //Test bare node has no particles:
        REQUIRE( testNode.getNumberOfParticles() == 0);
        REQUIRE( testNode.getParticle() == nullptr);

        //Test particle add increases the number of particles:
        BTree::Particle testIon1 = BTree::Particle(Core::Vector(1.0,1.0,1.0),1.0);
        BTree::Particle testIon2 = BTree::Particle(Core::Vector(2.0,1.0,1.0),1.0);
        BTree::Particle testIon3 = BTree::Particle(Core::Vector(3.0,1.0,1.0),1.0);
        testNode.insertParticle(&testIon1);
        testNode.insertParticle(&testIon2);
        testNode.insertParticle(&testIon3);

        //compute the charge distribution in the tree with testNode as root recursively:
        testNode.computeChargeDistributionRecursive();

        REQUIRE(testNode.getCenterOfCharge() == Core::Vector(2.0,1.0,1.0));
        BTree::AbstractNode* leafNode3 = testIon3.getHostNode();
        leafNode3->removeMyselfFromTree();
        REQUIRE( (testNode.getCenterOfCharge() - Core::Vector(1.5,1.0,1.0)).magnitude() < 1e-10 );
        REQUIRE(testNode.getCharge() == 2.0*Core::ELEMENTARY_CHARGE);

        BTree::AbstractNode* leafNode1 = testIon1.getHostNode();
        leafNode1->removeMyselfFromTree();
        REQUIRE( (testNode.getCenterOfCharge() - Core::Vector(2.0,1.0,1.0)).magnitude() < 1e-10 );
        REQUIRE(testNode.getCharge() == 1.0*Core::ELEMENTARY_CHARGE);
    }

    SECTION( "Particles should be mass insertable and removeable from parallel node (memory leak test)"){
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(0.0,0.0,0.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );

        int nNodesOriginal= testNode.getNumberOfNodes();
        const int nIons = 100;
        BTree::Particle testIons[nIons];
        double xPos;
        for (int i=0; i<nIons;i++){
            xPos = 1.99 / nIons *i;
            testIons[i] = BTree::Particle(Core::Vector(xPos,1.2,1.2),2.0);
            testNode.insertParticle(&testIons[i]);
        }
        REQUIRE(testNode.getNumberOfParticles() == nIons);
        REQUIRE( (testNode.getCharge() - nIons*Core::ELEMENTARY_CHARGE) < 1e-100);

        int maxRek = testNode.maximumRecursionDepth();
        REQUIRE(maxRek == 8);

        BTree::AbstractNode* hostNode;
        for (int i=0; i<nIons;i++){
            hostNode = testIons[i].getHostNode();
            hostNode->removeMyselfFromTree();
            if (hostNode != &testNode){
                delete (hostNode);
            }
        }
        REQUIRE(testNode.getNumberOfParticles() == 0);
        REQUIRE( (testNode.getCharge() - 0.0) < 1e-100);
        REQUIRE( testNode.getNumberOfNodes() == nNodesOriginal);
    }
}

TEST_CASE( "Test electric field calculation in parallel node","[Node]") {

    SECTION("Test point to point force calculation") {
        //Test static charge calculation methods:
        Core::Vector a = Core::Vector(1.0, 1.0, 1.0);
        Core::Vector b = Core::Vector(2.0, 1.0, 1.0);

        Core::Vector c = BTree::ParallelNode::calculateElectricField(a, b, 1.0);
        Core::Vector d = Core::Vector(-(1.0/(4*M_PI*8.854e-12)), 0.0, 0.0);

        REQUIRE(
                vectorApproxCompare(
                        BTree::ParallelNode::calculateElectricField(a, b, 1.0),
                        BTree::ParallelNode::calculateElectricField(b, a, 1.0)*(-1))
                        ==vectorsApproxEqual);

        REQUIRE(vectorApproxCompare(c, d)==vectorsApproxEqual);
    }
}