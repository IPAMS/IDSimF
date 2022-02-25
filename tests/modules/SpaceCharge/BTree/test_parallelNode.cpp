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
#include "Core_particle.hpp"
#include "Core_constants.hpp"
#include "catch.hpp"
#include "test_util.hpp"
#include <iostream>
#include <cmath>

TEST_CASE( "Test basic parallel node semantics", "[Node]"){

    SECTION( "Parallel nodes should initialize correctly", "[Node]") {
        //Test ion location:
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(1.0,1.0,1.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );

        //Test bare node is root:
        CHECK( testNode.isRoot() == true);

        //Test center is calculated correctly:
        CHECK( testNode.getCenter() == Core::Vector(1.5,1.5,1.5));

        //Test corners of the block are returned correctly:
        CHECK( testNode.getMin() == Core::Vector(1.0,1.0,1.0));
        CHECK( testNode.getMax() == Core::Vector(2.0,2.0,2.0));
    }

    SECTION( "Sub octants should have the correct geometry / locations", "[Node]") {
        //Test ion location:
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(-2.0,-2.0,-2.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );

        CHECK( testNode.getOctant(Core::Vector(1.0,1.0,1.0)) == BTree::ParallelNode::NET);
        CHECK( testNode.getOctant(Core::Vector(-1.0,1.0,1.0)) == BTree::ParallelNode::NWT);
        CHECK( testNode.getOctant(Core::Vector(1.0,-1.0,1.0)) == BTree::ParallelNode::SET);
        CHECK( testNode.getOctant(Core::Vector(-1.0,-1.0,1.0)) == BTree::ParallelNode::SWT);
        CHECK( testNode.getOctant(Core::Vector(1.0,1.0,-1.0)) == BTree::ParallelNode::NEB);
        CHECK( testNode.getOctant(Core::Vector(-1.0,1.0,-1.0)) == BTree::ParallelNode::NWB);
        CHECK( testNode.getOctant(Core::Vector(1.0,-1.0,-1.0)) == BTree::ParallelNode::SEB);
        CHECK( testNode.getOctant(Core::Vector(-1.0,-1.0,-1.0)) == BTree::ParallelNode::SWB);
    }

    SECTION( "All sub nodes / octants of a parallel node should be creatable and correct") {
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(1.0,1.0,1.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );

        //Explicitly test the generation of ALL subnodes:
        BTree::ParallelNode* subNode = testNode.createOctNode(BTree::ParallelNode::SWB);
        CHECK( subNode->getCenter() == Core::Vector(1.25,1.25,1.25));

        subNode = testNode.createOctNode(BTree::ParallelNode::SWT);
        CHECK( subNode->getCenter() == Core::Vector(1.25,1.25,1.75));

        subNode = testNode.createOctNode(BTree::ParallelNode::SEB);
        CHECK( subNode->getCenter() == Core::Vector(1.75,1.25,1.25));

        subNode = testNode.createOctNode(BTree::ParallelNode::SET);
        CHECK( subNode->getCenter() == Core::Vector(1.75,1.25,1.75));

        subNode = testNode.createOctNode(BTree::ParallelNode::NWB);
        CHECK( subNode->getCenter() == Core::Vector(1.25,1.75,1.25));

        subNode = testNode.createOctNode(BTree::ParallelNode::NWT);
        CHECK( subNode->getCenter() == Core::Vector(1.25,1.75,1.75));

        subNode = testNode.createOctNode(BTree::ParallelNode::NET);
        CHECK( subNode->getCenter() == Core::Vector(1.75,1.75,1.75));

        subNode = testNode.createOctNode(BTree::ParallelNode::NEB);
        CHECK( subNode->getCenter() == Core::Vector(1.75,1.75,1.25));
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
        CHECK( testNode.getNumberOfParticles() == 0);
        CHECK( testNode.getParticle() == nullptr);

        //Test particle add increases the number of particles:
        Core::Particle baseIon1(Core::Vector(1.2,1.2,1.2),2.0);
        Core::Particle baseIon2(Core::Vector(1.6,1.2,1.2),2.0);
        Core::Particle baseIon3(Core::Vector(1.65,1.2,1.2),2.0);
        Core::Particle baseIon4(Core::Vector(1.66,1.2,1.2),2.0);
        Core::Particle baseIon5(Core::Vector(1.665,1.2,1.2),2.0);

        BTree::TreeParticle testIon1(&baseIon1);
        BTree::TreeParticle testIon2(&baseIon2);
        BTree::TreeParticle testIon3(&baseIon3);
        BTree::TreeParticle testIon4(&baseIon4);
        BTree::TreeParticle testIon5(&baseIon5);

        testNode.insertParticle(&testIon1);
        CHECK( testNode.getNumberOfParticles() == 1);

        //Test particle is added as a reference:
        CHECK( testNode.getParticle() == &testIon1);


        testNode.insertParticle(&testIon2);
        CHECK( testNode.getNumberOfParticles() == 2);

        testNode.insertParticle(&testIon3);
        testNode.insertParticle(&testIon4);
        testNode.insertParticle(&testIon5);
        CHECK( testNode.getNumberOfParticles() == 5);

        //Test if particles host node is set correctly and the host node has correct characteristics:
        BTree::AbstractNode* hostNode5 = testIon5.getHostNode();
        Core::Vector pos = testIon5.get()->getLocation();
        Core::Vector min = hostNode5->getMin();
        Core::Vector max = hostNode5->getMax();

        CHECK(hostNode5->getNumberOfParticles() == 1);
        CHECK((min.x() < pos.x() && pos.x() < max.x()));
        CHECK((min.y() < pos.y() && pos.y() < max.y()));
        CHECK((min.z() < pos.z() && pos.z() < max.z()));
    }

    SECTION("Edge cases of particle insertion into parallel node should be handled correctly"){

        //test particle insertion with same particle position:
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(0.0,0.0,0.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );

        Core::Particle baseIon1 = Core::Particle(Core::Vector(1.0,1.0,1.0),2.0);
        Core::Particle baseIon2 = Core::Particle(Core::Vector(1.0,1.0,1.0),2.0);

        BTree::TreeParticle testIon1(&baseIon1);
        BTree::TreeParticle testIon2(&baseIon2);

        testNode.insertParticle(&testIon1);
        CHECK_THROWS(testNode.insertParticle(&testIon2));

    }

    SECTION( "Particles remove from parallel node correctly"){
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(0.0,0.0,0.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );


        //Test bare node has no particles:
        CHECK( testNode.getNumberOfParticles() == 0);
        CHECK( testNode.getParticle() == nullptr);

        //Test particle add increases the number of particles:
        Core::Particle baseIon1 = Core::Particle(Core::Vector(1.2,1.2,1.2),2.0);
        Core::Particle baseIon2 = Core::Particle(Core::Vector(1.6,1.2,1.2),2.0);
        Core::Particle baseIon3 = Core::Particle(Core::Vector(1.65,1.2,1.2),2.0);
        Core::Particle baseIon4 = Core::Particle(Core::Vector(1.66,1.2,1.2),2.0);
        Core::Particle baseIon5 = Core::Particle(Core::Vector(1.5,1.2,1.2),6.0);
        Core::Particle baseIon6 = Core::Particle(Core::Vector(1.6,1.2,1.2),2.0);

        BTree::TreeParticle testIon1(&baseIon1);
        BTree::TreeParticle testIon2(&baseIon2);
        BTree::TreeParticle testIon3(&baseIon3);
        BTree::TreeParticle testIon4(&baseIon4);
        BTree::TreeParticle testIon5(&baseIon5);
        BTree::TreeParticle testIon6(&baseIon6);

        testNode.insertParticle(&testIon5);
        testNode.insertParticle(&testIon6);

        BTree::AbstractNode* leafNode6 = testIon6.getHostNode();
        leafNode6->removeMyselfFromTree();
        CHECK(testNode.getNumberOfParticles() == 1);
        CHECK(isExactDoubleEqual(testNode.getCharge(), 6.0*Core::ELEMENTARY_CHARGE));

        BTree::AbstractNode* leafNode5 = testIon5.getHostNode();
        leafNode5->removeMyselfFromTree();
        CHECK(testNode.getNumberOfParticles() == 0);
        CHECK(isExactDoubleEqual(testNode.getCharge(), 0.0*Core::ELEMENTARY_CHARGE));
        testNode.insertParticle(&testIon1);
        testNode.insertParticle(&testIon2);
        testNode.insertParticle(&testIon3);
        testNode.insertParticle(&testIon4);
        testNode.insertParticle(&testIon5);
        //testNode.insertParticle(&testIon6);

        //testNode.computeChargeDistributionRecursive();
        CHECK(testNode.getNumberOfParticles() == 5);
        CHECK( (testNode.getCharge() - 16.0*Core::ELEMENTARY_CHARGE) < 1e-100);

        BTree::AbstractNode* leafNode1 = testIon1.getHostNode();
        leafNode1->removeMyselfFromTree();
        CHECK(testNode.getNumberOfParticles() == 4);
        CHECK( (testNode.getCharge() - 14.0*Core::ELEMENTARY_CHARGE) < 1e-100);

        BTree::AbstractNode* leafNode2 = testIon2.getHostNode();
        leafNode2->removeMyselfFromTree();
        CHECK(testNode.getNumberOfParticles() == 3);
        CHECK( (testNode.getCharge() - 12.0*Core::ELEMENTARY_CHARGE) < 1e-100);
    }

    SECTION("Particles should be removable from node, center of charge should update"){
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(0.0,0.0,0.0),
                Core::Vector(4.0,4.0,4.0),
                nullptr
        );

        //Test bare node has no particles:
        CHECK( testNode.getNumberOfParticles() == 0);
        CHECK( testNode.getParticle() == nullptr);

        //Test particle add increases the number of particles:
        Core::Particle baseIon1 = Core::Particle(Core::Vector(1.0,1.0,1.0),1.0);
        Core::Particle baseIon2 = Core::Particle(Core::Vector(2.0,1.0,1.0),1.0);
        Core::Particle baseIon3 = Core::Particle(Core::Vector(3.0,1.0,1.0),1.0);

        BTree::TreeParticle testIon1(&baseIon1);
        BTree::TreeParticle testIon2(&baseIon2);
        BTree::TreeParticle testIon3(&baseIon3);

        testNode.insertParticle(&testIon1);
        testNode.insertParticle(&testIon2);
        testNode.insertParticle(&testIon3);

        //compute the charge distribution in the tree with testNode as root recursively:
        testNode.computeChargeDistributionRecursive();

        CHECK(testNode.getCenterOfCharge() == Core::Vector(2.0,1.0,1.0));
        BTree::AbstractNode* leafNode3 = testIon3.getHostNode();
        leafNode3->removeMyselfFromTree();
        CHECK( (testNode.getCenterOfCharge() - Core::Vector(1.5,1.0,1.0)).magnitude() < 1e-10 );
        CHECK(isExactDoubleEqual(testNode.getCharge(), 2.0*Core::ELEMENTARY_CHARGE));

        BTree::AbstractNode* leafNode1 = testIon1.getHostNode();
        leafNode1->removeMyselfFromTree();
        CHECK( (testNode.getCenterOfCharge() - Core::Vector(2.0,1.0,1.0)).magnitude() < 1e-10 );
        CHECK(isExactDoubleEqual(testNode.getCharge(), 1.0*Core::ELEMENTARY_CHARGE));
    }

    SECTION( "Particles should be mass insertable and removeable from parallel node (memory leak test)"){
        BTree::ParallelNode testNode = BTree::ParallelNode(
                Core::Vector(0.0,0.0,0.0),
                Core::Vector(2.0,2.0,2.0),
                nullptr
        );

        int nNodesOriginal= testNode.getNumberOfNodes();
        const int nIons = 100;
        Core::Particle baseIons[nIons];
        std::vector<std::unique_ptr<BTree::TreeParticle>> treeParticles;
        double xPos;
        for (size_t i=0; i<nIons; i++){
            xPos = 1.99 / nIons *i;
            baseIons[i] = Core::Particle(Core::Vector(xPos,1.2,1.2),2.0);
            treeParticles.push_back(std::make_unique<BTree::TreeParticle>(&baseIons[i]));
            testNode.insertParticle(treeParticles[i].get());
        }
        CHECK(testNode.getNumberOfParticles() == nIons);
        CHECK( (testNode.getCharge() - nIons*Core::ELEMENTARY_CHARGE) < 1e-100);

        std::size_t maxRek = testNode.maximumRecursionDepth();
        CHECK(maxRek == 8);

        BTree::AbstractNode* hostNode;
        for (size_t i=0; i<nIons;i++){
            hostNode = treeParticles[i]->getHostNode();
            hostNode->removeMyselfFromTree();
            if (hostNode != &testNode){
                delete (hostNode);
            }
        }
        CHECK(testNode.getNumberOfParticles() == 0);
        CHECK( (testNode.getCharge() - 0.0) < 1e-100);
        CHECK( testNode.getNumberOfNodes() == nNodesOriginal);
    }
}

TEST_CASE( "Test electric field calculation in parallel node","[Node]") {

    SECTION("Test point to point force calculation") {
        //Test static charge calculation methods:
        Core::Vector a = Core::Vector(1.0, 1.0, 1.0);
        Core::Vector b = Core::Vector(2.0, 1.0, 1.0);

        Core::Vector c = BTree::ParallelNode::calculateElectricField(a, b, 1.0);
        Core::Vector d = Core::Vector(-(1.0/(4*M_PI*8.854e-12)), 0.0, 0.0);

        CHECK(
                vectorApproxCompare(
                        BTree::ParallelNode::calculateElectricField(a, b, 1.0),
                        BTree::ParallelNode::calculateElectricField(b, a, 1.0)*(-1))
                        ==vectorsApproxEqual);

        CHECK(vectorApproxCompare(c, d)==vectorsApproxEqual);
    }
}