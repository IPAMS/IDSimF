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
 test_node.cpp

 Testing of nodes of the serial version of the Barnes-Hut tree

 ****************************/
#include "Core_constants.hpp"
#include "BTree_node.hpp"
#include "BTree_tree.hpp"
#include "Core_particle.hpp"
#include "PSim_util.hpp"
#include "PSim_boxStartZone.hpp"
#include "catch.hpp"
#include "test_util.hpp"
#include <iostream>
#include <cmath>

TEST_CASE( "Test basic serial node semantics", "[Node]"){
    //Test ion location:
    BTree::Node testNode(
            Core::Vector(1.0,1.0,1.0),
            Core::Vector(2.0,2.0,2.0),
            nullptr
    );

    SECTION("Test node initialization"){
        BTree::Node** octs = testNode.getOctants();
        CHECK(octs[7] == nullptr);

        //Test bare node is root:
        CHECK(testNode.isRoot() == true);

        //Test center is calculated correctly:
        CHECK(testNode.getCenter() == Core::Vector(1.5,1.5,1.5));

        //Test corners of the block are returned correctly:
        CHECK(testNode.getMin() == Core::Vector(1.0,1.0,1.0));
        CHECK(testNode.getMax() == Core::Vector(2.0,2.0,2.0));
    }

    SECTION("Test node octant / sub nodes semantics") {
        SECTION("Test node octant determination") {
            BTree::Node testNode2(
                    Core::Vector(-2.0, -2.0, -2.0),
                    Core::Vector(2.0, 2.0, 2.0),
                    nullptr);

            CHECK(testNode2.getOctant(Core::Vector(1.0, 1.0, 1.0))==BTree::Node::NET);
            CHECK(testNode2.getOctant(Core::Vector(-1.0, 1.0, 1.0))==BTree::Node::NWT);
            CHECK(testNode2.getOctant(Core::Vector(1.0, -1.0, 1.0))==BTree::Node::SET);
            CHECK(testNode2.getOctant(Core::Vector(-1.0, -1.0, 1.0))==BTree::Node::SWT);
            CHECK(testNode2.getOctant(Core::Vector(1.0, 1.0, -1.0))==BTree::Node::NEB);
            CHECK(testNode2.getOctant(Core::Vector(-1.0, 1.0, -1.0))==BTree::Node::NWB);
            CHECK(testNode2.getOctant(Core::Vector(1.0, -1.0, -1.0))==BTree::Node::SEB);
            CHECK(testNode2.getOctant(Core::Vector(-1.0, -1.0, -1.0))==BTree::Node::SWB);
        }

        SECTION("Test node sub node / octant generation") {
            BTree::Node testNode2 = BTree::Node(
                    Core::Vector(1.0, 1.0, 1.0),
                    Core::Vector(2.0, 2.0, 2.0),
                    nullptr
            );

            //Explicitly test the generation of ALL subnodes:
            BTree::Node* subNode = testNode2.createOctNode(BTree::Node::SWB);
            CHECK(subNode->getCenter()==Core::Vector(1.25, 1.25, 1.25));

            subNode = testNode.createOctNode(BTree::Node::SWT);
            CHECK(subNode->getCenter()==Core::Vector(1.25, 1.25, 1.75));

            subNode = testNode.createOctNode(BTree::Node::SEB);
            CHECK(subNode->getCenter()==Core::Vector(1.75, 1.25, 1.25));

            subNode = testNode.createOctNode(BTree::Node::SET);
            CHECK(subNode->getCenter()==Core::Vector(1.75, 1.25, 1.75));

            subNode = testNode.createOctNode(BTree::Node::NWB);
            CHECK(subNode->getCenter()==Core::Vector(1.25, 1.75, 1.25));

            subNode = testNode.createOctNode(BTree::Node::NWT);
            CHECK(subNode->getCenter()==Core::Vector(1.25, 1.75, 1.75));

            subNode = testNode.createOctNode(BTree::Node::NET);
            CHECK(subNode->getCenter()==Core::Vector(1.75, 1.75, 1.75));

            subNode = testNode.createOctNode(BTree::Node::NEB);
            CHECK(subNode->getCenter()==Core::Vector(1.75, 1.75, 1.25));
        }
    }
}

TEST_CASE( "Test particle insertion and remove in serial node", "[Node]") {
    //Test ion location:
    BTree::Node testNode = BTree::Node(
              Core::Vector(0.0,0.0,0.0),
              Core::Vector(2.0,2.0,2.0),
              nullptr
        );
    
    SECTION("Test basic particle insertion") {
        //Test bare node has no particles:
        CHECK( testNode.getNumberOfParticles() == 0);
        CHECK( testNode.getParticle() == nullptr);

        //Test particle add increases the number of particles:
        Core::Particle baseIon1 = Core::Particle(Core::Vector(1.2,1.2,1.2),2.0);
        Core::Particle baseIon2 = Core::Particle(Core::Vector(1.6,1.2,1.2),2.0);
        Core::Particle baseIon3 = Core::Particle(Core::Vector(1.65,1.2,1.2),2.0);
        Core::Particle baseIon4 = Core::Particle(Core::Vector(1.66,1.2,1.2),2.0);
        Core::Particle baseIon5 = Core::Particle(Core::Vector(1.665,1.2,1.2),2.0);

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
        Core::Vector pos = testIon5.wrappedParticle->getLocation();
        Core::Vector min = hostNode5->getMin();
        Core::Vector max = hostNode5->getMax();

        CHECK(hostNode5->getNumberOfParticles() == 1);
        CHECK((min.x() < pos.x() && pos.x() < max.x()));
        CHECK((min.y() < pos.y() && pos.y() < max.y()));
        CHECK((min.z() < pos.z() && pos.z() < max.z()));
    }

    SECTION("Test particle insertion with special edge cases"){

        //test particle insertion with same particle position:
        //(currently the section is silently ignored)
        Core::Particle baseIon1 = Core::Particle(Core::Vector(1.0,1.0,1.0),2.0);
        Core::Particle baseIon2 = Core::Particle(Core::Vector(1.0,1.0,1.0),2.0);

        BTree::TreeParticle testIon1(&baseIon1);
        BTree::TreeParticle testIon2(&baseIon2);

        testNode.insertParticle(&testIon1);
        CHECK_THROWS(testNode.insertParticle(&testIon2));

        testNode.computeChargeDistributionRecursive();
        testNode.printTree(0);
        CHECK(testNode.getNumberOfParticles() == 1);
    }

    SECTION( "Test particle remove from Node"){

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

        testNode.computeChargeDistributionRecursive();
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

        testNode.computeChargeDistributionRecursive();
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

    SECTION("Test particle remove from Node with center of charge update"){
        BTree::Node testNode2 (
                Core::Vector(0.0,0.0,0.0),
                Core::Vector(4.0,4.0,4.0),
                nullptr
        );


        //Test bare node has no particles:
        CHECK( testNode2.getNumberOfParticles() == 0);
        CHECK( testNode2.getParticle() == nullptr);

        //Test particle add increases the number of particles:
        Core::Particle baseIon1(Core::Vector(1.0,1.0,1.0), 1.0);
        Core::Particle baseIon2(Core::Vector(2.0,1.0,1.0), 1.0);
        Core::Particle baseIon3(Core::Vector(3.0,1.0,1.0), 1.0);

        BTree::TreeParticle testIon1(&baseIon1);
        BTree::TreeParticle testIon2(&baseIon2);
        BTree::TreeParticle testIon3(&baseIon3);

        testNode2.insertParticle(&testIon1);
        testNode2.insertParticle(&testIon2);
        testNode2.insertParticle(&testIon3);
        testNode2.computeChargeDistributionRecursive();

        CHECK(testNode2.getCenterOfCharge() == Core::Vector(2.0,1.0,1.0));
        BTree::AbstractNode* leafNode3 = testIon3.getHostNode();
        leafNode3->removeMyselfFromTree();
        CHECK( (testNode2.getCenterOfCharge() - Core::Vector(1.5,1.0,1.0)).magnitude() < 1e-10 );
        CHECK(isExactDoubleEqual(testNode2.getCharge(), 2.0*Core::ELEMENTARY_CHARGE));

        BTree::AbstractNode* leafNode1 = testIon1.getHostNode();
        leafNode1->removeMyselfFromTree();
        CHECK( (testNode2.getCenterOfCharge() - Core::Vector(2.0,1.0,1.0)).magnitude() < 1e-10 );
        CHECK(isExactDoubleEqual(testNode2.getCharge(), 1.0*Core::ELEMENTARY_CHARGE));
    }

    SECTION("Test particle mass insertion and remove from node (memory leak test)"){

        int nNodesOriginal= testNode.getNumberOfNodes();
        const int nIons = 100;
        Core::Particle testIons[nIons];
        std::vector<std::unique_ptr<BTree::TreeParticle>> treeParticles;
        double xPos;
        for (int i=0; i<nIons;i++){
            xPos = 1.99 / nIons *i;
            testIons[i] = Core::Particle(Core::Vector(xPos,1.2,1.2),2.0);
            treeParticles.push_back(std::make_unique<BTree::TreeParticle>(&testIons[i]));
            testNode.insertParticle(treeParticles.back().get());
        }
        CHECK(testNode.getNumberOfParticles() == nIons);
        CHECK( (testNode.getCharge() - nIons*Core::ELEMENTARY_CHARGE) < 1e-100);


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

    SECTION( "Test particle remove transformation validity"){
        BTree::Node testNode2 = BTree::Node(
                Core::Vector(0.0,0.0,0.0),
                Core::Vector(20.0,20.0,20.0),
                nullptr
        );

        Core::Particle baseIon1(Core::Vector(9.2,0,0), 2.0);
        Core::Particle baseIon2(Core::Vector(9.201,0,0), 2.0);
        Core::Particle baseIon3(Core::Vector(9.3,0,0), 2.0);
        Core::Particle baseIon4(Core::Vector(9.301,0,0), 2.0);

        BTree::TreeParticle testIon1(&baseIon1);
        BTree::TreeParticle testIon2(&baseIon2);
        BTree::TreeParticle testIon3(&baseIon3);
        BTree::TreeParticle testIon4(&baseIon4);

        testNode2.insertParticle(&testIon1);
        testNode2.insertParticle(&testIon2);
        testNode2.insertParticle(&testIon3);
        testNode2.insertParticle(&testIon4);

        //testNode2.printTree(1);
        //std::cout << "------------------------------------------"<<std::endl;

        testIon2.getHostNode()->removeMyselfFromTree();
        testNode2.testNodeIntegrity(0);
        testNode2.testNodeParticleIntegrity();
        testNode2.testSpatialTreeIntegrity();

        //testNode2.printTree(1);
        //std::cout << "------------------------------------------"<<std::endl;
    }
}

TEST_CASE( "Test field calculation in serial node", "[Node]") {

    BTree::Tree testTree(
            Core::Vector(-1.5, -1.5, -1.5),
            Core::Vector(1.5, 1.5, 1.5)
    );

    BTree::Node testNode(
            Core::Vector(-10.0, -10.0, -10.0),
            Core::Vector(10.0, 10.0, 10.0),
            nullptr
    );

    SECTION( "Test point to point force calculation") {
        //Test static charge calculation methods:
        Core::Vector a = Core::Vector(1.0,1.0,1.0);
        Core::Vector b = Core::Vector(2.0,1.0,1.0);

        Core::Vector c = BTree::Node::calculateElectricField(a, b, 1.0);
        Core::Vector d = Core::Vector(-(1.0/(4*M_PI*8.854e-12)),0.0,0.0);

        CHECK(
                vectorApproxCompare(
                        BTree::Node::calculateElectricField(a, b, 1.0),
                        BTree::Node::calculateElectricField(b, a, 1.0)*(-1))
                        ==  vectorsApproxEqual);

        CHECK(vectorApproxCompare(c,d) == vectorsApproxEqual);
    }

    SECTION( "Test physical correctness of charge calculation") {

        //Test particle add increases the number of particles:
        Core::Particle testIon1(Core::Vector(0.0,0.0,0.0),1.0);
        Core::Particle testIon2(Core::Vector(1.0,0.0,0.0),1.0);
        Core::Particle testIon3(Core::Vector(0.0,2.0,0.0),1.0);
        Core::Particle testIon4(Core::Vector(0.0,0.0,-1.0),100.0);
        testTree.insertParticle(testIon2,2);
        testTree.insertParticle(testIon3,3);
        testTree.insertParticle(testIon4,4);

        testTree.computeChargeDistribution();
        Core::Vector testField1 = testTree.computeEFieldFromSpaceCharge(testIon1);
        CHECK( Approx(testField1.x()).epsilon(1e-4) == -1.4399645e-9);
        CHECK( Approx(testField1.y()).epsilon(1e-4) == -1.4399645e-9/4.0);
        CHECK( Approx(testField1.z()).epsilon(1e-4) == 1.4399645e-7);
    }

    SECTION( "Test charge distribution calculation symmetrical") {

        //Test particle add increases the number of particles:
        Core::Particle testIon1 = Core::Particle(Core::Vector(0.0,0.0,0.0),1.0);
        Core::Particle testIon2 = Core::Particle(Core::Vector(-1.0,-1.0,-1.0),1.0);
        Core::Particle testIon3 = Core::Particle(Core::Vector( 1.0,-1.0,-1.0),1.0);
        Core::Particle testIon4 = Core::Particle(Core::Vector(-1.0, 1.0,-1.0),1.0);
        Core::Particle testIon5 = Core::Particle(Core::Vector( 1.0, 1.0,-1.0),1.0);
        Core::Particle testIon6 = Core::Particle(Core::Vector(-1.0,-1.0, 1.0),1.0);
        Core::Particle testIon7 = Core::Particle(Core::Vector( 1.0,-1.0, 1.0),1.0);
        Core::Particle testIon8 = Core::Particle(Core::Vector(-1.0, 1.0, 1.0),1.0);
        Core::Particle testIon9 = Core::Particle(Core::Vector( 1.0, 1.0, 1.0),1.0);

        testTree.insertParticle(testIon2,2);
        testTree.insertParticle(testIon3,3);
        testTree.insertParticle(testIon4,4);
        testTree.insertParticle(testIon5,5);
        testTree.insertParticle(testIon6,6);
        testTree.insertParticle(testIon7,7);
        testTree.insertParticle(testIon8,8);
        testTree.insertParticle(testIon9,9);

        testTree.computeChargeDistribution();

        CHECK(testTree.getRoot()->getNumberOfParticles() == 8);
        CHECK(testTree.getRoot()->getCharge() - 8.0*Core::ELEMENTARY_CHARGE < 1e-30);

        Core::Vector testField1 = testTree.computeEFieldFromSpaceCharge(testIon1);
        CHECK( ((testField1.x() < 2e-25) && (testField1.y() < 2e-25) && (testField1.z() < 2e-25)) );
    }

    SECTION( "Test charge distribution calculation non symmetrical") {

        //Test particle add increases the number of particles:
        Core::Particle baseIon1 = Core::Particle(Core::Vector(0.0,0.0,0.0),1.0);
        Core::Particle baseIon2 = Core::Particle(Core::Vector( 1.0, 1.0, 1.0),1.0);

        BTree::TreeParticle testIon1(&baseIon1);
        BTree::TreeParticle testIon2(&baseIon2);

        testNode.insertParticle(&testIon1);
        testNode.insertParticle(&testIon2);
        testNode.computeChargeDistributionRecursive();

        Core::Vector testField1 = testNode.computeElectricFieldFromTree(baseIon1);
        CHECK( testField1.x() == Approx(testField1.y()));
        CHECK( testField1.x() == Approx(testField1.z()));
    }

    SECTION( "Test charge distribution calculation in all spatial directions with few particles") {

        //Test particle add increases the number of particles:
        Core::Particle baseIon1 = Core::Particle(Core::Vector(1.0,1.0,1.0),1.0);
        Core::Particle baseIon2 = Core::Particle(Core::Vector(0.5,1.0,1.0),1.0);
        Core::Particle baseIon3 = Core::Particle(Core::Vector(1.0,0.5,1.0),1.0);
        Core::Particle baseIon4 = Core::Particle(Core::Vector(1.0,1.0,0.5),1.0);

        BTree::TreeParticle testIon1(&baseIon1);
        BTree::TreeParticle testIon2(&baseIon2);
        BTree::TreeParticle testIon3(&baseIon3);
        BTree::TreeParticle testIon4(&baseIon4);

        testNode.insertParticle(&testIon1);
        testNode.insertParticle(&testIon2);
        testNode.insertParticle(&testIon3);
        testNode.insertParticle(&testIon4);

        testNode.computeChargeDistributionRecursive();

        Core::Vector testField1 = testNode.computeElectricFieldFromTree(baseIon1);

        CHECK( testField1.x() == Approx(testField1.y()).epsilon(0.001));
        CHECK( testField1.x() == Approx(testField1.z()).epsilon(0.001));
    }

    SECTION("Test charge distribution calculation in all spatial directions with many particles") {

        Core::Particle testIon1 = Core::Particle(Core::Vector(9.0,0.0,0.0),1.0);
        Core::Particle testIon2 = Core::Particle(Core::Vector(0.0,0.0,9.0),1.0);
        Core::Particle testIon3 = Core::Particle(Core::Vector(0.0,9.0,0.0),1.0);
        Core::Particle testIon4 = Core::Particle(Core::Vector(0.0,9.0,0.0),10.0);
        std::size_t nions = 10000;

        Core::Vector boxSize = Core::Vector(2.0,2.0,2.0);
        ParticleSimulation::BoxStartZone startZone(boxSize);
        std::vector<std::unique_ptr<Core::Particle>> ions= startZone.getRandomParticlesInStartZone(nions, 2.0);
        std::vector<std::unique_ptr<BTree::TreeParticle>> treeParticles;

        for (std::size_t i=0; i<nions; i++){
            treeParticles.push_back(std::make_unique<BTree::TreeParticle>(ions[i].get()));
            testNode.insertParticle(treeParticles.back().get());
        }
        testNode.computeChargeDistributionRecursive();

        Core::Vector testField1 = testNode.computeElectricFieldFromTree(testIon1);
        Core::Vector testField2 = testNode.computeElectricFieldFromTree(testIon2);
        Core::Vector testField3 = testNode.computeElectricFieldFromTree(testIon3);
        Core::Vector testField4 = testNode.computeElectricFieldFromTree(testIon4);

        CHECK( (testField1.x() - testField2.z()) < 1e-8);
        CHECK( (testField1.x() - testField3.y()) < 1e-8);
        CHECK( (testField3 - testField4).magnitude() < 1e-8);
    }
}