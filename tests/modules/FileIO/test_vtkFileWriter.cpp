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
 test_vtkFieldWriter.cpp

 Testing of vtk field file writer

 ****************************/

#include "PSim_simpleVTKwriter.hpp"
#include "BTree_particle.hpp"
#include "BTree_tree.hpp"

#include "catch.hpp"
#include <iostream>

TEST_CASE( "The VTK writer should at least write a file without exception",
        "[ParticleSimulation][VtkFieldReader][file writers]") {

    //Test file writer:
    Core::Vector min = Core::Vector(0.0,0.0,0.0);
    Core::Vector max = Core::Vector(2.0,2.0,2.0);
    BTree::Tree* testTree = new BTree::Tree(min,max);
    
    
    BTree::Particle testIon1 = BTree::Particle(Core::Vector(1.2,1.2,1.2),2.0);
    BTree::Particle testIon2 = BTree::Particle(Core::Vector(1.6,1.2,1.2),2.0);
    BTree::Particle testIon3 = BTree::Particle(Core::Vector(1.65,1.2,1.2),2.0);
    BTree::Particle testIon4 = BTree::Particle(Core::Vector(1.66,1.2,1.2),2.0);
    BTree::Particle testIon5 = BTree::Particle(Core::Vector(1.665,1.2,1.2),2.0);

    testTree->insertParticle(testIon1,1);
    testTree->insertParticle(testIon2,2);
    testTree->insertParticle(testIon3,3);
    testTree->insertParticle(testIon4,4);
    testTree->insertParticle(testIon5,5);
    //root->printNode(0);

    FileIO::SimpleVTKwriter* testWriter = new FileIO::SimpleVTKwriter("test");
    testWriter->write(*testTree,true);
    delete(testTree);

    CHECK(testIon1.getLocation() != testIon2.getLocation());
}
