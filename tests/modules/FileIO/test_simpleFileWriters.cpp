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
 test_simpleFileWriters.cpp

 Testing of simple / special purpose file writers

 ****************************/

#include "BTree_particle.hpp"
#include "BTree_tree.hpp"
#include "FileIO_scalar_writer.hpp"
#include "FileIO_averageChargePositionWriter.hpp"
#include "FileIO_idealizedQitFFTWriter.hpp"
#include "catch.hpp"
#include <memory>


TEST_CASE( "Scalar file writer should at least write a file without exception",
        "[ParticleSimulation][ScalarWriter][file writers]") {
    FileIO::Scalar_writer fw("scalar_file_writer_test.txt");

    fw.writeTimestep(10,0.1);
    fw.writeTimestep(20,0.2);
    fw.writeTimestep(1.23,0.3);
    fw.writeTimestep(0.99,10);
}

TEST_CASE( "Average ion position file writer should at least write a file without exception",
        "[ParticleSimulation][AverageChargePositionFileWriter][file writers]") {

    Core::Vector min = Core::Vector(-20.0,-20.0,-20.0);
    Core::Vector max = Core::Vector(20.0,20.0,20.0);
    BTree::Tree tree(min,max);

    std::vector<BTree::uniquePartPtr> particles;
    std::vector<BTree::Particle*> particlePtrs;
    unsigned int nParticles = 10;
    for (std::size_t i=0; i< nParticles; i++){
        std::unique_ptr<BTree::Particle> part = std::make_unique<BTree::Particle>(Core::Vector(0,0,static_cast<double>(i)*0.1),1);
        particlePtrs.push_back(part.get());
        particles.push_back(std::move(part));
        tree.insertParticle(*(particlePtrs[i]), i);
    }

    FileIO::AverageChargePositionWriter fw("average_position_test.txt");

    unsigned int nSteps = 10;
    double diff =0.1;
    for (unsigned int i=0; i< nSteps; ++i){
        for (std::size_t k=0; k < nParticles; ++k){
            Core::Vector newLocation = particlePtrs[k]->getLocation() + Core::Vector(diff,diff,diff);
            tree.updateParticleLocation(k, newLocation);
        }
        fw.writeTimestep(tree,0.1*i);
    }
}

TEST_CASE("Test FFT file writer", "[ParticleSimulation][ScalarWriter][file writers]") {

    std::vector<std::unique_ptr<BTree::Particle>> particles;
    std::vector<BTree::Particle*> particlePtrs;
    unsigned int nParticles_a = 10;
    unsigned int nParticles_b = 10;
    unsigned int nParticles = nParticles_a + nParticles_b;
    for (unsigned int i=0; i < nParticles_a; i++){
        std::unique_ptr<BTree::Particle> part = std::make_unique<BTree::Particle>(Core::Vector(0,0,i*0.1),1);
        part->setMassAMU(10);
        particlePtrs.push_back(part.get());
        particles.push_back(std::move(part));
    }
    for (unsigned int i=0; i< nParticles_b; i++){
        std::unique_ptr<BTree::Particle> part = std::make_unique<BTree::Particle>(Core::Vector(0.1,0,i*0.1),1);
        part->setMassAMU(20);
        particlePtrs.push_back(part.get());
        particles.push_back(std::move(part));
    }


    SECTION("Simple FFT file writer should at least write a file without exception when using non mass resolved mode"){

        FileIO::IdealizedQitFFTWriter fw(particlePtrs, "fft_file_writer_test.txt");

        unsigned int nSteps = 10;
        double diff = 0;
        for (unsigned int i=0; i< nSteps; i++){
            if (i% 2 ==0){
                diff = -1;
            }
            else{
                diff = 2;
            }
            for (std::size_t k=0; k < nParticles; k++){
                particles[i]->setVelocity(Core::Vector(0,0,diff));

            }
            fw.writeTimestep(0.1*i);
        }
    }

    SECTION("Simple Fft file writer should at least write a file without exception when using mass resolved mode"){

        FileIO::IdealizedQitFFTWriter fw(particlePtrs, "fft_file_writer_test_mass_resolved.txt");

        int nSteps = 10;
        double diff = 0;
        for (int i=0; i< nSteps; i++){
            if (i% 2 ==0){
                diff = -1;
            }
            else{
                diff = 2;
            }
            for (std::size_t k=0; k < nParticles_a; k++){
                particles[k]->setVelocity(Core::Vector(0, 0, diff));
            }
            for (std::size_t k=nParticles_a; k < nParticles; k++){
                particles[k]->setVelocity(Core::Vector(0, 0, 2.0*diff));
            }
            fw.writeTimestepMassResolved(0.1*i);
        }
    }
}
