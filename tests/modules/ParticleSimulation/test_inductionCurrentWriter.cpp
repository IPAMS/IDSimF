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
 test_inductionCurrentWriter.cpp

 Testing of file writer for electrode induction / mirror currents

 ****************************/

#include "BTree_particle.hpp"
#include "PSim_simionPotentialArray.hpp"
#include "PSim_inductionCurrentWriter.hpp"
#include "catch.hpp"
#include <memory>


TEST_CASE( "Test indcution current file writer", "[ParticleSimulation][file writers]") {

    SECTION("Induction current writer should be able to write file"){

        std::vector<std::unique_ptr<BTree::Particle>> particles;
        std::vector<BTree::Particle*> particlePtrs;
        unsigned int nParticles = 10;
        for (unsigned int i=0; i < nParticles; i++){
            std::unique_ptr<BTree::Particle> part = std::make_unique<BTree::Particle>(Core::Vector(0.01,0.01,i*0.001),1);
            particlePtrs.push_back(part.get());
            particles.push_back(std::move(part));
        }

        ParticleSimulation::SimionPotentialArray simPa1("cylinder_capacitor.pa1");
        ParticleSimulation::SimionPotentialArray simPa2("cylinder_capacitor.pa2");
        std::vector<ParticleSimulation::SimionPotentialArray*> detectionPAs= {&simPa1, &simPa2};
        std::vector<double> detectionScalingFactors = {1.0,-1.0};

        ParticleSimulation::InductionCurrentWriter inductionWriter(particlePtrs, "induction_current_file_writer_test.txt",
                detectionPAs,detectionScalingFactors, 1.0);

        int nSteps = 15;
        double diff =0;
        for (int i=0; i< nSteps; i++){
            if (i% 2 ==0){
                diff = -0.001;
            }
            else{
                diff = 0.002;
            }
            for (std::size_t k=0; k < nParticles; k++){
                particles[k]->setVelocity(Core::Vector(diff,0,0));

            }
            inductionWriter.writeTimestep(0.1*i);
        }
    }
}