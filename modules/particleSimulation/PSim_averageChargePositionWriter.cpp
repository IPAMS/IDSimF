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
 ****************************/

#include "PSim_averageChargePositionWriter.hpp"
#include "BTree_tree.hpp"
#include "BTree_particle.hpp"
#include <vector>

/**
 * Constructor
 * @param transientFilename a filename to write to
 */
ParticleSimulation::AverageChargePositionWriter::AverageChargePositionWriter(std::string transientFilename){
    transientFile_ = new std::ofstream();
    transientFile_->open(transientFilename);
}

/**
 * Destructor
 */
ParticleSimulation::AverageChargePositionWriter::~AverageChargePositionWriter(){
    if(transientFile_ != nullptr){
        transientFile_->flush();
        delete (transientFile_);
    }
}

/**
 * Writes a timestep to the file
 * @param tree a valid Core which contains the charged particles
 * @param time the time of the current time step
 */
void ParticleSimulation::AverageChargePositionWriter::writeTimestep(BTree::Tree &tree, double time){
    *transientFile_<<time<<" "<<tree.getRoot()->getCenterOfCharge()<<std::endl;
}
