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

#include "PSim_idealizedQitFFTWriter.hpp"
#include "Core_vector.hpp"
#include "BTree_particle.hpp"
#include <map>

/**
 * Constructor: Creates a filewriter with a given particle cloud and a file to write to
 * @param particles vector with links to the partilces in a particle cloud to write to a file
 * @param transientFilename the filename of the file to write to
 */
ParticleSimulation::IdealizedQitFFTWriter::IdealizedQitFFTWriter(std::vector<BTree::Particle*> particles,
                                                                 std::string transientFilename)
{
    particles_ = particles;
    transientFile_ = std::make_unique<std::ofstream>();
    transientFile_->open(transientFilename);
}

/**
 * Destructor
 */
ParticleSimulation::IdealizedQitFFTWriter::~IdealizedQitFFTWriter(){
    transientFile_->flush();
}


/**
 * Writes a time step to the file: The current state of the particles in the particle cloud is used to
 * determine the average velocity in z direction for all particles and the result is written to
 * the result file.
 * @param time the time of the current time step
 */
void ParticleSimulation::IdealizedQitFFTWriter::writeTimestep(double time){
    double avgVelocityZ = 0;

    for (const auto& part: particles_) {
        avgVelocityZ = avgVelocityZ + part->getVelocity().z();
    }
    avgVelocityZ = avgVelocityZ / particles_.size();
    *transientFile_ << time << " " << avgVelocityZ << std::endl;
}

/**
 * Writes a time step to the file in a mass resolved way: The current state of the particles in the particle cloud is
 * used to determine the average velocity in z direction individually for groups containing particles with the same mass.
 * The result is written to the result file.
 * @param time the time of the current time step
 */
void ParticleSimulation::IdealizedQitFFTWriter::writeTimestepMassResolved(double time){
    std::map<double,double> avgVelocityZ;
    std::map<double,int> nParticles;

    for (const auto& part: particles_) {
        double mass = part->getMass();
        if (avgVelocityZ.count(mass) == 0){
            avgVelocityZ[mass] = 0;
            nParticles[mass] = 0;
        }
        avgVelocityZ[mass] = avgVelocityZ[mass] + part->getVelocity().z();
        nParticles[mass] += 1;
    }

    *transientFile_ << time <<" ";
    for ( const auto &pair : avgVelocityZ) {
        auto key = pair.first;
        avgVelocityZ[key] = avgVelocityZ[key] / nParticles[key];

        *transientFile_ << avgVelocityZ[key] <<" ";
    }
    *transientFile_<< std::endl;
}

