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

#include "PSim_inductionCurrentWriter.hpp"
#include "Core_vector.hpp"
#include "BTree_particle.hpp"
#include "PSim_simionPotentialArray.hpp"

/**
 * Constructor: Creates a induction current writer
 *
 * @param particles vector with links to the partilces in a particle cloud to write to a file
 * @param transientFilename the filename of the file to write to
 * @param weightFields a vector of the weight potential arrays for the individual detection electrodes. A weight field
 * of a detection electrode is a potential array in which the detection electrode has 1 V potential while all other
 * electrodes are on ground
 * @param weightFactors scaling factors for the individual weight fields
 * @param scale_mm_per_gu scaling factor between mm and grid units of the given weight potential arrays
 */
FileIO::InductionCurrentWriter::InductionCurrentWriter(
        std::vector<BTree::Particle*> particles,
        std::string transientFilename,
        const std::vector<ParticleSimulation::SimionPotentialArray*>& weightFields,
        std::vector<double> weightFactors,
        double scale_mm_per_gu) {
    scale_mm_per_gu_ = scale_mm_per_gu;
    particles_ = particles;

    std::size_t nFields = weightFields.size();
    if (nFields != weightFactors.size()){
        throw std::invalid_argument("Field and field factor vectors differ in size");
    }

    for (std::size_t i=0; i<nFields; ++i){
        weightFields_.emplace_back(
                std::make_pair(weightFields[i],weightFactors[i]));
    }

    transientFile_ = std::make_unique<std::ofstream>();
    transientFile_->open(transientFilename);
}

/**
 * Destructor
 */
FileIO::InductionCurrentWriter::~InductionCurrentWriter(){
    transientFile_->flush();
}


/**
 * Writes a time step to the file
 * @param time the time of the current time step
 */
void FileIO::InductionCurrentWriter::writeTimestep(double time){
    double inducedCurrent = 0;

    // iterate through all particles:
    for (const auto& part: particles_) {
        if (part->isActive()) {
            Core::Vector pos = part->getLocation();
            Core::Vector velocity = part->getVelocity();
            Core::Vector weightField(0, 0, 0);

            //get local weigh field (interpolated value of the weight PA * weight factor:
            for (std::pair<ParticleSimulation::SimionPotentialArray*, double> &wPA: weightFields_) {
                weightField = weightField +
                        (wPA.first->getField(pos.x(), pos.y(), pos.z()) * wPA.second);
            }
            weightField = weightField * 1.0 / scale_mm_per_gu_;

            //calculate induction current of this particle:
            double partInducedCurrent = part->getCharge() *
                                        (velocity.x() * weightField.x() +
                                         velocity.y() * weightField.y() +
                                         velocity.z() * weightField.z());

            inducedCurrent += partInducedCurrent;
        }
    }

    *transientFile_ << time << " " << inducedCurrent << std::endl;
}