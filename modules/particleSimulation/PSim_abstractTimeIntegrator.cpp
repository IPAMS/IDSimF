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

#include "PSim_abstractTimeIntegrator.hpp"


ParticleSimulation::AbstractTimeIntegrator::AbstractTimeIntegrator() :
time_(0.0),
timestep_(0),
particles_(std::vector<BTree::Particle *>()),
nParticles_(0),
particleTOBs_(std::vector<pTobPair_t>()),
particlesBornIdx_(0)
{}

ParticleSimulation::AbstractTimeIntegrator::AbstractTimeIntegrator(std::vector<BTree::Particle*> particles):
time_(0.0),
timestep_(0),
particles_(std::vector<BTree::Particle *>()),
nParticles_(0),
particleTOBs_(std::vector<pTobPair_t>()),
particlesBornIdx_(0)
{
    //init velocities and accelerations:

    //particleTOBs_ = std::vector<std::pair<double,Core::Particle*>>();
    for (const auto &part: particles){
        particleTOBs_.emplace_back(
                part->getTimeOfBirth(),
                part);
    }

    // + sort particleTOBs_ according to time of birth
    std::sort(particleTOBs_.begin(), particleTOBs_.end(),
            [](const pTobPair_t &p1, const pTobPair_t &p2) {return p1.first < p2.first;});
}

/**
 * Generate ions in the simulation which are born up to the time "time"
 * The time of birth in the particles is set according to the actual birth time in the simulation
 *
 * @param time the time until particles are born
 */
void ParticleSimulation::AbstractTimeIntegrator::bearParticles_(double time){

    if (particlesBornIdx_ < particleTOBs_.size()) {
        while (particlesBornIdx_ < particleTOBs_.size() && particleTOBs_[particlesBornIdx_].first <= time) {
            BTree::Particle *part = particleTOBs_[particlesBornIdx_].second;
            part->setTimeOfBirth(time); //set particle TOB to the actual TOB in the simulation
            addParticle(part);
            ++particlesBornIdx_;
        }
    }
}