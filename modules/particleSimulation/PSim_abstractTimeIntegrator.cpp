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


ParticleSimulation::AbstractTimeIntegrator::AbstractTimeIntegrator(particleStartMonitoringFctType ionStartMonitorFct) :
        particleStartMonitorFct_(std::move(ionStartMonitorFct))
{}

ParticleSimulation::AbstractTimeIntegrator::AbstractTimeIntegrator(std::vector<BTree::Particle*> particles,
                                                                   particleStartMonitoringFctType ionStartMonitorFct):
        particleStartMonitorFct_(std::move(ionStartMonitorFct))
{
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
 * Indicates that the time integration should be terminated at the next possible time
 */
void ParticleSimulation::AbstractTimeIntegrator::setTerminationState() {
    runState_ = ParticleSimulation::AbstractTimeIntegrator::IN_TERMINATION;
}


/**
 * Gets the simulated time
 */
double ParticleSimulation::AbstractTimeIntegrator::time() {
    return time_;
}

/**
 * Gets the current time step
 */
int ParticleSimulation::AbstractTimeIntegrator::timeStep() {
    return timestep_;
}

/**
 * Generate ions in the simulation which are born up to the time "time"
 * The time of birth in the particles is set according to the actual birth time in the simulation
 *
 * @param time the time until particles are born
 * @return true if particles were born
 */
bool ParticleSimulation::AbstractTimeIntegrator::bearParticles_(double time){

    if (particlesBornIdx_ < particleTOBs_.size()) {
        int oldParticlesBornIdx = particlesBornIdx_;

        while (particlesBornIdx_ < particleTOBs_.size() && particleTOBs_[particlesBornIdx_].first <= time) {
            BTree::Particle *part = particleTOBs_[particlesBornIdx_].second;
            part->setTimeOfBirth(time); //set particle TOB to the actual TOB in the simulation
            addParticle(part);
            if (particleStartMonitorFct_ != nullptr){
                particleStartMonitorFct_(part, time);
            }
            ++particlesBornIdx_;
        }

        return particlesBornIdx_ > oldParticlesBornIdx;
    } else{
        return false;
    }
}