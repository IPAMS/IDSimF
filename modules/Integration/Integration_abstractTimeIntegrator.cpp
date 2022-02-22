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

#include "Integration_abstractTimeIntegrator.hpp"


Integration::AbstractTimeIntegrator::AbstractTimeIntegrator(particleStartMonitoringFctType ionStartMonitorFct) :
        particleStartMonitorFct_(std::move(ionStartMonitorFct))
{}

Integration::AbstractTimeIntegrator::AbstractTimeIntegrator(const std::vector<BTree::Particle*>& particles,
                                                                   particleStartMonitoringFctType ionStartMonitorFct):
        particleStartMonitorFct_(std::move(ionStartMonitorFct))
{
    //particleTOBs_ = std::vector<std::pair<double,Core::Particle*>>();
    bool tobLargerZeroFound = false;
    for (const auto &part: particles){
        if (part->getTimeOfBirth() > 0.0){
            tobLargerZeroFound = true;
        }
        particleTOBs_.emplace_back(
                part->getTimeOfBirth(),
                part);
    }

    // + sort particleTOBs_ according to time of birth
    if (tobLargerZeroFound) {
        std::sort(particleTOBs_.begin(), particleTOBs_.end(),
                [](const pTobPair_t& p1, const pTobPair_t& p2) { return p1.first<p2.first; });
    }
}

/**
 * Indicates that the time integration should be terminated at the next possible time
 */
void Integration::AbstractTimeIntegrator::setTerminationState() {
    runState_ = Integration::AbstractTimeIntegrator::IN_TERMINATION;
}

/**
 * Returns the run state of the integrator
 */
Integration::AbstractTimeIntegrator::RunState Integration::AbstractTimeIntegrator::runState() const{
    return runState_;
}


/**
 * Gets the simulated time
 */
double Integration::AbstractTimeIntegrator::time() const{
    return time_;
}

/**
 * Gets the current time step
 */
unsigned int Integration::AbstractTimeIntegrator::timeStep() const{
    return timestep_;
}

/**
 * Generate ions in the simulation which are born up to the time "time"
 * The time of birth in the particles is set according to the actual birth time in the simulation
 *
 * @param time the time until particles are born
 * @return true if particles were born
 */
bool Integration::AbstractTimeIntegrator::bearParticles_(double time){

    if (particlesBornIdx_ < particleTOBs_.size()) {
        std::size_t oldParticlesBornIdx = particlesBornIdx_;

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