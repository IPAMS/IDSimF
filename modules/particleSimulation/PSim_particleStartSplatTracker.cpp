/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2021 - Physical and Theoretical Chemistry /
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

#include "PSim_particleStartSplatTracker.hpp"
#include "BTree_particle.hpp"
#include <algorithm>

ParticleSimulation::ParticleStartSplatTracker::ParticleStartSplatTracker():
    pMap_(),
    pInsertIndex_(0)
{}

void ParticleSimulation::ParticleStartSplatTracker::particleStart(BTree::Particle* particle, double time) {

    // Existing particle entries are silently overwritten:
    pMapEntry entry;
    entry.globalIndex = pInsertIndex_;
    entry.startLocation = particle->getLocation();
    entry.startTime = time;
    entry.state = STARTED;

    auto emplaceResult = pMap_.try_emplace(particle, entry);
    if (emplaceResult.second){ //emplace was successful, key was not already existing
        particle->setIntegerAttribute("global index", (int) entry.globalIndex);
        pInsertIndex_++;
    } else {
        throw (std::invalid_argument("Illegal double insert into start splat tracker: Particle is already existing"));
    }

}

void ParticleSimulation::ParticleStartSplatTracker::particleRestart(BTree::Particle* particle,
                                                                    Core::Vector oldPosition, Core::Vector newPosition,
                                                                    double time) {
    try{
        pMapEntry& entry = pMap_.at(particle);
        entry.splatLocation = oldPosition;
        entry.splatTime = time;
        entry.state = SPLATTED_AND_RESTARTED;
        restartedParticlesData_.emplace_back(entry);

        pMapEntry newEntry;
        newEntry.globalIndex = pInsertIndex_;
        newEntry.startLocation = newPosition;
        newEntry.startTime = time;
        newEntry.state = RESTARTED;

        pMap_.at(particle) = newEntry;
        particle->setIntegerAttribute("global index", (int) newEntry.globalIndex);
        pInsertIndex_++;
    }
    catch (const std::out_of_range& exception) {
        throw (std::invalid_argument("Particle to restart was not registered as started before"));
    }
}


void ParticleSimulation::ParticleStartSplatTracker::particleSplat(BTree::Particle* particle, double time) {
    try{
        pMapEntry& entry = pMap_.at(particle);
        entry.splatLocation = particle->getLocation();
        entry.splatTime = time;
        entry.state = SPLATTED;
    }
    catch (const std::out_of_range& exception) {
        throw (std::invalid_argument("Particle to splat was not registered as started before"));
    }
}

ParticleSimulation::ParticleStartSplatTracker::pMapEntry
        ParticleSimulation::ParticleStartSplatTracker::get(BTree::Particle* particle) const{
    return pMap_.at(particle);
}

void ParticleSimulation::ParticleStartSplatTracker::sortStartSplatData() {
    sortedParticleData_ = std::vector<pMapEntry>();


    for (auto const &mapEntry : pMap_){
        sortedParticleData_.emplace_back(mapEntry.second);
    }

    for (auto const &entry : restartedParticlesData_){
        sortedParticleData_.emplace_back(entry);
    }

    // Sort according to global index:
    std::sort(sortedParticleData_.begin(), sortedParticleData_.end(),
            [](const pMapEntry &e1, const pMapEntry &e2) {return e1.globalIndex < e2.globalIndex;});
}


std::vector<ParticleSimulation::ParticleStartSplatTracker::pMapEntry>
        ParticleSimulation::ParticleStartSplatTracker::getStartSplatData() const{
    return sortedParticleData_;
}

std::vector<int> ParticleSimulation::ParticleStartSplatTracker::getSplatState() const{
    std::vector<int> result;
    for (auto const &entry: sortedParticleData_){
        result.emplace_back(entry.state);
    }
    return result;
}

std::vector<double> ParticleSimulation::ParticleStartSplatTracker::getStartTimes() const{
    std::vector<double> result;
    for (auto const &entry: sortedParticleData_){
        result.emplace_back(entry.startTime);
    }
    return result;
}

std::vector<double> ParticleSimulation::ParticleStartSplatTracker::getSplatTimes() const{
    std::vector<double> result;
    for (auto const &entry: sortedParticleData_){
        result.emplace_back(entry.splatTime);
    }
    return result;
}

std::vector<Core::Vector> ParticleSimulation::ParticleStartSplatTracker::getStartLocations() const{
    std::vector<Core::Vector> result;
    for (auto const &entry: sortedParticleData_){
        result.emplace_back(entry.startLocation);
    }
    return result;
}

std::vector<Core::Vector> ParticleSimulation::ParticleStartSplatTracker::getSplatLocations() const{
    std::vector<Core::Vector> result;
    for (auto const &entry: sortedParticleData_){
        result.emplace_back(entry.splatLocation);
    }
    return result;
}