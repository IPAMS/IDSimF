/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2022 - Physical and Theoretical Chemistry /
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

#include "SC_generic.hpp"


/**
 * Default Constructor
 */
SpaceCharge::GenericSpaceChargeSolver::GenericSpaceChargeSolver() {
    iVec_ = std::make_unique<particlePtrList>();
    iMap_ = std::make_unique<std::unordered_map<std::size_t, particlePtrList::const_iterator>>();
    pMap_ = std::make_unique<std::unordered_map<Core::Particle*, particlePtrList::const_iterator>>();
}

/**
 Insert particle into the field solver

 \param particle the particle to insert
 \param ext_index an external index number for the particle / numerical particle id (most likely from simion)
 */
void SpaceCharge::GenericSpaceChargeSolver::insertParticle(Core::Particle &particle, std::size_t ext_index){
    iVec_->push_front({&particle, Core::Vector(), 0.0});
    iMap_->insert({ext_index, iVec_->cbegin()});
    pMap_->insert({&particle, iVec_->cbegin()});
}

/**
 Removes a particle with a given particle index / particle id and its hosting leaf node
 from the fmm solver

 \param ext_index the external numerical particle id
 */
void SpaceCharge::GenericSpaceChargeSolver::removeParticle(std::size_t ext_index){
    auto iter =(*iMap_)[ext_index];
    Core::Particle* particle = iter->particle;
    iMap_->erase(ext_index);
    pMap_->erase(particle);
    iVec_->erase(iter);
}

/**
 Gets the number of particles in the fmm solver
 */
std::size_t SpaceCharge::GenericSpaceChargeSolver::getNumberOfParticles() const{
    return(iVec_->size());
}

Core::Vector SpaceCharge::GenericSpaceChargeSolver::getEFieldFromSpaceCharge(Core::Particle& particle) {
    auto iter =(*pMap_)[&particle];
    return iter->gradient/NEGATIVE_ELECTRIC_CONSTANT;
}

