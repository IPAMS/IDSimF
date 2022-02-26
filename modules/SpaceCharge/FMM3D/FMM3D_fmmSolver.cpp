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

#include "FMM3D_fmmSolver.hpp"
extern "C" {
    #include "FMM3D_fmmSolver_C_interface.h"
}
#include "Core_randomGenerators.hpp"
#include <vector>
#include <iostream>

/**
 * Default Constructor
 */
FMM3D::FMMSolver::FMMSolver() {
    iVec_ = std::make_unique<particlePtrList>();
    iMap_ = std::make_unique<std::unordered_map<std::size_t, particlePtrList::const_iterator>>();
    pMap_ = std::make_unique<std::unordered_map<Core::Particle*, particlePtrList::const_iterator>>();
}

/**
 Insert particle into the field solver

 \param particle the particle to insert
 \param ext_index an external index number for the particle / numerical particle id (most likely from simion)
 */
void FMM3D::FMMSolver::insertParticle(Core::Particle &particle, std::size_t ext_index){
    iVec_->push_front({&particle, Core::Vector(), 0.0});
    iMap_->insert({ext_index, iVec_->cbegin()});
    pMap_->insert({&particle, iVec_->cbegin()});
}

/**
 Removes a particle with a given particle index / particle id and its hosting leaf node
 from the fmm solver

 \param ext_index the external numerical particle id
 */
void FMM3D::FMMSolver::removeParticle(std::size_t ext_index){
    auto iter =(*iMap_)[ext_index];
    Core::Particle* particle = iter->particle;
    iMap_->erase(ext_index);
    pMap_->erase(particle);
    iVec_->erase(iter);
}

/**
 Gets the number of particles in the fmm solver
 */
std::size_t FMM3D::FMMSolver::getNumberOfParticles() const{
    return(iVec_->size());
}

Core::Vector FMM3D::FMMSolver::getEFieldFromSpaceCharge(Core::Particle& particle) {
    auto iter =(*pMap_)[&particle];
    return iter->gradient;
}


void FMM3D::FMMSolver::computeChargeDistribution() {
    std::size_t nParticles = getNumberOfParticles();

    std::vector<double> sources(3*nParticles);
    std::vector<double> charges(nParticles);
    std::vector<double> potentials(nParticles);
    std::vector<double> gradients(3*nParticles);

    std::size_t i = 0;
    Core::Vector particlePos;
    for (const FMM3D::particleListEntry& pListEntry : *iVec_) {
        particlePos = pListEntry.particle -> getLocation();
        sources[3*i] = particlePos.x();
        sources[3*i+1] = particlePos.y();
        sources[3*i+2] = particlePos.z();

        charges[i] = pListEntry.particle -> getCharge();
        i++;
    }

    double eps = 0.5e-6;

    // call the actual fmm routine
    int ier =0;
    int nP = nParticles;
    lfmm3d_s_c_g_wrapper(&eps, &nP, sources.data(), charges.data(),
                         potentials.data(), gradients.data(), &ier);

    i = 0;
    for (FMM3D::particleListEntry& pListEntry : *iVec_) {
        pListEntry.potential = potentials[i];
        pListEntry.gradient = {gradients[3*i], gradients[3*i+1], gradients[3*i+2]};
        i++;
    }
}