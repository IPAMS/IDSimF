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
#include "SC_fullSumSolver.hpp"
#include "BTree_abstractNode.hpp"
#include <algorithm>

Core::Vector SpaceCharge::FullSumSolver::getEFieldFromSpaceCharge(Core::Particle& particle) {
    Core::Vector particleLocation = particle.getLocation();
    std::size_t nParticles = getNumberOfParticles();

    std::vector<const particleListEntry*> particleVector;
    std::transform((*iVec_).begin(), (*iVec_).end(), std::back_inserter(particleVector),
            [](const particleListEntry& pe) -> const particleListEntry* { return &pe; });


    Core::Vector sum(0.0, 0.0, 0.0);
    #pragma omp parallel for reduction(+:sum) default(none) shared(particleLocation, nParticles, particleVector)
    for (std::size_t i=0; i<nParticles; ++i){
        Core::Vector field = BTree::AbstractNode::calculateElectricField(
                particleLocation,
                particleVector[i]->particle->getLocation(),
                particleVector[i]->particle->getCharge()
                );
        sum += field;
    }
    return sum;
}

void SpaceCharge::FullSumSolver::computeChargeDistribution() {}
/*
void SPa::computeChargeDistribution() {
    std::size_t nParticles = getNumberOfParticles();

    std::vector<double> sources(3*nParticles);
    std::vector<double> charges(nParticles);
    std::vector<double> potentials(nParticles);
    std::vector<double> gradients(3*nParticles);

    std::size_t i = 0;
    Core::Vector particlePos;
    for (const SpaceCharge::particleListEntry& pListEntry : *iVec_) {
        particlePos = pListEntry.particle -> getLocation();
        sources[3*i] = particlePos.x();
        sources[3*i+1] = particlePos.y();
        sources[3*i+2] = particlePos.z();

        charges[i] = pListEntry.particle -> getCharge();
        i++;
    }

    //double eps = 0.5e-6;
    double eps = 0.5e-5;

    // call the actual fmm routine
    int ier =0;
    int nP = nParticles;
    lfmm3d_s_c_g_wrapper(&eps, &nP, sources.data(), charges.data(),
            potentials.data(), gradients.data(), &ier);

    if (ier != 0){
        throw (std::runtime_error("FMM Calculation failed"));
    }

    i = 0;
    for (SpaceCharge::particleListEntry& pListEntry : *iVec_) {
        pListEntry.potential = potentials[i];
        pListEntry.gradient = {gradients[3*i], gradients[3*i+1], gradients[3*i+2]};
        i++;
    }
}
*/