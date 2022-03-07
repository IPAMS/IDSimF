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

 ------------
 test_particleStarting.hpp

 Description

 ****************************/
#ifndef IDSIMF_TEST_PARTICLEINIT_HPP
#define IDSIMF_TEST_PARTICLEINIT_HPP

#include "PSim_boxStartZone.hpp"


//FIXME: get rid of this method: Rewrite tests with start zone directly
inline std::vector<std::unique_ptr<Core::Particle>> getRandomIonsInBox(std::size_t numIons, Core::Vector corner, Core::Vector boxSize){
    Core::Vector center = corner + boxSize/2.0;
    ParticleSimulation::BoxStartZone startZone(boxSize, center);

    return startZone.getRandomParticlesInStartZone(numIons, 1.0);
}

inline std::vector<std::unique_ptr<Core::Particle>> getIonsInLattice(unsigned int nPerDirection){

    std::vector<std::unique_ptr<Core::Particle>> particles;
    unsigned int nTotal = 0;
    for (unsigned int i=0; i<nPerDirection; i++){
        double pX = i*1.0/nPerDirection;
        for (unsigned int j=0; j<nPerDirection; j++){
            double pY = j*1.0/nPerDirection;
            for (unsigned int k=0; k<nPerDirection; k++){
                double pZ = k*1.0/nPerDirection;
                std::unique_ptr<Core::Particle> newIon = std::make_unique<Core::Particle>(Core::Vector(pX,pY,pZ), 1.0);
                newIon -> setMassAMU(100);
                particles.push_back(std::move(newIon));
                nTotal++;
            }
        }
    }
    return particles;
}

#endif //IDSIMF_TEST_PARTICLEINIT_HPP
