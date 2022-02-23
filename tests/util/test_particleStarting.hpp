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
#ifndef IDSIMF_TEST_PARTICLESTARTING_HPP
#define IDSIMF_TEST_PARTICLESTARTING_HPP

#include "PSim_boxStartZone.hpp"


//FIXME: get rid of this method: Rewrite tests with start zone directly
inline std::vector<std::unique_ptr<Core::Particle>> getRandomIonsInBox(std::size_t numIons, Core::Vector corner, Core::Vector boxSize){
    Core::Vector center = corner + boxSize/2.0;
    ParticleSimulation::BoxStartZone startZone(boxSize, center);

    return startZone.getRandomParticlesInStartZone(numIons, 1.0);
}

#endif //IDSIMF_TEST_PARTICLESTARTING_HPP
