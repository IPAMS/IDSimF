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

#include "PSim_boxStartZone.hpp"
#include "Core_math.hpp"

/**
 * Constructs a box start zone
 *
 * @param size Size of the box in x,y,z, direction
 * @param centerPosition Position of the box center
 */
ParticleSimulation::BoxStartZone::BoxStartZone(Core::Vector size, Core::Vector centerPosition) {

    Core::Vector cornerLower = centerPosition - (size*0.5);
    Core::Vector cornerUpper = cornerLower + size;

    rnd_x = Core::globalRandomGeneratorPool->getUniformDistribution(cornerLower.x(), cornerUpper.x());
    rnd_y = Core::globalRandomGeneratorPool->getUniformDistribution(cornerLower.y(), cornerUpper.y());
    rnd_z = Core::globalRandomGeneratorPool->getUniformDistribution(cornerLower.z(), cornerUpper.z());
}

/**
 * Gets a random position in the box particle start zone.
 */
Core::Vector ParticleSimulation::BoxStartZone::getRandomParticlePosition() {
    return {rnd_x->rndValue(), rnd_y->rndValue(), rnd_z->rndValue()};
}
