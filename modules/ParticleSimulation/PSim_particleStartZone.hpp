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

 ------------
 BTree_particleStartZone.hpp

 Abstract particle start zone class: Defines a generalized start zone for particles in the simulation

 ****************************/

#ifndef Particle_simulation_particle_start_zone
#define Particle_simulation_particle_start_zone

#include "Core_vector.hpp"
#include "Core_randomGenerators.hpp"
#include "Core_particle.hpp"
#include <memory>
#include <vector>


namespace ParticleSimulation{

    /**
     * Abstract particle start zone class
     */
    class ParticleStartZone {

    public:
        virtual ~ParticleStartZone() = default;
        virtual Core::Vector getRandomParticlePosition() = 0;

        std::vector<std::unique_ptr<Core::Particle>> getRandomParticlesInStartZone(
                std::size_t numIons, double charge, double timeOfBirthRange=0.0);
    };
}

#endif //Particle_simulation_particle_start_zone
