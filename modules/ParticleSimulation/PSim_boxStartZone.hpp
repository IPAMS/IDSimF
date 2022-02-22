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
 PSim_BoxStartZone.hpp

 Box particle start zone

 ****************************/
#ifndef Particle_simulation_box_start_zone
#define Particle_simulation_box_start_zone

#include "Core_vector.hpp"
#include "PSim_particleStartZone.hpp"

namespace ParticleSimulation{

    /**
     * A box start zone with faces parallel to the axes of the coordinate system
     */
    class BoxStartZone : public ParticleStartZone {

    public:
        explicit BoxStartZone(Core::Vector size, Core::Vector centerPosition = {0.0, 0.0, 0.0});

        [[nodiscard]] Core::Vector getRandomParticlePosition() override;

    private:
        Core::RndDistPtr rnd_x; ///< random distribution of positions in x direction
        Core::RndDistPtr rnd_y; ///< random distribution of positions in x direction
        Core::RndDistPtr rnd_z; ///< random distribution of positions in x direction
    };
}

#endif //Particle_simulation_box_start_zone
