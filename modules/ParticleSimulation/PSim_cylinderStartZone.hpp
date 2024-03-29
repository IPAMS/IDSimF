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
 PSim_cylinderStartZone.hpp

 Cylindrical particle start zone

 ****************************/
#ifndef Particle_simulation_cylinder_start_zone
#define Particle_simulation_cylinder_start_zone

#include "Core_vector.hpp"
#include "PSim_particleStartZone.hpp"

namespace ParticleSimulation{

    /**
     * Cylindrical particle start zone which can be rotated and translated in space.
     */
    class CylinderStartZone : public ParticleStartZone {

    public:
        CylinderStartZone(double radius, double length,
                          Core::Vector normalVector = {0.0, 0.0, 0.0},
                          Core::Vector baseVector = {0.0, 0.0, 0.0});

        [[nodiscard]] Core::Vector getRandomParticlePosition() override;

    private:
        double radius_ = 0.0;
        double length_ = 0.0;

        bool isRotated_ = false;
        bool isShifted_ = false;
        double azimuth_ = 0.0;
        double elevation_ = 0.0;
        Core::Vector baseVector_;

        Core::RndDistPtr rnd_x_;
        Core::RndDistPtr rnd_R_ = Core::globalRandomGeneratorPool->getUniformDistribution(0,1);
        Core::RndDistPtr rnd_phi_ = Core::globalRandomGeneratorPool->getUniformDistribution(0,2*M_PI);
    };
}

#endif //Particle_simulation_cylinder_start_zone
