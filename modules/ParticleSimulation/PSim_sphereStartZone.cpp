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

#include "PSim_sphereStartZone.hpp"
#include "Core_math.hpp"

/**
 * Constructs a spherical start zone
 *
 * @param radius Radius of the sphere
 * @param centerPosition Center position of the sphere
 */
ParticleSimulation::SphereStartZone::SphereStartZone(
    double radius, Core::Vector centerPosition) :
    radius_(radius),
    centerPosition_(centerPosition)
{

    rnd_theta_ = Core::globalRandomGeneratorPool->getUniformDistribution(0,2*M_PI);
    rnd_a_ = Core::globalRandomGeneratorPool->getUniformDistribution(-1, 1);
    rnd_R_ = Core::globalRandomGeneratorPool->getUniformDistribution(0,1);

}

/**
 * Returns a new random particle position in the spherical start zone
 */
Core::Vector ParticleSimulation::SphereStartZone::getRandomParticlePosition() {

    double phi = acos(rnd_a_->rndValue());
    double theta = rnd_theta_->rndValue();
    double R = cbrt(rnd_R_->rndValue()) * radius_;

    Core::Vector x(sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi));

    Core::Vector x_scaled(x.x()*R,x.y()*R,x.z()*R);
    Core::Vector particleCoordinates(centerPosition_.x()+x_scaled.x(),centerPosition_.y()+x_scaled.y(),centerPosition_.z()+x_scaled.z());

    return particleCoordinates;
}
