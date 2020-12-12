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

#include "PSim_cylinderStartZone.hpp"
#include "Core_math.hpp"

/**
 * Constructs a cylinder start zone
 *
 * @param radius radius of the cylinder
 * @param length length of the cylinder
 * @param normalVector normal vector of the
 */
ParticleSimulation::CylinderStartZone::CylinderStartZone(
    double radius, double length, Core::Vector normalVector, Core::Vector baseVector) :
radius_(radius),
length_(length),
baseVector_(baseVector)
{
    if (normalVector.magnitude() > 0.0){
        normalVector_ = normalVector/normalVector.magnitude();
        isRotated_ = true;
        Core::Vector normalPolarCoords = Core::cartesianToPolar(normalVector);
        azimuth_ = normalPolarCoords.y();
        elevation_ = normalPolarCoords.z();
    }

    if (baseVector.magnitude() > 0.0){
        isShifted_ = true;
    }

    rnd_x_ = Core::globalRandomGenerator->getUniformDistribution(0.0, length_);
}

Core::Vector ParticleSimulation::CylinderStartZone::getRandomParticlePosition() {

    double R = sqrt(rnd_R_->rndValue()) * radius_;
    double phi = rnd_phi_->rndValue();

    // base cartesian coordinates in base cylinder (in x-direction)
    Core::Vector particleCoordinates(rnd_x_->rndValue(), sin(phi)*R, cos(phi)*R);

    if (isRotated_){
        particleCoordinates = Core::elevationRotate(particleCoordinates, elevation_);
        particleCoordinates = Core::azimuthRotate(particleCoordinates, azimuth_);
    }

    if (isShifted_){
        particleCoordinates = particleCoordinates + baseVector_;
    }

    return particleCoordinates;
}

/**
 * Generates a set of random ions in the cylindrical ion start zone
 *
 * @param numIons number of ions
 * @param charge charge of the generated ions
 * @param timeOfBirthRange ions are generated with times of birth uniformly distributed in this range
 * @return vector of random ions in cylindrical start zone
 */
std::vector<std::unique_ptr<BTree::Particle>> ParticleSimulation::CylinderStartZone::getRandomParticlesInStartZone(
        int numIons, double charge, double timeOfBirthRange) {

    Core::RndDistPtr rnd_tob = Core::globalRandomGenerator->getUniformDistribution(0,timeOfBirthRange);

    std::vector<std::unique_ptr<BTree::Particle>> result;
    //Core::Particle* result = new Core::Particle[numIons];
    for (int i=0; i<numIons; i++){
        std::unique_ptr<BTree::Particle> newIon = std::make_unique<BTree::Particle>(
                getRandomParticlePosition(),
                charge);

        newIon -> setTimeOfBirth(rnd_tob->rndValue());
        result.push_back(std::move(newIon));
    }
    return result;
}


