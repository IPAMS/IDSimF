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

#include "RS_ReactiveParticle.hpp"


RS::ReactiveParticle::ReactiveParticle(Substance* species){
    this->setSpecies(species);
}

RS::ReactiveParticle::ReactiveParticle(Substance* species, Core::Vector location):
    Core::Particle(location, species->charge())
{
    this->setSpecies(species);
}

RS::ReactiveParticle::ReactiveParticle(Substance* species, Core::Vector location, double charge):
    Core::Particle(location,charge)
{
    this->setSpecies(species);
}

void RS::ReactiveParticle::setSpecies(RS::Substance* species){
    species_ = species;
    this->updateParticleParametersFromSpecies_();
}

RS::Substance* RS::ReactiveParticle::getSpecies() const{
    return species_;
}

void RS::ReactiveParticle::updateParticleParametersFromSpecies_() {
    RS::Substance* species = this->getSpecies();
    this->setMobility(species->lowFieldMobility());
    this->setMassAMU(species->mass());
    this->setDiameter(species->collisionDiameter());
    this->setChargeElementary(species->charge());
}

std::ostream& operator<<(std::ostream& os, const RS::ReactiveParticle& particle)
{
    os <<" | "<< particle.species_->name();
    return os;
}