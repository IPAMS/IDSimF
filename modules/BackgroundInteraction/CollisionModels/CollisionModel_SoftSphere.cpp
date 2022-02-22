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

#include "CollisionModel_SoftSphere.hpp"

void CollisionModel::SoftSphereModel::initializeModelParameters(BTree::Particle &ion) const {

}

void CollisionModel::SoftSphereModel::updateModelParameters(BTree::Particle &ion) const {

}

void CollisionModel::SoftSphereModel::modifyAcceleration(Core::Vector &acceleration, BTree::Particle &particle,
                                                         double dt) {

}

void CollisionModel::SoftSphereModel::modifyVelocity(BTree::Particle &particle, double dt) {

}

void CollisionModel::SoftSphereModel::modifyPosition(Core::Vector &position, BTree::Particle &particle, double dt) {

}