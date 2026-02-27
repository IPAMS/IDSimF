/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2024 - Physical and Theoretical Chemistry /
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

#include "CollisionModel_AbstractMDForceField.hpp"

void CollisionModel::AbstractMDForceField::calculateTorque(std::vector<Core::Vector>& positions, 
                    Core::Vector& forceMolecules, 
                    std::vector<Core::Vector>& torqueMolecules) {
    
    torqueMolecules[0].x(torqueMolecules[0].x() + (positions[0].y()*forceMolecules.z() - positions[0].z()*forceMolecules.y())); 
    torqueMolecules[0].y(torqueMolecules[0].y() + (positions[0].z()*forceMolecules.x() - positions[0].x()*forceMolecules.z()));     
    torqueMolecules[0].z(torqueMolecules[0].z() + (positions[0].x()*forceMolecules.y() - positions[0].y()*forceMolecules.x()));  
    torqueMolecules[1].x(torqueMolecules[0].x() - (positions[1].y()*forceMolecules.z() - positions[1].z()*forceMolecules.y())); 
    torqueMolecules[1].y(torqueMolecules[0].y() - (positions[1].z()*forceMolecules.x() - positions[1].x()*forceMolecules.z()));     
    torqueMolecules[1].z(torqueMolecules[0].z() - (positions[1].x()*forceMolecules.y() - positions[1].y()*forceMolecules.x()));               

}