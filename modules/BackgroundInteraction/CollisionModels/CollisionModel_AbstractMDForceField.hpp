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

 ------------
 CollisionModel_AbstractMDForceField.hpp

 Abstract base class for force fields for MD collision modelling

 ****************************/

#ifndef IDSIMF_COLLISIONMODEL_ABSTRACTMDFORCEFIELD_HPP
#define IDSIMF_COLLISIONMODEL_ABSTRACTMDFORCEFIELD_HPP

#include "CollisionModel_Molecule.hpp"

namespace CollisionModel {

    class AbstractMDForceField {

    public:
        virtual ~AbstractMDForceField() = default;

        virtual void calculateForceField(std::vector<CollisionModel::Molecule*>& moleculesPtr, std::vector<Core::Vector>& forceMolecules) = 0;
    };
}

#endif //IDSIMF_COLLISIONMODEL_ABSTRACTMDFORCEFIELD_HPP
