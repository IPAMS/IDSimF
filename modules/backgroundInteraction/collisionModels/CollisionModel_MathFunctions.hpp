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
 CollisionModel_MathFunctions.hpp

 Additional mathematical support functions used in collision models

 ****************************/

#ifndef Collision_MathFunctions_hpp
#define Collision_MathFunctions_hpp

#include "Core_vector.hpp"
#include "Core_constants.hpp"
#include <cmath>

namespace CollisionModel {

    double degToRad(double phi);
    double radToDeg(double phi);

    Core::Vector cartesianToPolar(Core::Vector vec);
    Core::Vector elevationRotate(Core::Vector vec, double angle);
    Core::Vector azimuthRotate(Core::Vector vec, double angle);

    Core::Vector sphereRand(double r);
}

#endif /* Collision_MathFunctions_hpp */
