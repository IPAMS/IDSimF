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

#include "CollisionModel_SpatialFieldFunctions.hpp"

/**
 * Generates a constant scalar function
 *
 * @param constantValue the constant value which is returned by the function
 * @return spatial scalar function always returning the given constant value
 */
std::function<double(const Core::Vector&)> CollisionModel::getConstantDoubleFunction(double constantValue){
    
    return [=](const Core::Vector&)->double{ return constantValue;};
}

/**
 * Generates a constant vector valued function
 *
 * @param constantValue the constant vector which is returned by the function
 * @return a vector valued function always returning the given constant vector
 */
std::function<Core::Vector(const Core::Vector&)> CollisionModel::getConstantVectorFunction(Core::Vector constantValue){
    
    return [=](const Core::Vector&)->Core::Vector{ return constantValue;};
}