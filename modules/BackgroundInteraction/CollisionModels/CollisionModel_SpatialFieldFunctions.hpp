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
 CollisionModel_SpatialFieldFunctions.hpp

 Helper / convenience functions defining simple spatial fields
 (e.g. Constant functions returning a static value)

 ****************************/

#ifndef Collision_SpatialFieldFunctions_hpp
#define Collision_SpatialFieldFunctions_hpp

#include "Core_vector.hpp"
#include "PSim_interpolatedField.hpp"
#include "PSim_simionPotentialArray.hpp"
#include <functional>

namespace CollisionModel{
    
    std::function<double(const Core::Vector&)>getConstantScalarFunction(double constantValue);
    std::function<Core::Vector(const Core::Vector&)>getConstantVectorFunction(Core::Vector constantValue);
    std::function<double(const Core::Vector&)>getVariableScalarFunction(ParticleSimulation::InterpolatedField &interpolatedField, std::size_t fieldIndex = 0);
    std::function<Core::Vector(const Core::Vector&)>getVariableVectorFunction(ParticleSimulation::InterpolatedField &interpolatedField, std::size_t fieldIndex = 0);
    std::function<double(const Core::Vector&)>getVariableScalarFunction(ParticleSimulation::SimionPotentialArray &simionPA);

}

#endif /* Collision_SpatialFieldFunctions_hpp */
