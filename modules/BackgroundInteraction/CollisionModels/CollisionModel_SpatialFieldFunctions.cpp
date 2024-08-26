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
std::function<double(const Core::Vector&)> CollisionModel::getConstantScalarFunction(double constantValue){
    
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

std::function<double(const Core::Vector&)> CollisionModel::getVariableScalarFunction(
    ParticleSimulation::InterpolatedField &interpolatedField, std::size_t fieldIndex) {

    return [&interpolatedField,fieldIndex](const Core::Vector& position)->double
    {
        return interpolatedField.getInterpolatedScalar(position.x(), position.y(), position.z(), fieldIndex);
    };
}

std::function<Core::Vector(const Core::Vector&)> CollisionModel::getVariableVectorFunction(
    ParticleSimulation::InterpolatedField &interpolatedField, std::size_t fieldIndex) {

    return [&interpolatedField,fieldIndex](const Core::Vector& position)->Core::Vector
    {
        return interpolatedField.getInterpolatedVector(position.x(), position.y(), position.z(), fieldIndex);
    };
}

std::function<double(const Core::Vector&)> CollisionModel::getVariableScalarFunction(
    ParticleSimulation::SimionPotentialArray &simionPA) {

    return [&simionPA](const Core::Vector& position)->double
    {
        return simionPA.getInterpolatedPotential(position.x(), position.y(), position.z());
    };
}

std::function<Core::Vector(const Core::Vector&)> CollisionModel::getVariableVectorFunction(
    ParticleSimulation::SimionPotentialArray& pa_x,
    ParticleSimulation::SimionPotentialArray& pa_y,
    ParticleSimulation::SimionPotentialArray& pa_z) {

    return [&pa_x, &pa_y, &pa_z](const Core::Vector& position)->Core::Vector
    {
        double x = position.x();
        double y = position.y();
        double z = position.z();
        return {
            pa_x.getInterpolatedPotential(x,y,z),
            pa_y.getInterpolatedPotential(x,y,z),
            pa_z.getInterpolatedPotential(x,y,z)
        };
    };
}
