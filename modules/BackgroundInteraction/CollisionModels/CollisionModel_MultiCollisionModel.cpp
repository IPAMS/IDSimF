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

#include "CollisionModel_MultiCollisionModel.hpp"

/**
 * Constructs a multi collision model from a vector of collision models.
 * Note that the multi collision model does *not* check in any way if the combination
 * of collision models is physically reasonable at all.
 *
 * @param models a vector of collision models to be combined
 */
CollisionModel::MultiCollisionModel::MultiCollisionModel(
        std::vector<std::unique_ptr<CollisionModel::AbstractCollisionModel>> models) :
        models_(std::move(models)){}

/**
 * Calls updateModelParticleParameters for all combined sub models
 */
void CollisionModel::MultiCollisionModel::updateModelParticleParameters(Core::Particle &ion) const{
    for(const auto &model: models_){
        model->updateModelParticleParameters(ion);
    }
}

/**
 * Calls initializeModelParticleParameters for all combined sub models
 */
void CollisionModel::MultiCollisionModel::initializeModelParticleParameters(Core::Particle &ion) const{
    for(const auto &model: models_){
        model->initializeModelParticleParameters(ion);
    }
}

/**
 * Calls updateModelTimestepParameters for all combined sub models
 */
void CollisionModel::MultiCollisionModel::updateModelTimestepParameters(int timestep, double time) {
    for(const auto &model: models_){
        model->updateModelTimestepParameters(timestep, time);
    }
}

/**
 * Calls modifyVelocity for all combined sub models
 */
void CollisionModel::MultiCollisionModel::modifyVelocity(Core::Particle &ion, double dt){
    for(const auto &model: models_){
        model->modifyVelocity(ion,dt);
    }
}

/**
 * Calls modifyPosition for all combined sub models
 */
void CollisionModel::MultiCollisionModel::modifyPosition(Core::Vector &position, Core::Particle &ion, double dt){
    for(const auto &model: models_){
        model->modifyPosition(position,ion,dt);
    }
}


/**
 * Calls modifyAcceleration for all combined sub models
 */
void CollisionModel::MultiCollisionModel::modifyAcceleration(Core::Vector &acceleration, Core::Particle &ion,
                                                             double dt){
    for(const auto &model: models_){
        model->modifyAcceleration(acceleration,ion,dt);
    }
}
