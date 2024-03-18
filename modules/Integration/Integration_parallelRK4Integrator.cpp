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

#include "Integration_parallelRK4Integrator.hpp"
#include <utility>
#include <algorithm>

Integration::ParallelRK4Integrator::ParallelRK4Integrator(
        const std::vector<Core::Particle *>& particles,
        Integration::accelerationFctType accelerationFunction,
        Integration::accelerationFctSpaceChargeType spaceChargeAccelerationFunction,
        Integration::timestepWriteFctType timestepWriteFunction,
        Integration::otherActionsFctType otherActionsFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel) :
    AbstractTimeIntegrator(particles, ionStartMonitoringFunction),
    collisionModel_(collisionModel),
    accelerationFunction_(std::move(accelerationFunction)),
    timestepWriteFunction_(std::move(timestepWriteFunction)),
    otherActionsFunction_(std::move(otherActionsFunction))
{}

Integration::ParallelRK4Integrator::ParallelRK4Integrator(
        Integration::accelerationFctType accelerationFunction,
        Integration::accelerationFctSpaceChargeType spaceChargeAccelerationFunction,
        Integration::timestepWriteFctType timestepWriteFunction,
        Integration::otherActionsFctType otherActionsFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel) :
    AbstractTimeIntegrator(ionStartMonitoringFunction),
    collisionModel_(collisionModel),
    accelerationFunction_(std::move(accelerationFunction)),
    timestepWriteFunction_(std::move(timestepWriteFunction)),
    otherActionsFunction_(std::move(otherActionsFunction))
{
    initInternalState_();
}


/**
 * Adds a particle to the verlet integrator (required if particles are generated in the course of the simulation
 * @param particle the particle to add to the verlet integration
 */
void Integration::ParallelRK4Integrator::addParticle(Core::Particle *particle){
    particles_.push_back(particle);
    newPos_.emplace_back(Core::Vector(0,0,0));
    spaceChargeAcceleration_.emplace_back(Core::Vector(0,0,0));

    tree_.insertParticle(*particle, nParticles_);
    ++nParticles_;
}

void Integration::ParallelRK4Integrator::bearParticles_(double time) {
    Integration::AbstractTimeIntegrator::bearParticles_(time);
    initInternalState_();
}

void Integration::ParallelRK4Integrator::initInternalState_(){
    tree_.init();
}

void Integration::ParallelRK4Integrator::run(unsigned int nTimesteps, double dt) {

    // run init:
    this->runState_ = RUNNING;
    bearParticles_(0.0);

    if (timestepWriteFunction_ !=nullptr) {
        timestepWriteFunction_(particles_, time_, timestep_, false);
    }

    // run:
    for (unsigned int step=0; step< nTimesteps; step++){
        runSingleStep(dt);
        if (this->runState_ == IN_TERMINATION){
            break;
        }
    }
    this->finalizeSimulation();
    this->runState_ = STOPPED;
}

void Integration::ParallelRK4Integrator::runSingleStep(double dt){

    //first: Generate new particles if necessary
    bearParticles_(time_);

    if (collisionModel_ !=nullptr){
        collisionModel_->updateModelTimestepParameters(timestep_, time_);
    }
    std::size_t i;
    #pragma omp parallel \
            default(none) shared(newPos_, spaceChargeAcceleration_, dt, particles_) \
            private(i)
    {

        #pragma omp for schedule(dynamic, 40)
        for (i=0; i<nParticles_; i++){

            if (particles_[i]->isActive()){
                if (collisionModel_ != nullptr) {
                    collisionModel_->updateModelParticleParameters(*(particles_[i]));
                }


                // calculation of acceleration



                //acceleration changes due to background interaction:

                if (collisionModel_ != nullptr) {
                    collisionModel_->modifyAcceleration(a_tdt_[i], *(particles_[i]), dt);
                }


                //velocity changes due to background interaction:
                if (collisionModel_ != nullptr) {
                    //std::cout << "before:" << particles_[i]->getVelocity() << std::endl;
                    collisionModel_->modifyVelocity(*(particles_[i]),dt);
                    //std::cout << "after:" << particles_[i]->getVelocity() << std::endl;
                }
            }
        }
    }

    // First find all new positions, then perform otherActions then update tree.
    // This ensures that all new particle positions are found with the state from
    // last time step. No particle positions are found with a partly updated tree.
    // Additionally, the update step is not (yet) parallel.
    for (std::size_t i=0; i<nParticles_; i++){
        if (particles_[i]->isActive()){
            //position changes due to background interaction:
            if (collisionModel_ != nullptr) {
                collisionModel_->modifyPosition(newPos_[i], *(particles_[i]), dt);
            }

            if (otherActionsFunction_ != nullptr) {
                otherActionsFunction_(newPos_[i], particles_[i], i, time_, timestep_);
            }
            tree_.updateParticleLocation(i, newPos_[i], &ver);
        }
    }

    // Update serialized tree structure:
    tree_.updateNodes(ver);
    time_ = time_ + dt;
    timestep_++;
    if (timestepWriteFunction_ != nullptr) {
        timestepWriteFunction_(particles_, time_, timestep_, false);
    }
}

/**
 * Finalizes the RK4 integration run (should be called after the last time step).
 */
void Integration::ParallelRK4Integrator::finalizeSimulation(){
    if (timestepWriteFunction_ != nullptr){
        timestepWriteFunction_(particles_, time_, timestep_, true);
    }
}


