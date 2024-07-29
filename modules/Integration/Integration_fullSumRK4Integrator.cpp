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

#include "Integration_fullSumRK4Integrator.hpp"
#include <utility>
#include <algorithm>

Integration::FullSumRK4Integrator::FullSumRK4Integrator(
        const std::vector<Core::Particle *>& particles,
        Integration::accelerationFctType accelerationFunction,
        Integration::accelerationFctSpaceChargeType spaceChargeAccelerationFunction,
        Integration::postTimestepFctType timestepWriteFunction,
        Integration::otherActionsFctType otherActionsFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel) :
        AbstractTimeIntegrator(particles, ionStartMonitoringFunction),
        collisionModel_(collisionModel),
        accelerationFunction_(std::move(accelerationFunction)),
        spaceChargeAccelerationFunction_(std::move(spaceChargeAccelerationFunction)),
        postTimestepFunction_(std::move(timestepWriteFunction)),
        otherActionsFunction_(std::move(otherActionsFunction))
{}

Integration::FullSumRK4Integrator::FullSumRK4Integrator(
        Integration::accelerationFctType accelerationFunction,
        Integration::accelerationFctSpaceChargeType spaceChargeAccelerationFunction,
        Integration::postTimestepFctType timestepWriteFunction,
        Integration::otherActionsFctType otherActionsFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel) :
        AbstractTimeIntegrator(ionStartMonitoringFunction),
        collisionModel_(collisionModel),
        accelerationFunction_(std::move(accelerationFunction)),
        spaceChargeAccelerationFunction_(std::move(spaceChargeAccelerationFunction)),
        postTimestepFunction_(std::move(timestepWriteFunction)),
        otherActionsFunction_(std::move(otherActionsFunction))
{}


/**
 * Adds a particle to the Full Sum RK4 integrator (required if particles are generated in the course of the simulation
 * @param particle the particle to add to the verlet integration
 */
void Integration::FullSumRK4Integrator::addParticle(Core::Particle *particle){
    particles_.push_back(particle);
    newPos_.emplace_back(Core::Vector(0,0,0));
    spaceChargeAcceleration_.emplace_back(Core::Vector(0,0,0));

    fullSumSolver_.insertParticle(*particle, nParticles_);
    ++nParticles_;
}

void Integration::FullSumRK4Integrator::bearParticles_(double time) {
    Integration::AbstractTimeIntegrator::bearParticles_(time);
}

void Integration::FullSumRK4Integrator::run(unsigned int nTimesteps, double dt) {

    // run init:
    this->runState_ = RUNNING;
    bearParticles_(0.0);

    if (postTimestepFunction_ !=nullptr) {
        postTimestepFunction_(this, particles_, time_, timestep_, false);
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

void Integration::FullSumRK4Integrator::runSingleStep(double dt){

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
                Core::Vector spaceChargeAcceleration = spaceChargeAccelerationFunction_(particles_[i], i, fullSumSolver_, time_, timestep_);

                // actual Runge Kutta 4 steps
                Core::Particle* particle = particles_[i];

                Core::Vector x_n = particle->getLocation();
                Core::Vector v_n = particle->getVelocity();

                Core::Vector k1_v = evaluateAccelerationFunction_(particle, x_n, v_n, spaceChargeAcceleration, time_, timestep_, dt)*dt;
                Core::Vector k1_x = v_n*dt;

                Core::Vector k2_v = evaluateAccelerationFunction_(
                        particle,
                        x_n + k1_x/2.0,
                        v_n + k1_v/2.0,
                        spaceChargeAcceleration,
                        time_+ dt/2.0,
                        timestep_,
                        dt) * dt;
                Core::Vector k2_x = (v_n + k1_v/2.0) * dt;

                Core::Vector k3_v = evaluateAccelerationFunction_(
                        particle,
                        x_n + k2_x/2.0,
                        v_n + k2_v/2.0,
                        spaceChargeAcceleration,
                        time_ + dt/2.0,
                        timestep_,
                        dt) * dt;
                Core::Vector k3_x = (v_n + k2_v/2.0) * dt;

                Core::Vector k4_v = evaluateAccelerationFunction_(
                        particle,
                        x_n + k3_x,
                        v_n + k3_v,
                        spaceChargeAcceleration,
                        time_ + dt,
                        timestep_,
                        dt) * dt;
                Core::Vector k4_x = (v_n + k3_v) * dt;

                Core::Vector v_new = v_n + 1.0/6.0 * (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v);
                newPos_[i] = x_n + 1.0/6.0 * (k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);

                particle->setVelocity(v_new);
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
            particles_[i]->setLocation(newPos_[i]);
        }
    }

    time_ = time_ + dt;
    timestep_++;
    if (postTimestepFunction_ != nullptr) {
        postTimestepFunction_(this, particles_, time_, timestep_, false);
    }
}

/**
 * Finalizes the RK4 integration run (should be called after the last time step).
 */
void Integration::FullSumRK4Integrator::finalizeSimulation(){
    if (postTimestepFunction_ != nullptr){
        postTimestepFunction_(this, particles_, time_, timestep_, true);
    }
}

Core::Vector Integration::FullSumRK4Integrator::evaluateAccelerationFunction_(Core::Particle* particle,
                                                                               Core::Vector position,
                                                                               Core::Vector velocity,
                                                                               Core::Vector spaceChargeAcceleration,
                                                                               double time,
                                                                               unsigned int timestep,
                                                                               double dt) {
    Core::Vector acceleration = accelerationFunction_(particle, position, velocity, time, timestep) + spaceChargeAcceleration;
    if (collisionModel_ != nullptr) {
        collisionModel_->modifyAcceleration(acceleration, *particle, dt);
    }
    return acceleration;
}

