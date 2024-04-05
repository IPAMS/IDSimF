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

#include "Integration_verletIntegrator.hpp"
#include "Core_particle.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include <utility>
#include <iostream>

/**
 * Creates a verlet integrator
 * @param particles group / cloud of charged particles to be integrated
 * @param accelerationFunction a function to calculate the acceleration to the individual particles
 * @param timestepWriteFunction  a function to export data from the simulation
 * @param otherActionsFunction  a function to perform arbitrary other actions in every time step of the simulation
 * @param collisionModel a collision model, modeling the interaction between charged particles and background gas
 */
Integration::VerletIntegrator::VerletIntegrator(
        std::vector<Core::Particle *> particles,
        Integration::accelerationFctSingleStepType accelerationFunction,
        Integration::postTimestepFctType timestepWriteFunction,
        Integration::otherActionsFctType otherActionsFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel) :
        Integration::AbstractTimeIntegrator(particles, ionStartMonitoringFunction),
        collisionModel_(collisionModel),
        accelerationFunction_(std::move(accelerationFunction)),
        postTimestepFunction_(std::move(timestepWriteFunction)),
        otherActionsFunction_(std::move(otherActionsFunction))
{}

/**
 * Creates an empty verlet integrator (without simulated particles)
 *
 * @param accelerationFunction a function to calculate the acceleration to the individual particles
 * @param timestepWriteFunction  a function to export data from the simulation
 * @param otherActionsFunction  a function to perform arbitrary other actions in every time step of the simulation
 * @param collisionModel a collision model, modeling the interaction between charged particles and background gas
 */
Integration::VerletIntegrator::VerletIntegrator(
        Integration::accelerationFctSingleStepType accelerationFunction,
        Integration::postTimestepFctType timestepWriteFunction,
        Integration::otherActionsFctType otherActionsFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel) :
        AbstractTimeIntegrator(ionStartMonitoringFunction),
        collisionModel_(collisionModel),
        accelerationFunction_(std::move(accelerationFunction)),
        postTimestepFunction_(std::move(timestepWriteFunction)),
        otherActionsFunction_(std::move(otherActionsFunction))
{}

/**
 * Adds a particle to the verlet integrator (required if particles are generated in the course of the simulation)
 * @param particle the particle to add to the verlet integration
 * @param extIndex an external numerical index / key for the particle to add
 */
void Integration::VerletIntegrator::addParticle(Core::Particle *particle){
    particles_.push_back(particle);
    newPos_.push_back(Core::Vector(0,0,0));
    a_t_.push_back(Core::Vector(0,0,0));
    a_tdt_.push_back(Core::Vector(0,0,0));

    tree_.insertParticle(*particle, nParticles_);
    ++nParticles_;
}

/**
 * Run the verlet integrator for a number of time steps and terminates the simulation.
 * @param nTimesteps number of time steps to run the verlet integration
 * @param dt time step length
 */
void Integration::VerletIntegrator::run(unsigned int nTimesteps, double dt) {

    this->runState_ = RUNNING;
    bearParticles_(0.0);
    if (postTimestepFunction_ !=nullptr) {
        postTimestepFunction_(this, particles_, time_, timestep_, false);
    }

    for (unsigned int step=0; step< nTimesteps; ++step){
        runSingleStep(dt);
        if (this->runState_ == IN_TERMINATION){
            break;
        }
    }
    this->finalizeSimulation();
    this->runState_ = STOPPED;
}

/**
 * Run the verlet integrator for a single time step
 * @param dt time step length
 */
void Integration::VerletIntegrator::runSingleStep(double dt) {

    bearParticles_(time_);

    if (collisionModel_ !=nullptr){
        collisionModel_->updateModelTimestepParameters(timestep_, time_);
    }

    for (std::size_t i=0; i<nParticles_; ++i){
        if (particles_[i]->isActive() == true){

            if (collisionModel_ !=nullptr) {
                collisionModel_->updateModelParticleParameters(*(particles_[i]));
            }

            newPos_[i] = particles_[i]->getLocation() + particles_[i]->getVelocity() * dt + a_t_[i]*(1.0/2.0*dt*dt);
            a_tdt_[i] = accelerationFunction_(particles_[i],i,tree_,time_,timestep_);
            //acceleration changes due to background interaction:
            if (collisionModel_ !=nullptr) {
                collisionModel_->modifyAcceleration(a_tdt_[i], *(particles_[i]), dt);
            }

            particles_[i]->setVelocity( particles_[i]->getVelocity() + ((a_t_[i]+ a_tdt_[i])*1.0/2.0 *dt) );
            a_t_[i] = a_tdt_[i];

            if (collisionModel_ !=nullptr){
                //velocity changes due to background interaction:
                collisionModel_->modifyVelocity(*(particles_[i]),dt);
                //position changes due to background interaction:
                collisionModel_->modifyPosition(newPos_[i],*(particles_[i]),dt);
            }
        }
    }

    // First find all new positions, then perform arbitrary otherActions and update tree.
    // This ensures that all new particle positions are found with the state from
    // last time step. No particle positions are found with a partly updated tree.
    for (std::size_t i=0; i<nParticles_; ++i){
        if (particles_[i]->isActive() == true){

            if (otherActionsFunction_ != nullptr) {
                otherActionsFunction_(newPos_[i], particles_[i], i, time_, timestep_);
            }
            tree_.updateParticleLocation(i,newPos_[i]);
        }
    }
    timestep_++;
    time_ = time_ + dt;
    if (postTimestepFunction_ != nullptr){
        postTimestepFunction_(this, particles_, time_, timestep_, false);
    }
}

/**
 * Finalizes the verlet integration run (should be called after the last time step).
 */
void Integration::VerletIntegrator::finalizeSimulation(){
    if (postTimestepFunction_ !=nullptr) {
        postTimestepFunction_(this, particles_, time_, timestep_, true);
    }
}

void Integration::VerletIntegrator::bearParticles_(double time) {
    bool particlesCreated = Integration::AbstractTimeIntegrator::bearParticles_(time);
    if (particlesCreated){
        tree_.computeChargeDistribution();
    }
}