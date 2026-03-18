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

#include "Integration_parallelVerletIntegrator.hpp"
#include <utility>
#include <algorithm>

Integration::ParallelVerletIntegrator::ParallelVerletIntegrator(
        const std::vector<Core::Particle *>& particles,
        Integration::accelerationFctSingleStepType accelerationFunction,
        Integration::postTimestepFctType postTimestepFunction,
        Integration::otherActionsFctType otherActionsFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel) :
        AbstractTimeIntegrator(particles, ionStartMonitoringFunction),
        collisionModel_(collisionModel),
        accelerationFunction_(std::move(accelerationFunction)),
        postTimestepWriteFunction_(std::move(postTimestepFunction)),
        otherActionsFunction_(std::move(otherActionsFunction))
{}

Integration::ParallelVerletIntegrator::ParallelVerletIntegrator(
        Integration::accelerationFctSingleStepType accelerationFunction,
        Integration::postTimestepFctType timestepWriteFunction,
        Integration::otherActionsFctType postTimestepFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel) :
        AbstractTimeIntegrator(ionStartMonitoringFunction),
        collisionModel_(collisionModel),
        accelerationFunction_(std::move(accelerationFunction)),
        postTimestepWriteFunction_(std::move(timestepWriteFunction)),
        otherActionsFunction_(std::move(postTimestepFunction))
{
    initInternalState_();
}


/**
 * Adds a particle to the verlet integrator (required if particles are generated in the course of the simulation
 * @param particle the particle to add to the verlet integration
 */
void Integration::ParallelVerletIntegrator::addParticle(Core::Particle *particle){
    particles_.push_back(particle);
    a_t_.emplace_back(Core::Vector(0,0,0));
    a_tdt_.emplace_back(Core::Vector(0,0,0));

    tree_.insertParticle(*particle, nParticles_);
    ++nParticles_;
}

void Integration::ParallelVerletIntegrator::bearParticles_(double time) {
    Integration::AbstractTimeIntegrator::bearParticles_(time);
    initInternalState_();
}

void Integration::ParallelVerletIntegrator::initInternalState_(){
    tree_.init();
}

/**
 * Runs the integration
 * @param nTimesteps number of time steps to run
 * @param dt time step length
 */
void Integration::ParallelVerletIntegrator::run(unsigned int nTimesteps, double dt) {

    // run init:
    this->runState_ = RUNNING;
    bearParticles_(0.0);

    if (postTimestepWriteFunction_ !=nullptr) {
        postTimestepWriteFunction_(this, particles_, time_, timestep_, false);
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

/**
 * Runs a single step of the integration
 * @param dt time step length
 */
void Integration::ParallelVerletIntegrator::runSingleStep(double dt){

    bearParticles_(time_);

    int ver=0;

    if (collisionModel_ !=nullptr){
        collisionModel_->updateModelTimestepParameters(timestep_, time_);
    }
    std::size_t i;
    #pragma omp parallel \
            default(none) shared(a_t_, dt, particles_) \
            private(i)
    {
        #pragma omp for schedule(dynamic, 40)
        for (i=0; i<nParticles_; ++i){

            if (particles_[i]->isActive()){

                particles_[i]->setLocation(
                    particles_[i]->getLocation() + particles_[i]->getVelocity() * dt + a_t_[i]*(1.0/2.0*dt*dt));

                if (collisionModel_ != nullptr) {
                    collisionModel_->updateModelParticleParameters(*(particles_[i]));
                }
            }
        }
    }

    // Update particle positions and serialized structure of the tree
    for (std::size_t i=0; i<nParticles_; i++){
        if (particles_[i]->isActive()){
            tree_.updateParticleLocation(i, &ver);
        }
    }

    // Update serialized tree structure:
    tree_.updateNodes(ver);

    // Now calculate acceleration and velocity and perform particle modification in a parallel way
    #pragma omp parallel \
    default(none) shared(a_tdt_, a_t_, dt, particles_) \
    private(i)
    {
        #pragma omp for schedule(dynamic, 40)
        for (i=0; i<nParticles_; ++i){

            if (particles_[i]->isActive()){
                a_tdt_[i] = accelerationFunction_(particles_[i], i, tree_, time_, timestep_);
                //acceleration changes due to background interaction:

                if (collisionModel_ != nullptr) {
                    collisionModel_->modifyAcceleration(a_tdt_[i], *(particles_[i]), dt);
                }

                particles_[i]->setVelocity( particles_[i]->getVelocity() + ((a_t_[i]+ a_tdt_[i])*1.0/2.0 *dt) );
                a_t_[i] = a_tdt_[i];

                //velocity changes due to background interaction:
                if (collisionModel_ != nullptr) {
                    collisionModel_->modifyVelocity(*(particles_[i]),dt);
                    collisionModel_->modifyPosition(*(particles_[i]), dt);
                }

                // other actions:
                if (otherActionsFunction_ != nullptr) {
                    otherActionsFunction_(particles_[i], i, time_, timestep_);
                }
            }
        }
    }

    time_ = time_ + dt;
    timestep_++;
    if (postTimestepWriteFunction_ != nullptr) {
        postTimestepWriteFunction_(this, particles_, time_, timestep_, false);
    }
}

/**
 * Finalizes the verlet integration run (should be called after the last time step).
 */
void Integration::ParallelVerletIntegrator::finalizeSimulation(){
    if (postTimestepWriteFunction_ != nullptr){
        postTimestepWriteFunction_(this, particles_, time_, timestep_, true);
    }
}

/**
 * Sets the theta value (multipole acceptance criterion) for the tree used by the integrator
 */
void Integration::ParallelVerletIntegrator::setTheta(double newTheta) {
    tree_.getRoot()->setTheta(newTheta);
}



