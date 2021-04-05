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

#include "PSim_velocityIntegrator.hpp"
#include "BTree_particle.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include <utility>

/**
 * Creates a velocity integrator
 *
 * @param particles group / cloud of charged particles to be integrated
 * @param accelerationFunction A function to calculate the acceleration to the individual particles
 * @param timestepWriteFunction A function to export data from the simulation (can be nullptr to flag non usage)
 * @param otherActionsFunction A function to perform arbitrary other actions in every time step of the simulation
 *  (can be nullptr to flag non usage)

 */
ParticleSimulation::VelocityIntegrator::VelocityIntegrator(
        std::vector<BTree::Particle *> particles,
        ParticleSimulation::VelocityIntegrator::velocityFctType velocityFunction,
        ParticleSimulation::VelocityIntegrator::timestepWriteFctType timestepWriteFunction,
        ParticleSimulation::VelocityIntegrator::otherActionsFctType otherActionsFunction):
velocityFunction_(std::move(velocityFunction)),
timestepWriteFunction_(std::move(timestepWriteFunction)),
otherActionsFunction_(std::move(otherActionsFunction))
{
    //init velocities and accelerations:
    for (const auto& particle : particles){
        addParticle_(particle); //Do not call virtual functions in constructor, call private non virtual function
    }
}

/**
 * Adds a particle to the velocity integrator
 * @param particle the particle to add
 */
void ParticleSimulation::VelocityIntegrator::addParticle(BTree::Particle *particle){
    addParticle_(particle);
}

void ParticleSimulation::VelocityIntegrator::addParticle_(BTree::Particle *particle){
    particles_.push_back(particle);
    ++nParticles_;
}

/**
 * Run the integrator for a number of time steps and terminates the simulation.
 * @param nTimesteps number of time steps to run the integration
 * @param dt time step length
 */
void ParticleSimulation::VelocityIntegrator::run(int nTimesteps, double dt) {
    this->runState_ = RUNNING;
    for (int step=0; step< nTimesteps; step++){
        runSingleStep(dt);
        if (this->runState_ == IN_TERMINATION){
            break;
        }
    }
    this->finalizeSimulation();
    this->runState_ = STOPPED;
}

/**
 * Run the integrator for a single time step
 * @param dt time step length
 */
void ParticleSimulation::VelocityIntegrator::runSingleStep(double dt) {
    for (int i=0; i<nParticles_; i++){
        if (particles_[i]->isActive()){
            Core::Vector velocity = velocityFunction_(particles_[i],i,time_,timestep_);
            particles_[i]->setVelocity(velocity);
            Core::Vector newPos = particles_[i]->getLocation() + particles_[i]->getVelocity() * dt;
            particles_[i]->setLocation(newPos);
            if (otherActionsFunction_ !=nullptr) {
                otherActionsFunction_(newPos, particles_[i], i, time_, timestep_);
            }
        }
    }

    if (timestepWriteFunction_ != nullptr) {
        timestepWriteFunction_(particles_, time_, timestep_, false);
    }
    timestep_++;
    time_ = time_ + dt;

}

/**
 * Finalizes the verlet integration run (should be called after the last time step).
 */
void ParticleSimulation::VelocityIntegrator::finalizeSimulation(){
    if (timestepWriteFunction_ != nullptr) {
        timestepWriteFunction_(particles_, time_, timestep_, true);
    }
}