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
 * @param particles group / cloud of charged particles to be integrated
 * @param accelerationFunction a function to calculate the acceleration to the individual particles
 * @param timestepWriteFunction  a function to export data from the simulation
 * @param otherActionsFunction  a function to perform arbitrary other actions in every time step of the simulation
 * @param collisionModel a collision model, modeling the interaction between charged particles and background gas
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
    for (int i=0; i<particles.size(); ++i){
        addParticle(particles[i]);
    }
}

/**
 * Adds a particle to the velocity integrator
 * @param particle the particle to add
 */
void ParticleSimulation::VelocityIntegrator::addParticle(BTree::Particle *particle){
    particles_.push_back(particle);
    ++nParticles_;
}

/**
 * Run the integrator for a number of time steps and terminates the simulation.
 * @param nTimesteps number of time steps to run the integration
 * @param dt time step length
 */
void ParticleSimulation::VelocityIntegrator::run(int nTimesteps, double dt) {
    for (int step=0; step< nTimesteps; step++){
        runSingleStep(dt);
    }
    this->finalizeSimulation();
}

/**
 * Run the integrator for a single time step
 * @param dt time step length
 */
void ParticleSimulation::VelocityIntegrator::runSingleStep(double dt) {
    for (int i=0; i<nParticles_; i++){
        if (particles_[i]->isActive() == true){
            Core::Vector velocity = velocityFunction_(particles_[i],i,time_,timestep_);
            particles_[i]->setVelocity(velocity);
            particles_[i]->setLocation(particles_[i]->getLocation() + particles_[i]->getVelocity() * dt);
        }
    }

    timestepWriteFunction_(particles_,time_,timestep_,false);
    timestep_++;
    time_ = time_ + dt;

}

/**
 * Finalizes the verlet integration run (should be called after the last time step).
 */
void ParticleSimulation::VelocityIntegrator::finalizeSimulation(){
    timestepWriteFunction_(particles_,time_,timestep_,true);
}