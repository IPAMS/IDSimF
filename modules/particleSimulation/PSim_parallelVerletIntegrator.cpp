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

#include "PSim_parallelVerletIntegrator.hpp"
#include <utility>
#include <iostream>
#include <algorithm>
#include <iostream>

ParticleSimulation::ParallelVerletIntegrator::ParallelVerletIntegrator(std::vector<BTree::Particle *> particles,
                               ParticleSimulation::ParallelVerletIntegrator::accelerationFctType accelerationFunction,
                               ParticleSimulation::ParallelVerletIntegrator::timestepWriteFctType timestepWriteFunction,
                               ParticleSimulation::ParallelVerletIntegrator::otherActionsFctType otherActionsFunction,
                               CollisionModel::AbstractCollisionModel &collisionModel):
AbstractTimeIntegrator(particles),
accelerationFunction_(std::move(accelerationFunction)),
timestepWriteFunction_(std::move(timestepWriteFunction)),
otherActionsFunction_(std::move(otherActionsFunction)),
collisionModel_(&collisionModel)
{}

ParticleSimulation::ParallelVerletIntegrator::ParallelVerletIntegrator(
            ParticleSimulation::ParallelVerletIntegrator::accelerationFctType accelerationFunction,
            ParticleSimulation::ParallelVerletIntegrator::timestepWriteFctType timestepWriteFunction,
            ParticleSimulation::ParallelVerletIntegrator::otherActionsFctType otherActionsFunction,
            CollisionModel::AbstractCollisionModel &collisionModel) :
        accelerationFunction_(std::move(accelerationFunction)),
        timestepWriteFunction_(std::move(timestepWriteFunction)),
        otherActionsFunction_(std::move(otherActionsFunction)),
        collisionModel_(&collisionModel)
{
    initInternalState_();
}

/**
 * Adds a particle to the verlet integrator (required if particles are generated in the course of the simulation
 * @param particle the particle to add to the verlet integration
 */
void ParticleSimulation::ParallelVerletIntegrator::addParticle(BTree::Particle *particle){
    particles_.push_back(particle);
    newPos_.push_back(Core::Vector(0,0,0));
    a_t_.push_back(Core::Vector(0,0,0));
    a_tdt_.push_back(Core::Vector(0,0,0));

    tree_.insertParticle(*particle, nParticles_);
    ++nParticles_;
}

void ParticleSimulation::ParallelVerletIntegrator::bearParticles_(double time) {
    ParticleSimulation::AbstractTimeIntegrator::bearParticles_(time);
    initInternalState_();
}

void ParticleSimulation::ParallelVerletIntegrator::initInternalState_(){
    numberOfNodes_ = tree_.init();
}

void ParticleSimulation::ParallelVerletIntegrator::run(int nTimesteps, double dt) {

    //write initial state to the trajectory:
    timestepWriteFunction_(particles_, tree_, time_, timestep_, false);

    for (int step=0; step< nTimesteps; step++){
        runSingleStep(dt);
    }
    this->terminateSimulation();
}

void ParticleSimulation::ParallelVerletIntegrator::runSingleStep(double dt){

    //std::cout << "runSingleStep "<<dt<<" "<<time_<<std::endl;
    //first: Generate new particles if necessary
    bearParticles_(time_);

    int ver=0;

    // Feld, das für die Kraftberechnung benutzt wird
    // Vector of tree nodes which is used for the serialized, non recursive force calculation
    // (
    //std::vector<BTree::ParallelNode*> MyNod(numberOfNodes_, nullptr);

    //
    int i;
    #pragma omp parallel \
            default(none) shared(newPos_, a_tdt_, a_t_, dt, particles_) \
            private(i) //firstprivate(MyNod)
    {

        #pragma omp for
        for (i=0; i<nParticles_; i++){

            if (particles_[i]->isActive() == true){

                collisionModel_->updateModelParameters(*(particles_[i]));

                newPos_[i] = particles_[i]->getLocation() + particles_[i]->getVelocity() * dt + a_t_[i]*(1.0/2.0*dt*dt);
                a_tdt_[i] = accelerationFunction_(particles_[i],i,tree_,time_,timestep_);
                //acceleration changes due to background interaction:
                collisionModel_->modifyAcceleration(a_tdt_[i],*(particles_[i]),dt);

                particles_[i]->setVelocity( particles_[i]->getVelocity() + ((a_t_[i]+ a_tdt_[i])*1.0/2.0 *dt) );
                a_t_[i] = a_tdt_[i];

                //velocity changes due to background interaction:
                collisionModel_->modifyVelocity(*(particles_[i]),dt);
            }
        }
    }

    // First find all new positions, then perform otherActions then update tree.
    // This ensures that all new particle positions are found with the state from
    // last time step. No particle positions are found with a partly updated tree.
    // Additionally, the update step is not (yet) parallel.
    for (int i=0; i<nParticles_; i++){
        if (particles_[i]->isActive() == true){
            //position changes due to background interaction:
            collisionModel_->modifyPosition(newPos_[i],*(particles_[i]),dt);

            otherActionsFunction_(newPos_[i],particles_[i],i,tree_,time_,timestep_);
            tree_.updateParticleLocation(i, newPos_[i], &ver);
        }
    }

    // Update serialized tree structure:
    numberOfNodes_ = tree_.updateNodes(ver);
    time_ = time_ + dt;
    timestep_++;
    timestepWriteFunction_(particles_,tree_,time_,timestep_,false);
}



void ParticleSimulation::ParallelVerletIntegrator::terminateSimulation(){
    timestepWriteFunction_(particles_,tree_,time_,timestep_,true);
}


