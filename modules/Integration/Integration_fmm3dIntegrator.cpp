/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2022 - Physical and Theoretical Chemistry /
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

#include "Integration_fmm3dIntegrator.hpp"

Integration::FMMVerletIntegrator::FMMVerletIntegrator(
        const std::vector<Core::Particle*>& particles,
        accelerationFctType accelerationFunction,
        timestepWriteFctType timestepWriteFunction,
        otherActionsFctType otherActionsFunction,
        AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel):
    AbstractTimeIntegrator(particles, ionStartMonitoringFunction),
    collisionModel_(collisionModel),
    accelerationFunction_(std::move(accelerationFunction)),
    timestepWriteFunction_(std::move(timestepWriteFunction)),
    otherActionsFunction_(std::move(otherActionsFunction))
{}

Integration::FMMVerletIntegrator::FMMVerletIntegrator(
        Integration::accelerationFctType accelerationFunction,
        Integration::timestepWriteFctType timestepWriteFunction,
        Integration::otherActionsFctType otherActionsFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel) :
    AbstractTimeIntegrator(ionStartMonitoringFunction),
    collisionModel_(collisionModel),
    accelerationFunction_(std::move(accelerationFunction)),
    timestepWriteFunction_(std::move(timestepWriteFunction)),
    otherActionsFunction_(std::move(otherActionsFunction))
{}

void Integration::FMMVerletIntegrator::addParticle(Core::Particle* particle) {
    particles_.push_back(particle);
    solver_.insertParticle(*particle, nParticles_);
    ++nParticles_;
}

void Integration::FMMVerletIntegrator::run(unsigned int nTimesteps, double dt) {
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

void Integration::FMMVerletIntegrator::runSingleStep(double dt) {
    bearParticles_(time_);

    // first: Calculate charge distribution with particle positions and charges at the beginning of the time step
    solver_.computeChargeDistribution();

    // then: perform particle update / velocity verlet time integration
    std::size_t i;
    #pragma omp parallel \
            default(none) shared(newPos_, a_tdt_, a_t_, dt, particles_) \
            private(i) //firstprivate(MyNod)
    {

        #pragma omp for
        for (i=0; i<nParticles_; i++){

            if (particles_[i]->isActive()){

                if (collisionModel_ != nullptr) {
                    collisionModel_->updateModelParameters(*(particles_[i]));
                }

                newPos_[i] = particles_[i]->getLocation() + particles_[i]->getVelocity() * dt + a_t_[i]*(1.0/2.0*dt*dt);
                a_tdt_[i] = accelerationFunction_(particles_[i], i, solver_, time_, timestep_);
                //acceleration changes due to background interaction:

                if (collisionModel_ != nullptr) {
                    collisionModel_->modifyAcceleration(a_tdt_[i], *(particles_[i]), dt);
                }

                particles_[i]->setVelocity( particles_[i]->getVelocity() + ((a_t_[i]+ a_tdt_[i])*1.0/2.0 *dt) );
                a_t_[i] = a_tdt_[i];

                //velocity changes due to background interaction:
                if (collisionModel_ != nullptr) {
                    collisionModel_->modifyVelocity(*(particles_[i]),dt);
                    collisionModel_->modifyPosition(newPos_[i], *(particles_[i]), dt);
                }

                if (otherActionsFunction_ != nullptr) {
                    otherActionsFunction_(newPos_[i], particles_[i], i, time_, timestep_);
                }

            }
        }
    }

    // Update time and timestep, if existting: Call write timestep function
    time_ = time_ + dt;
    timestep_++;
    if (timestepWriteFunction_ != nullptr) {
        timestepWriteFunction_(particles_, time_, timestep_, false);
    }
}

void Integration::FMMVerletIntegrator::finalizeSimulation() {
    if (timestepWriteFunction_ != nullptr){
        timestepWriteFunction_(particles_, time_, timestep_, true);
    }
}